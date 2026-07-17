"""Tests for aperture_snv.models.gbm — mechanics on synthetic data (no real data yet).

These validate encoding, training/prediction, calibration, save/load, and the PPM
weight seam. They do not assess real-world accuracy — that needs labeled fragments.
"""

import numpy as np

from aperture_snv.extract import FEATURE_COLS, FragmentFeatures
from aperture_snv.models.gbm import (
    MONOTONE_CONSTRAINTS,
    GBMModel,
    fragments_to_matrix,
)
from aperture_snv.score import sum_supporting_reads


def _frag(rng=None, *, signal: bool, supports_alt: bool = True) -> FragmentFeatures:
    """A fragment whose features carry (noisy) signal for the `signal` class.

    signal=True  -> high quality, well aligned, concordant (ctDNA-like)
    signal=False -> low quality, mismatched, 3'-end, discordant (artifact-like)
    """
    r = rng.normal if rng is not None else (lambda loc, scale: loc)
    if signal:
        return FragmentFeatures(
            chrom="chr1", pos=100, ref="C", alt="T", frag_id="f",
            supports_alt=supports_alt, read_base="T" if supports_alt else "C",
            vbq=int(r(37, 2)), mrbq=float(r(34, 2)), n_low_bq=max(0, int(r(1, 1))),
            mapq=int(r(60, 1)), edit_distance=max(0, int(r(0, 0.5))), pir=float(r(0.5, 0.1)),
            dist_3prime=int(r(65, 8)), concordance="concordant",
            fragment_length=int(r(166, 10)), softclip_len=max(0, int(r(0, 0.5))),
            is_proper_pair=True, n_mates_at_site=2,
        )
    return FragmentFeatures(
        chrom="chr1", pos=100, ref="C", alt="T", frag_id="f",
        supports_alt=supports_alt, read_base="T" if supports_alt else "C",
        vbq=int(r(22, 3)), mrbq=float(r(21, 3)), n_low_bq=max(0, int(r(8, 2))),
        mapq=int(r(35, 5)), edit_distance=max(0, int(r(4, 1))), pir=float(r(0.9, 0.05)),
        dist_3prime=int(r(6, 3)), concordance="discordant",
        fragment_length=int(r(120, 15)), softclip_len=max(0, int(r(6, 2))),
        is_proper_pair=False, n_mates_at_site=2,
    )


def _trained_model(seed=0):
    rng = np.random.RandomState(seed)
    pos = [_frag(rng, signal=True) for _ in range(300)]
    neg = [_frag(rng, signal=False) for _ in range(300)]
    return GBMModel.train(pos, neg), rng


class TestEncoding:
    def test_matrix_shape_and_order(self):
        x = fragments_to_matrix([_frag(signal=True), _frag(signal=False)])
        assert x.shape == (2, len(FEATURE_COLS))

    def test_empty(self):
        assert fragments_to_matrix([]).shape == (0, len(FEATURE_COLS))

    def test_concordance_encoded_numerically(self):
        x = fragments_to_matrix([_frag(signal=True)])  # concordant -> 0
        assert x[0, FEATURE_COLS.index("concordance")] == 0.0

    def test_monotone_constraints_align_with_features(self):
        assert len(MONOTONE_CONSTRAINTS) == len(FEATURE_COLS)
        assert MONOTONE_CONSTRAINTS[FEATURE_COLS.index("edit_distance")] == -1
        assert MONOTONE_CONSTRAINTS[FEATURE_COLS.index("vbq")] == 1


class TestTrainPredict:
    def test_separates_classes(self):
        model, rng = _trained_model()
        p_pos = model.predict_proba([_frag(rng, signal=True) for _ in range(50)])
        p_neg = model.predict_proba([_frag(rng, signal=False) for _ in range(50)])
        assert p_pos.mean() > p_neg.mean()

    def test_probabilities_in_unit_interval(self):
        model, rng = _trained_model()
        p = model.predict_proba([_frag(rng, signal=bool(i % 2)) for i in range(40)])
        assert p.min() >= 0.0 and p.max() <= 1.0

    def test_predict_empty(self):
        model, _ = _trained_model()
        assert model.predict_proba([]).shape == (0,)


class TestCalibration:
    def test_calibrated_probabilities_valid(self):
        model, rng = _trained_model()
        val = [_frag(rng, signal=bool(i % 2)) for i in range(200)]
        labels = [i % 2 for i in range(200)]
        model.calibrate(val, labels)
        p = model.predict_proba(val)
        assert p.min() >= 0.0 and p.max() <= 1.0
        assert model.calibrator is not None


class TestPersistence:
    def test_save_load_roundtrip(self, tmp_path):
        model, rng = _trained_model()
        frags = [_frag(rng, signal=True) for _ in range(20)]
        before = model.predict_proba(frags)
        path = tmp_path / "gbm.joblib"
        model.save(path)
        after = GBMModel.load(path).predict_proba(frags)
        np.testing.assert_allclose(before, after)


class TestPPMWeightSeam:
    def test_weights_feed_sum_supporting_reads(self):
        """predict_proba supplies the per-fragment weight for the PPM numerator."""
        model, rng = _trained_model()
        frags = [_frag(rng, signal=True, supports_alt=True) for _ in range(10)]
        frags += [_frag(rng, signal=True, supports_alt=False) for _ in range(5)]  # ref reads
        probs = {id(f): p for f, p in zip(frags, model.predict_proba(frags), strict=True)}
        weighted = sum_supporting_reads(frags, weight_fn=lambda f: float(probs[id(f)]))
        # only the 10 alt-supporting fragments contribute; each weight in [0, 1]
        assert 0.0 <= weighted <= 10.0
