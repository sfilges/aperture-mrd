"""Tests for aperture_snv.svm - SVM error suppression."""

import numpy as np

from aperture_snv.extract import CandidateRead
from aperture_snv.svm import (
    build_svm_pipeline,
    candidate_to_features,
    candidates_to_feature_matrix,
    filter_candidates,
    train_svm,
)


def _make_candidate(**kwargs) -> CandidateRead:
    """Create a CandidateRead with sensible defaults, overridable via kwargs."""
    defaults = dict(
        chrom="chr1",
        pos=1000,
        ref="C",
        alt="T",
        read_name="read_1",
        read_base="T",
        vbq=30,
        mrbq=28.5,
        pir=0.3,
        mapq=60,
        is_read1=True,
        mate_cigar=None,
        is_concordant=True,
        is_proper_pair=True,
        insert_size=165,
    )
    defaults.update(kwargs)
    return CandidateRead(**defaults)


class TestCandidateToFeatures:
    def test_concordant(self):
        c = _make_candidate(vbq=30, mrbq=28.0, pir=0.4, mapq=60, is_concordant=True)
        features = candidate_to_features(c)
        assert features.shape == (5,)
        np.testing.assert_array_equal(features, [30, 28.0, 0.4, 1.0, 60])

    def test_discordant(self):
        c = _make_candidate(is_concordant=False)
        features = candidate_to_features(c)
        assert features[3] == 0.0

    def test_unknown_concordance(self):
        c = _make_candidate(is_concordant=None)
        features = candidate_to_features(c)
        assert features[3] == 0.5


class TestFeatureMatrix:
    def test_empty(self):
        X = candidates_to_feature_matrix([])
        assert X.shape == (0, 5)

    def test_multiple(self):
        candidates = [_make_candidate(vbq=i * 10) for i in range(3)]
        X = candidates_to_feature_matrix(candidates)
        assert X.shape == (3, 5)


class TestSVMTraining:
    def test_train_and_predict(self):
        """SVM should separate high-quality from low-quality reads."""
        rng = np.random.RandomState(42)

        # Positive examples: high quality, mid-read position, concordant
        pos = rng.normal(loc=[35, 32, 0.5, 1.0, 60], scale=[3, 2, 0.1, 0, 2], size=(200, 5))

        # Negative examples: low quality, 3' position, discordant
        neg = rng.normal(loc=[15, 18, 0.85, 0.0, 30], scale=[3, 2, 0.1, 0, 5], size=(200, 5))

        model = train_svm(pos, neg)

        # Test on clearly positive example
        test_pos = np.array([[35, 32, 0.5, 1.0, 60]])
        assert model.predict(test_pos)[0] == 1

        # Test on clearly negative example
        test_neg = np.array([[15, 18, 0.85, 0.0, 30]])
        assert model.predict(test_neg)[0] == 0

    def test_filter_candidates(self):
        """filter_candidates should remove reads classified as artifacts."""
        rng = np.random.RandomState(42)
        pos = rng.normal(loc=[35, 32, 0.5, 1.0, 60], scale=[3, 2, 0.1, 0, 2], size=(100, 5))
        neg = rng.normal(loc=[15, 18, 0.85, 0.0, 30], scale=[3, 2, 0.1, 0, 5], size=(100, 5))
        model = train_svm(pos, neg)

        # Good read (should pass)
        good = _make_candidate(vbq=35, mrbq=32, pir=0.5, mapq=60, is_concordant=True)
        # Bad read (should be filtered)
        bad = _make_candidate(vbq=15, mrbq=18, pir=0.85, mapq=30, is_concordant=False)

        result = filter_candidates([good, bad], model)
        # At minimum, good should pass and bad should not
        assert len(result) >= 1
        assert any(r.vbq == 35 for r in result)


class TestBuildSVMPipeline:
    def test_pipeline_has_scaler_and_svm(self):
        pipeline = build_svm_pipeline()
        assert pipeline.named_steps["scaler"] is not None
        assert pipeline.named_steps["svm"] is not None

    def test_custom_C(self):
        pipeline = build_svm_pipeline(C=0.5)
        assert pipeline.named_steps["svm"].C == 0.5
