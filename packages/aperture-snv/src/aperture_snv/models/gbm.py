"""Gradient-boosted per-fragment classifier: P(tumor-derived | fragment).

A LightGBM model over the generic fragment features (extract.FEATURE_COLS). This is
the v1 model — sample-efficient, natively handles the mixed tabular features, calibrates
cleanly, and (unlike the fragment-tensor CNN, deferred to v2) needs no large labeled set.

It is trained on fragment-intrinsic features only and never on tumor-locus membership,
so the same model serves the agnostic channel later without retraining (see DESIGN.md).
Its calibrated output is the per-fragment weight for the PPM numerator
(score.sum_supporting_reads).

Labels (compendium sites only): positives are alt-supporting fragments in high-TF /
tumor-confirmed plasma; negatives are alt-supporting fragments at the same sites in TF=0
control plasma (artifacts by construction).
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import joblib
import numpy as np
from lightgbm import LGBMClassifier

from aperture_snv.extract import FEATURE_COLS

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from numpy.typing import ArrayLike, NDArray

    from aperture_snv.calibrate import IsotonicCalibrator
    from aperture_snv.extract import FragmentFeatures


# concordance is categorical; encode as small integer codes. With only three levels a
# tree can partition them with numeric splits, so we do not mark it categorical (which
# avoids a pandas dependency and the categorical_feature deprecation). Its monotone
# constraint is 0 — the code order is not a trust order.
_CONCORDANCE_CODES = {"concordant": 0, "discordant": 1, "single": 2}

# Monotone constraints where the direction of evidence is unambiguous, to suppress
# spurious (batch-driven) fits. Features not listed are unconstrained (0).
_MONOTONE_CONSTRAINTS = {
    "vbq": 1,  # higher variant base quality -> more trustworthy
    "mrbq": 1,  # higher mean base quality -> more trustworthy
    "mapq": 1,  # higher mapping quality -> more trustworthy
    "dist_3prime": 1,  # further from the 3' end -> less error-prone
    "n_low_bq": -1,  # more low-quality bases -> less trustworthy
    "edit_distance": -1,  # more other mismatches -> likely misaligned
    "softclip_len": -1,  # more soft-clipping -> less trustworthy
    # pir, concordance, fragment_length, is_proper_pair: non-monotonic / categorical -> 0
}

MONOTONE_CONSTRAINTS: list[int] = [_MONOTONE_CONSTRAINTS.get(col, 0) for col in FEATURE_COLS]


def _encode(fragment: FragmentFeatures) -> list[float]:
    """Encode one fragment into the FEATURE_COLS-ordered numeric feature vector."""
    row: list[float] = []
    for col in FEATURE_COLS:
        value = getattr(fragment, col)
        if col == "concordance":
            row.append(float(_CONCORDANCE_CODES[value]))
        else:
            row.append(float(value))  # ints and bools coerce cleanly
    return row


def fragments_to_matrix(fragments: Sequence[FragmentFeatures]) -> NDArray[np.float64]:
    """Encode fragments into a (n_fragments, n_features) matrix in FEATURE_COLS order."""
    if not fragments:
        return np.empty((0, len(FEATURE_COLS)), dtype=np.float64)
    return np.asarray([_encode(f) for f in fragments], dtype=np.float64)


def default_params() -> dict:
    """Default LightGBM hyperparameters for the fragment classifier."""
    return {
        "objective": "binary",
        "n_estimators": 300,
        "learning_rate": 0.05,
        "num_leaves": 31,
        "min_child_samples": 50,
        "subsample": 0.8,
        "subsample_freq": 1,
        "colsample_bytree": 0.8,
        "reg_lambda": 1.0,
        "class_weight": "balanced",
        "monotone_constraints": MONOTONE_CONSTRAINTS,
        "random_state": 42,
        "n_jobs": -1,
        "verbose": -1,
    }


class GBMModel:
    """Trained per-fragment classifier with optional probability calibration."""

    def __init__(
        self,
        classifier: LGBMClassifier,
        calibrator: IsotonicCalibrator | None = None,
    ) -> None:
        self.classifier = classifier
        self.calibrator = calibrator

    @classmethod
    def train(
        cls,
        positive_fragments: Sequence[FragmentFeatures],
        negative_fragments: Sequence[FragmentFeatures],
        **param_overrides,
    ) -> GBMModel:
        """Train on labeled fragments.

        Args:
            positive_fragments: True ctDNA fragments (label 1).
            negative_fragments: Artifact fragments (label 0).
            **param_overrides: LightGBM hyperparameters overriding `default_params`.

        Returns:
            An (uncalibrated) trained GBMModel. Call `calibrate` on a held-out fold
            before using the probabilities in the PPM numerator.
        """
        if not positive_fragments or not negative_fragments:
            raise ValueError("both positive and negative fragments are required to train")

        x_pos = fragments_to_matrix(positive_fragments)
        x_neg = fragments_to_matrix(negative_fragments)
        features = np.vstack([x_pos, x_neg])
        labels = np.concatenate([np.ones(len(x_pos)), np.zeros(len(x_neg))])

        classifier = LGBMClassifier(**(default_params() | param_overrides))
        classifier.fit(features, labels)
        return cls(classifier)

    def raw_scores(self, fragments: Sequence[FragmentFeatures]) -> NDArray[np.float64]:
        """Uncalibrated P(tumor) from the classifier."""
        if not fragments:
            return np.empty(0, dtype=np.float64)
        with warnings.catch_warnings():
            # We predict on unnamed numpy arrays; the trained booster carries LightGBM's
            # auto-generated feature names, which triggers a benign sklearn warning.
            warnings.filterwarnings("ignore", message="X does not have valid feature names")
            return self.classifier.predict_proba(fragments_to_matrix(fragments))[:, 1]

    def predict_proba(self, fragments: Sequence[FragmentFeatures]) -> NDArray[np.float64]:
        """Calibrated P(tumor) per fragment (raw scores if not yet calibrated).

        This is the per-fragment weight for the PPM numerator: pass fragments in a batch
        and multiply/sum against `supports_alt` rather than calling per fragment.
        """
        scores = self.raw_scores(fragments)
        if self.calibrator is not None and len(scores):
            return self.calibrator.apply(scores)
        return scores

    def calibrate(
        self,
        fragments: Sequence[FragmentFeatures],
        labels: ArrayLike,
    ) -> None:
        """Fit the probability calibrator on a held-out validation fold."""
        from aperture_snv.calibrate import IsotonicCalibrator

        self.calibrator = IsotonicCalibrator.fit(self.raw_scores(fragments), labels)

    def save(self, path: str | Path) -> None:
        """Serialize the model (classifier + calibrator) to disk."""
        joblib.dump(self, path)

    @staticmethod
    def load(path: str | Path) -> GBMModel:
        """Load a serialized GBMModel from disk."""
        return joblib.load(path)
