"""Probability calibration for per-fragment classifier scores.

The PPM numerator is a *sum of per-fragment probabilities* (score.sum_supporting_reads),
so those probabilities must be honest: if the model assigns 0.3 to a set of fragments,
about 30% of them should be true ctDNA. Raw gradient-boosting scores are not calibrated
in that sense, so we fit an isotonic regression on a held-out validation fold mapping raw
score -> empirical P(tumor).

Isotonic (monotone, non-parametric) is preferred over Platt scaling here because the
score->probability relationship need not be sigmoidal and we usually have enough
validation fragments to support it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from sklearn.isotonic import IsotonicRegression

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray


class IsotonicCalibrator:
    """Monotone mapping from raw classifier score to calibrated probability."""

    def __init__(self, isotonic: IsotonicRegression) -> None:
        self._isotonic = isotonic

    @classmethod
    def fit(cls, scores: ArrayLike, labels: ArrayLike) -> IsotonicCalibrator:
        """Fit the calibrator on held-out (score, label) pairs.

        Args:
            scores: Raw model scores on the validation fold.
            labels: Binary ground-truth labels (1 = ctDNA, 0 = artifact).

        Returns:
            A fitted IsotonicCalibrator.
        """
        isotonic = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        isotonic.fit(np.asarray(scores, dtype=np.float64), np.asarray(labels, dtype=np.float64))
        return cls(isotonic)

    def apply(self, scores: ArrayLike) -> NDArray[np.float64]:
        """Map raw scores to calibrated probabilities in [0, 1]."""
        return self._isotonic.predict(np.asarray(scores, dtype=np.float64))
