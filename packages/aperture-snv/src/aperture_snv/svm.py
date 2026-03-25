"""Read-centric SVM error suppression for cfDNA SNV candidates.

Classifies individual reads as likely true somatic ctDNA vs sequencing artifact
using a 5-feature linear SVM. This is fundamentally different from standard
locus-centric variant calling: at TFs below 10⁻⁴, at most one supporting read
per site is expected, so read-level classification is the appropriate paradigm.

Features:
    1. VBQ  — Variant Base Quality at the SNV position
    2. MRBQ — Mean Read Base Quality (correlates with PIR; 3' quality drop)
    3. PIR  — Position in Read (0–1 fraction; artifacts cluster at 3' end)
    4. CONC — R1/R2 paired-end concordance (1.0 if concordant, 0.5 if unknown, 0.0 if discordant)
    5. MAPQ — Mapping Quality

References:
    Zviran et al. Nature Medicine (2020), Methods: "Sequencing error suppression"
    Extended Data Fig. 2d: SVM outperforms random forest on this feature set.
"""

from __future__ import annotations

import pickle
from typing import TYPE_CHECKING

import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from numpy.typing import NDArray

    from aperture_snv.extract import CandidateRead


FEATURE_NAMES = ("vbq", "mrbq", "pir", "concordance", "mapq")


def candidate_to_features(candidate: CandidateRead) -> NDArray[np.float64]:
    """Extract the 5-element feature vector from a CandidateRead.

    Args:
        candidate: A CandidateRead with alignment features.

    Returns:
        Array of shape (5,) with [VBQ, MRBQ, PIR, concordance, MAPQ].
    """
    if candidate.is_concordant is True:
        conc = 1.0
    elif candidate.is_concordant is None:
        conc = 0.5
    else:
        conc = 0.0

    return np.array(
        [candidate.vbq, candidate.mrbq, candidate.pir, conc, candidate.mapq],
        dtype=np.float64,
    )


def candidates_to_feature_matrix(
    candidates: Sequence[CandidateRead],
) -> NDArray[np.float64]:
    """Convert a sequence of CandidateReads to a feature matrix.

    Args:
        candidates: Sequence of CandidateRead objects.

    Returns:
        Array of shape (n_reads, 5).
    """
    if not candidates:
        return np.empty((0, 5), dtype=np.float64)
    return np.stack([candidate_to_features(c) for c in candidates])


def build_svm_pipeline(C: float = 1.0) -> Pipeline:
    """Create an untrained SVM pipeline with standardization.

    Architecture matches Zviran et al.: linear SVM with hinge loss,
    L2 regularization, C=1.0.

    Args:
        C: Regularization parameter (default 1.0 per paper).

    Returns:
        Untrained sklearn Pipeline (StandardScaler → LinearSVC).
    """
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "svm",
                LinearSVC(
                    C=C,
                    loss="hinge",
                    penalty="l2",
                    max_iter=10000,
                    class_weight="balanced",
                    random_state=42,
                ),
            ),
        ]
    )


def train_svm(
    positive_features: NDArray[np.float64],
    negative_features: NDArray[np.float64],
    C: float = 1.0,
) -> Pipeline:
    """Train the SVM classifier on labeled read features.

    Positive examples: reads at known true somatic SNV positions (from high-TF
    plasma or tumor WGS at compendium sites).
    Negative examples: reads at compendium positions in control plasma (TF=0),
    representing sequencing artifacts.

    Training set per Zviran et al.: 8 control samples, 10,000 variants per class.

    Args:
        positive_features: Feature matrix of true ctDNA reads, shape (n, 5).
        negative_features: Feature matrix of artifact reads, shape (m, 5).
        C: Regularization parameter.

    Returns:
        Trained Pipeline ready for prediction.
    """
    X = np.vstack([positive_features, negative_features])
    y = np.concatenate([np.ones(len(positive_features)), np.zeros(len(negative_features))])
    pipeline = build_svm_pipeline(C=C)
    pipeline.fit(X, y)
    return pipeline


def filter_candidates(
    candidates: Sequence[CandidateRead],
    model: Pipeline,
) -> list[CandidateRead]:
    """Apply the trained SVM to filter candidate reads.

    Reads classified as artifacts (class 0) are removed.
    Reads classified as likely true ctDNA (class 1) are retained.

    Args:
        candidates: Sequence of CandidateRead objects.
        model: Trained SVM pipeline.

    Returns:
        List of CandidateReads that pass the SVM filter.
    """
    if not candidates:
        return []

    X = candidates_to_feature_matrix(candidates)
    predictions = model.predict(X)
    return [c for c, pred in zip(candidates, predictions, strict=True) if pred == 1]


def save_model(model: Pipeline, path: str | Path) -> None:
    """Serialize a trained SVM model to disk."""
    with open(path, "wb") as f:
        pickle.dump(model, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_model(path: str | Path) -> Pipeline:
    """Load a serialized SVM model from disk."""
    with open(path, "rb") as f:
        return pickle.load(f)  # noqa: S301 — trusted internal model files only
