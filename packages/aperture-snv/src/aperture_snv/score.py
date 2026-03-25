"""Binomial model for tumor fraction estimation and detection z-scoring.

After SVM filtering, the number of detected variants at compendium sites
in plasma follows a binomial model. This module computes the tumor fraction
(TF) and a detection z-score against a control noise distribution.

The key insight: at low TF, each compendium site has independent probability
of being sampled, governed by coverage depth. Integrating signal across
thousands of sites overcomes the per-site sampling limitation.

References:
    Zviran et al. Nature Medicine (2020), Methods: Equations 4 and 5.
"""

from __future__ import annotations

import dataclasses


@dataclasses.dataclass(frozen=True, slots=True)
class SNVScore:
    """Result of SNV-based ctDNA scoring for a single plasma sample."""

    detected_variants: int  # M: number of compendium sites with supporting reads
    total_sites: int  # N: total sites in patient compendium
    mean_coverage: float  # cov: mean coverage at compendium sites in plasma
    noise_rate: float  # μ: mean artifact noise rate (errors per read evaluated)
    total_reads_evaluated: int  # R: total reads covering compendium sites
    tumor_fraction: float  # estimated TF
    detection_rate: float  # M / R
    z_score: float  # z-score against control noise distribution
    is_detected: bool  # whether z_score exceeds threshold


def estimate_tumor_fraction(
    M: int,
    N: int,
    cov: float,
    mu: float,
    R: int,
) -> float:
    """Estimate tumor fraction from detected variant count.

    Implements Eq. 5 from Zviran et al.:
        TF = 1 - (1 - [M - μ*R] / N) ^ (1/cov)

    At very low TF, the relationship between detected mutations and TF
    is approximately linear: M ≈ N * TF * cov + μ*R

    Args:
        M: Number of compendium sites with detected variants in plasma.
        N: Total number of SNVs in the patient compendium.
        cov: Mean coverage depth at compendium sites in the plasma sample.
        mu: Noise rate (errors per read evaluated) from control samples.
        R: Total reads evaluated at compendium sites.

    Returns:
        Estimated tumor fraction (0–1). Returns 0.0 if signal is below noise.
    """
    if N == 0 or cov == 0:
        return 0.0

    signal = M - mu * R
    if signal <= 0:
        return 0.0

    ratio = signal / N
    if ratio >= 1.0:
        return 1.0

    tf = 1.0 - (1.0 - ratio) ** (1.0 / cov)
    return max(0.0, min(1.0, tf))


def compute_detection_rate(M: int, R: int) -> float:
    """Compute the detection rate: fraction of evaluated reads supporting a variant.

    Args:
        M: Number of detected variants.
        R: Total reads evaluated.

    Returns:
        Detection rate M/R.
    """
    if R == 0:
        return 0.0
    return M / R


def compute_z_score(
    detection_rate: float,
    noise_mean: float,
    noise_std: float,
) -> float:
    """Compute z-score of detection rate against control noise distribution.

    z = (det_rate - μ) / σ

    The noise distribution is estimated by applying the patient-specific
    compendium to control plasma samples (cross-patient analysis) or to
    other patients' plasma samples.

    Args:
        detection_rate: Observed detection rate in patient plasma.
        noise_mean: Mean detection rate in control samples (μ).
        noise_std: Standard deviation of detection rate in controls (σ).

    Returns:
        z-score. Values > 1.2 typically indicate detection at >95% specificity.
    """
    if noise_std == 0:
        return 0.0 if detection_rate <= noise_mean else float("inf")
    return (detection_rate - noise_mean) / noise_std


def score_sample(
    detected_variants: int,
    total_sites: int,
    mean_coverage: float,
    noise_rate: float,
    total_reads_evaluated: int,
    noise_mean: float,
    noise_std: float,
    z_threshold: float = 1.2,
) -> SNVScore:
    """Compute full SNV-based ctDNA score for a plasma sample.

    Combines tumor fraction estimation and z-score detection into a
    single result object.

    Args:
        detected_variants: Number of compendium sites with supporting reads post-SVM.
        total_sites: Total SNVs in patient compendium.
        mean_coverage: Mean coverage at compendium sites in plasma.
        noise_rate: Per-read noise rate (μ) from control samples.
        total_reads_evaluated: Total reads covering compendium sites.
        noise_mean: Mean detection rate across control samples.
        noise_std: Std deviation of detection rate across controls.
        z_threshold: z-score threshold for positive detection (default 1.2).

    Returns:
        SNVScore with TF estimate, z-score, and detection call.
    """
    tf = estimate_tumor_fraction(
        M=detected_variants,
        N=total_sites,
        cov=mean_coverage,
        mu=noise_rate,
        R=total_reads_evaluated,
    )
    det_rate = compute_detection_rate(detected_variants, total_reads_evaluated)
    z = compute_z_score(det_rate, noise_mean, noise_std)

    return SNVScore(
        detected_variants=detected_variants,
        total_sites=total_sites,
        mean_coverage=mean_coverage,
        noise_rate=noise_rate,
        total_reads_evaluated=total_reads_evaluated,
        tumor_fraction=tf,
        detection_rate=det_rate,
        z_score=z,
        is_detected=z > z_threshold,
    )
