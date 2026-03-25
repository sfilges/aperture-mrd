"""Noise model construction from control plasma samples.

The noise model captures the expected background detection rate when a
patient-specific compendium is applied to plasma samples with no tumor
(TF=0). This is done via cross-patient analysis: applying one patient's
compendium to other patients' control plasma, or to a dedicated panel of
healthy donor plasma samples.

The resulting noise distribution (mean μ, std σ) defines the z-score
denominator for detection calls.

References:
    Zviran et al. Nature Medicine (2020), Methods:
        "Sequencing error suppression" and "Plasma SNV-based ctDNA detection"
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np

from aperture_snv.extract import (
    CandidateRead,
    CompendiumSite,
    extract_all_candidates,
    load_compendium,
    resolve_paired_end_concordance,
)
from aperture_snv.svm import filter_candidates

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

    from sklearn.pipeline import Pipeline


@dataclasses.dataclass(frozen=True, slots=True)
class NoiseProfile:
    """Noise statistics from control plasma samples for a patient compendium."""

    compendium_size: int  # N: number of sites in the compendium
    n_controls: int  # number of control samples evaluated
    detection_rates: tuple[float, ...]  # per-control detection rates
    noise_mean: float  # μ: mean detection rate across controls
    noise_std: float  # σ: std of detection rate across controls
    noise_rate_per_read: float  # μ_read: mean errors per read evaluated


def compute_noise_for_control(
    candidates: dict[CompendiumSite, list[CandidateRead]],
    svm_model: Pipeline,
) -> tuple[int, int, int]:
    """Compute noise statistics for a single control sample.

    Args:
        candidates: Extracted candidates at compendium sites (already concordance-resolved).
        svm_model: Trained SVM pipeline for filtering.

    Returns:
        Tuple of (detected_variants, total_reads_evaluated, total_sites_checked).
    """
    total_detected = 0
    total_reads = 0
    total_sites = 0

    for site, reads in candidates.items():
        total_sites += 1

        # Count reads that support the alt allele
        alt_reads = [r for r in reads if r.read_base == site.alt]

        if alt_reads:
            passed = filter_candidates(alt_reads, svm_model)
            total_detected += len(passed)

        total_reads += len(reads)

    return total_detected, total_reads, total_sites


def build_noise_profile(
    compendium_path: str | Path,
    control_cram_paths: Sequence[str | Path],
    reference_path: str | Path,
    svm_model: Pipeline,
) -> NoiseProfile:
    """Build a noise profile by applying a compendium to control plasma samples.

    For each control sample (healthy donor or cross-patient plasma, TF=0):
    1. Extract candidates at all compendium sites
    2. Resolve paired-end concordance
    3. Apply SVM filter
    4. Compute detection rate

    The resulting distribution of detection rates defines the noise model.

    Args:
        compendium_path: Path to the patient SNV compendium VCF.
        control_cram_paths: Paths to control plasma CRAM/BAM files (n≥8 recommended).
        reference_path: Path to reference FASTA.
        svm_model: Trained SVM pipeline.

    Returns:
        NoiseProfile with noise statistics.

    Raises:
        ValueError: If fewer than 2 control samples are provided.
    """
    if len(control_cram_paths) < 2:
        raise ValueError(
            f"At least 2 control samples required, got {len(control_cram_paths)}. "
            "Recommend n≥8 for stable noise estimates (n≥30 for publication)."
        )

    compendium = load_compendium(compendium_path)
    compendium_size = len(compendium)

    detection_rates: list[float] = []
    noise_rates_per_read: list[float] = []

    for cram_path in control_cram_paths:
        candidates = extract_all_candidates(
            cram_path=cram_path,
            compendium_path=compendium_path,
            reference_path=reference_path,
        )
        candidates = resolve_paired_end_concordance(candidates)

        detected, total_reads, _total_sites = compute_noise_for_control(candidates, svm_model)

        det_rate = detected / total_reads if total_reads > 0 else 0.0
        noise_per_read = detected / total_reads if total_reads > 0 else 0.0

        detection_rates.append(det_rate)
        noise_rates_per_read.append(noise_per_read)

    rates_array = np.array(detection_rates)

    return NoiseProfile(
        compendium_size=compendium_size,
        n_controls=len(control_cram_paths),
        detection_rates=tuple(detection_rates),
        noise_mean=float(np.mean(rates_array)),
        noise_std=float(np.std(rates_array, ddof=1)),
        noise_rate_per_read=float(np.mean(noise_rates_per_read)),
    )
