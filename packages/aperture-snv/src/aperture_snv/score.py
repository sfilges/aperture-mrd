"""PPM-based mutational-load estimation and detection z-scoring.

Instead of inverting a coverage/binomial model to a "true" tumor fraction (fragile at
low WGS depth and undefined in the tumor-agnostic case), we report an interpretable
mutational-load proxy modeled on Caris Assure's "tumor parts per million" (PPM):

    PPM = tumor-supporting reads / total aligned reads * 1e6

The numerator is the model signal — for a calibrated per-fragment classifier, the
expected count ``sum(p_i)`` over alt-supporting fragments; with no model yet, a plain
count of alt-supporting fragments. The denominator is a pure sequencing-depth
normalizer, so the same formula works in both tumor-informed (compendium sites) and
agnostic (called candidates) settings — only the site set feeding the numerator changes.

Detection is a separate step: z-score the sample PPM against a control-panel PPM
distribution (see noise.py). PPM itself needs no controls.

References:
    Abraham et al. Sci Rep (2025), Caris Assure: "Mutationome for ABCDai" (tumor PPM).
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from aperture_snv.extract import FragmentFeatures

PPM_SCALE = 1_000_000


@dataclasses.dataclass(frozen=True, slots=True)
class SNVScore:
    """Result of PPM-based ctDNA scoring for a single plasma sample."""

    supporting_reads: float  # expected tumor-supporting reads (sum of per-fragment weights)
    total_aligned_reads: int  # depth normalizer
    ppm: float  # supporting_reads / total_aligned_reads * 1e6
    z_score: float  # PPM vs control-panel noise distribution
    is_detected: bool  # whether z_score exceeds the threshold


def sum_supporting_reads(
    fragments: Iterable[FragmentFeatures],
    weight_fn: Callable[[FragmentFeatures], float] | None = None,
) -> float:
    """Expected tumor-supporting reads across alt-supporting fragments.

    Args:
        fragments: Candidate fragments at the sites of interest.
        weight_fn: Maps a fragment to its P(tumor) in [0, 1]. Defaults to a weight of
            1.0 per alt-supporting fragment (a plain count) until a calibrated model
            exists; pass the model's per-fragment probability once available.

    Returns:
        Sum of weights over alt-supporting fragments.
    """
    total = 0.0
    for fragment in fragments:
        if fragment.supports_alt:
            total += 1.0 if weight_fn is None else weight_fn(fragment)
    return total


def compute_ppm(supporting_reads: float, total_aligned_reads: int) -> float:
    """Tumor-supporting reads per million aligned reads.

    Args:
        supporting_reads: Expected tumor-supporting reads (numerator).
        total_aligned_reads: Total mapped reads in the sample (denominator).

    Returns:
        PPM, or 0.0 if the sample has no aligned reads.
    """
    if total_aligned_reads <= 0:
        return 0.0
    return supporting_reads / total_aligned_reads * PPM_SCALE


def compute_z_score(ppm: float, noise_mean: float, noise_std: float) -> float:
    """Z-score of a sample PPM against the control-panel noise distribution.

    Args:
        ppm: Observed PPM in the patient plasma.
        noise_mean: Mean PPM across control (TF=0) samples.
        noise_std: Standard deviation of PPM across controls.

    Returns:
        Z-score. If the control std is zero, returns +inf above the mean else 0.0.
    """
    if noise_std == 0:
        return 0.0 if ppm <= noise_mean else float("inf")
    return (ppm - noise_mean) / noise_std


def score_sample(
    supporting_reads: float,
    total_aligned_reads: int,
    noise_mean: float,
    noise_std: float,
    z_threshold: float = 1.2,
) -> SNVScore:
    """Compute the PPM mutational-load score and detection call for a plasma sample.

    Args:
        supporting_reads: Expected tumor-supporting reads (from sum_supporting_reads).
        total_aligned_reads: Total mapped reads in the plasma sample.
        noise_mean: Mean PPM across control samples.
        noise_std: Std of PPM across controls.
        z_threshold: Z-score threshold for a positive detection call.

    Returns:
        SNVScore with PPM, z-score, and detection call.
    """
    ppm = compute_ppm(supporting_reads, total_aligned_reads)
    z = compute_z_score(ppm, noise_mean, noise_std)
    return SNVScore(
        supporting_reads=supporting_reads,
        total_aligned_reads=total_aligned_reads,
        ppm=ppm,
        z_score=z,
        is_detected=z > z_threshold,
    )
