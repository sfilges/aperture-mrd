"""Noise model: control-panel PPM distribution for detection z-scoring.

Applying the patient compendium (tumor-informed) — or the same candidate generation
(agnostic) — to control plasma with no tumor (TF=0) yields a background PPM per control.
The mean and standard deviation of that distribution are the denominator of the
detection z-score (see score.compute_z_score).

The fragments are filtered and weighted with the *same* config and model used to score
the patient sample, so the PPMs are comparable.

References:
    Zviran et al. Nature Medicine (2020), Methods: cross-patient control analysis.
    Abraham et al. Sci Rep (2025), Caris Assure: tumor PPM.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np

from aperture_snv.extract import count_mapped_reads, iter_fragments
from aperture_snv.filters import DEFAULT_FILTERS
from aperture_snv.score import compute_ppm, sum_supporting_reads

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from pathlib import Path

    from aperture_snv.extract import FragmentFeatures
    from aperture_snv.filters import FilterConfig


@dataclasses.dataclass(frozen=True, slots=True)
class NoiseProfile:
    """Background PPM statistics from control plasma samples for a compendium."""

    n_controls: int  # number of control samples evaluated
    control_ppms: tuple[float, ...]  # per-control background PPM
    noise_mean: float  # μ: mean control PPM
    noise_std: float  # σ: std of control PPM


def _control_ppm(
    control_path: str | Path,
    compendium_path: str | Path,
    reference_path: str | Path | None,
    weight_fn: Callable[[FragmentFeatures], float] | None,
    config: FilterConfig,
) -> float:
    """Background PPM for one control sample (TF=0)."""
    supporting = 0.0
    for _site, fragments in iter_fragments(
        alignment_file_path=control_path,
        compendium_vcf_path=compendium_path,
        reference_fasta_path=reference_path,
        config=config,
    ):
        supporting += sum_supporting_reads(fragments, weight_fn)
    total_reads = count_mapped_reads(control_path, reference_path)
    return compute_ppm(supporting, total_reads)


def build_noise_profile(
    compendium_path: str | Path,
    control_paths: Sequence[str | Path],
    reference_path: str | Path | None = None,
    weight_fn: Callable[[FragmentFeatures], float] | None = None,
    config: FilterConfig = DEFAULT_FILTERS,
) -> NoiseProfile:
    """Build the background PPM distribution from control (TF=0) plasma samples.

    Args:
        compendium_path: Path to the patient SNV compendium VCF.
        control_paths: Control plasma CRAM/BAM files (healthy donor or cross-patient).
        reference_path: Reference FASTA (required for CRAM).
        weight_fn: Per-fragment P(tumor) weighter; must match what is used to score the
            patient sample. Defaults to a plain count of alt-supporting fragments.
        config: Fragment-stage filter thresholds; must match the sample's.

    Returns:
        NoiseProfile with the per-control PPMs and their mean/std.

    Raises:
        ValueError: If fewer than 2 control samples are provided.
    """
    if len(control_paths) < 2:
        raise ValueError(
            f"At least 2 control samples required, got {len(control_paths)}. "
            "Recommend n≥8 for stable noise estimates (n≥30 for publication)."
        )

    ppms = [
        _control_ppm(path, compendium_path, reference_path, weight_fn, config)
        for path in control_paths
    ]
    arr = np.asarray(ppms)
    return NoiseProfile(
        n_controls=len(control_paths),
        control_ppms=tuple(ppms),
        noise_mean=float(arr.mean()),
        noise_std=float(arr.std(ddof=1)),
    )
