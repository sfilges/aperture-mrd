"""Fragment-stage hard filters, applied inline during extraction.

These are the *permissive* floors that remove unusable reads/fragments before the
model sees them — structurally broken alignments and clearly untrustworthy bases —
not aggressive cuts on the quantities we also feed the model as features. The model
does the fine discrimination in the gray zone; these filters just remove garbage that
would otherwise be training noise.

The single FilterConfig must be applied identically across the scored plasma, the
control plasma (noise model), and the training extraction, or the supporting-read
counts (and the PPM z-score built on them) are not comparable across samples.

Filters live here as pure predicates and are called inline by `extract` so that
failing reads never materialize into fragments — important for genome-wide agnostic
scanning where the candidate set is large.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from aperture_snv.extract import FragmentFeatures


@dataclasses.dataclass(frozen=True, slots=True)
class FilterConfig:
    """Permissive fragment-stage thresholds.

    Attributes:
        min_mapq: Minimum mapping quality (enforced by the pileup for efficiency).
        min_bq: Minimum base quality at the variant position (enforced by the pileup).
        min_insert: Minimum fragment insert size (bp).
        max_insert: Maximum fragment insert size (bp); cfDNA is nucleosome-sized.
    """

    min_mapq: int = 20
    min_bq: int = 20
    min_insert: int = 40
    max_insert: int = 240


DEFAULT_FILTERS = FilterConfig()


class _Alignment(Protocol):
    """Structural flags read by `read_flags_ok` (a subset of pysam.AlignedSegment)."""

    is_unmapped: bool
    is_secondary: bool
    is_supplementary: bool
    is_duplicate: bool
    is_qcfail: bool


def read_flags_ok(alignment: _Alignment) -> bool:
    """True if the alignment is structurally usable.

    Drops unmapped, secondary, supplementary, duplicate, and QC-fail reads. MAPQ and
    base-quality floors are enforced at the pileup (see `extract_fragments_at_site`),
    so they are not re-checked here.
    """
    return not (
        alignment.is_unmapped
        or alignment.is_secondary
        or alignment.is_supplementary
        or alignment.is_duplicate
        or alignment.is_qcfail
    )


def fragment_passes(fragment: FragmentFeatures, config: FilterConfig) -> bool:
    """True if an assembled fragment passes the fragment-level floors (insert size)."""
    return config.min_insert <= fragment.fragment_length <= config.max_insert
