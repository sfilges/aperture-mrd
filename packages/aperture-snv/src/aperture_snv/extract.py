"""Extract candidate SNV-supporting fragments from cfDNA at compendium loci.

Extracts reads from WGS of a liquid biopsy sample at mutations identified in
tumor-normal somatic variant calling from the same patient (SNV compendium) for
ultrasensitive detection of ctDNA.

The unit is the *fragment* (both mates of a read pair collapsed), not the single
read: for ~165 bp cfDNA inserts with 150 bp paired-end reads R1 and R2 overlap
heavily, and whether both overlapping mates agree at the site (concordance) is one
of the strongest generic signals separating a true mutation from a sequencing
error. Concordance is therefore computed here as a per-fragment *feature* rather
than a separate post-hoc pass.

At the tumor fractions relevant to MRD, at most one supporting fragment per locus
is expected, so the downstream model is read/fragment-centric ("is this fragment
tumor-derived?") rather than locus-centric ("do enough reads support a variant?").

Requires:
    - BAM/CRAM file (liquid biopsy)
    - Reference FASTA (required for CRAM)
    - Compendium VCF

Inspired by MRDetect (Zviran et al.) and MRD-EDGE (Widman et al.); see DESIGN.md.

References:
    Widman et al. Nature Medicine (2024), machine-learning-guided signal enrichment
    Zviran et al. Nature Medicine (2020), Methods: "Plasma cfDNA SNV identification"
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np
import pysam

from aperture_snv.filters import DEFAULT_FILTERS, fragment_passes, read_flags_ok

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from pathlib import Path

    from aperture_snv.filters import FilterConfig


# Bases at/below this Phred score count toward `n_low_bq`.
LOW_BQ_THRESHOLD = 20
_BAM_CSOFT_CLIP = 4

# Per-fragment feature columns consumed by the classifier. Allele identity
# (`supports_alt`, `read_base`) selects which fragments to score and is not itself
# a discriminative feature, so it is excluded here. Sequence-context and
# site-noise features are added in a later step (see DESIGN.md).
FEATURE_COLS: tuple[str, ...] = (
    "vbq",
    "mrbq",
    "n_low_bq",
    "mapq",
    "edit_distance",
    "pir",
    "dist_3prime",
    "concordance",
    "fragment_length",
    "softclip_len",
    "is_proper_pair",
)


@dataclasses.dataclass(frozen=True, slots=True)
class CompendiumSite:
    """A single SNV site from the patient compendium."""

    chrom: str
    pos: int  # 0-based
    ref: str
    alt: str


@dataclasses.dataclass(frozen=True, slots=True)
class MateObservation:
    """Values pulled from one alignment (mate) that covers a compendium site.

    This is the intermediate the pysam layer produces and the pure `build_fragment`
    logic consumes, so fragment assembly can be unit-tested without a BAM.
    """

    read_name: str
    is_read1: bool
    is_reverse: bool
    base_at_site: str
    bq_at_site: int
    query_position: int
    query_length: int
    mean_bq: float
    n_low_bq: int
    mapq: int
    edit_distance: int  # NM already adjusted to exclude the candidate's own mismatch
    softclip_len: int
    is_proper_pair: bool
    reference_start: int
    reference_end: int
    template_length: int


@dataclasses.dataclass(frozen=True, slots=True)
class FragmentFeatures:
    """A cfDNA fragment overlapping a compendium site, with per-fragment features."""

    # --- identity ---
    chrom: str
    pos: int  # 0-based
    ref: str
    alt: str
    frag_id: str

    # --- allele call (selects which fragments to score; not a model feature) ---
    supports_alt: bool
    read_base: str

    # --- generic per-fragment features (model inputs) ---
    vbq: int  # variant base quality at the site (max across mates covering it)
    mrbq: float  # mean base quality across the fragment
    n_low_bq: int  # bases <= LOW_BQ_THRESHOLD across the fragment
    mapq: int  # min mapping quality across mates
    edit_distance: int  # mismatches on the fragment EXCLUDING the candidate site
    pir: float  # fractional position in read (mean across mates, 0-1)
    dist_3prime: int  # bp from the 3' end, min across mates (error enriches at 3')
    concordance: int  # 0=concordant, 1=discordant, 2=single
    fragment_length: int  # insert size
    softclip_len: int  # total soft-clipped bases on the fragment
    is_proper_pair: bool
    n_mates_at_site: int  # 1 or 2


def load_compendium(vcf_path: str | Path) -> list[CompendiumSite]:
    """Load unique biallelic SNV sites from a compendium VCF.

    Indels and multi-nucleotide ALTs are skipped. Exact-duplicate sites are collapsed:
    variant annotation (one record per transcript/consequence) can emit the same
    variant on multiple rows, and loading it twice would double-count its supporting
    reads downstream.

    VCF (Variant Call Format) files use a 1-based, fully-closed coordinate system!
    record.pos: The record start position on chrom/contig (1-based inclusive).
    record.start: The record start position on chrom/contig (0-based inclusive).

    Args:
        vcf_path: Path to the compendium VCF (.vcf.gz with .tbi index).

    Returns:
        List of unique CompendiumSite objects, in first-seen order.
    """
    sites: list[CompendiumSite] = []
    seen: set[CompendiumSite] = set()
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            if len(record.ref) != 1:
                continue
            for alt in record.alts or []:
                if len(alt) != 1:
                    continue
                site = CompendiumSite(
                    chrom=record.chrom,
                    pos=record.start,  # 0-based; record.pos is 1-based (VCF POS)
                    ref=record.ref,
                    alt=alt,
                )
                if site not in seen:
                    seen.add(site)
                    sites.append(site)
    return sites


def _dist_from_3prime(mate: MateObservation) -> int:
    """Distance in bp of the variant base from the read's 3' end.

    Reverse-strand reads are stored reverse-complemented, so query position 0 is
    the 3' end; forward reads have their 3' end at the high query coordinate.
    """
    if mate.query_length <= 0:
        return 0
    if mate.is_reverse:
        return mate.query_position
    return mate.query_length - 1 - mate.query_position


def _fragment_length(mates: Sequence[MateObservation]) -> int:
    """Fragment insert size: prefer TLEN, fall back to the mates' reference span."""
    tlens = [abs(m.template_length) for m in mates if m.template_length]
    if tlens:
        return max(tlens)
    starts = [m.reference_start for m in mates]
    ends = [m.reference_end for m in mates]
    return max(ends) - min(starts)


def build_fragment(site: CompendiumSite, mates: Sequence[MateObservation]) -> FragmentFeatures:
    """Collapse the mate(s) covering a site into one fragment feature record.

    Args:
        site: The compendium SNV site.
        mates: One or two MateObservations sharing a read name at this site.

    Returns:
        A FragmentFeatures record.
    """
    if not mates:
        raise ValueError("build_fragment requires at least one mate")

    n_mates = len(mates)
    bases = {m.base_at_site for m in mates}
    if n_mates == 1:
        concordance = 2  # single
    elif len(bases) == 1:
        concordance = 0  # concordant
    else:
        concordance = 1  # discordant

    # Consensus base = the highest-quality observation (settles discordant pairs).
    best = max(mates, key=lambda m: m.bq_at_site)
    read_base = best.base_at_site
    supports_alt = read_base == site.alt

    covered = [m.query_position / m.query_length for m in mates if m.query_length > 0]
    pir = float(np.mean(covered)) if covered else 0.5

    return FragmentFeatures(
        chrom=site.chrom,
        pos=site.pos,
        ref=site.ref,
        alt=site.alt,
        frag_id=mates[0].read_name,
        supports_alt=supports_alt,
        read_base=read_base,
        vbq=max(m.bq_at_site for m in mates),
        mrbq=float(np.mean([m.mean_bq for m in mates])),
        n_low_bq=sum(m.n_low_bq for m in mates),
        mapq=min(m.mapq for m in mates),
        edit_distance=max(m.edit_distance for m in mates),
        pir=pir,
        dist_3prime=min(_dist_from_3prime(m) for m in mates),
        concordance=concordance,
        fragment_length=_fragment_length(mates),
        softclip_len=sum(m.softclip_len for m in mates),
        is_proper_pair=any(m.is_proper_pair for m in mates),
        n_mates_at_site=n_mates,
    )


def _raw_edit_distance(alignment: pysam.AlignedSegment) -> int:
    """Number of mismatches+indels vs reference (NM tag, or computed from alignment)."""
    if alignment.has_tag("NM"):
        return int(alignment.get_tag("NM"))
    # Fall back to counting substitutions from the aligned pairs (requires the
    # AlignmentFile to have been opened with a reference for MD/seq resolution).
    mismatches = 0
    for _qpos, _rpos, ref_base in alignment.get_aligned_pairs(with_seq=True):
        if _qpos is None or _rpos is None or ref_base is None:
            continue
        if ref_base.islower():  # lowercase ref base marks a mismatch in pysam
            mismatches += 1
    return mismatches


def _softclip_len(alignment: pysam.AlignedSegment) -> int:
    """Total soft-clipped bases across the alignment's CIGAR."""
    if alignment.cigartuples is None:
        return 0
    return sum(length for op, length in alignment.cigartuples if op == _BAM_CSOFT_CLIP)


def _observe_mate(
    site: CompendiumSite,
    pileup_read: pysam.PileupRead,
) -> MateObservation | None:
    """Build a MateObservation for one read covering the site, or None to skip.

    Structurally unusable reads (unmapped/secondary/supplementary/duplicate/QC-fail)
    are dropped here so they never become fragments; MAPQ and base-quality floors are
    enforced upstream by the pileup.
    """
    if pileup_read.is_del or pileup_read.is_refskip:
        return None
    alignment = pileup_read.alignment
    if not read_flags_ok(alignment):
        return None
    query_pos = pileup_read.query_position
    if query_pos is None or alignment.query_sequence is None:
        return None

    base_at_site = alignment.query_sequence[query_pos]
    qualities = alignment.query_qualities
    if qualities is not None:
        bq_at_site = int(qualities[query_pos])
        mean_bq = float(np.mean(qualities))
        n_low_bq = int(np.sum(np.asarray(qualities) <= LOW_BQ_THRESHOLD))
    else:
        bq_at_site, mean_bq, n_low_bq = 0, 0.0, 0

    # Exclude the candidate variant's own mismatch so a true variant read does not
    # look poorly aligned: NM counts the site mismatch whenever the base != ref.
    raw_nm = _raw_edit_distance(alignment)
    edit_distance = max(0, raw_nm - (1 if base_at_site != site.ref else 0))

    return MateObservation(
        read_name=alignment.query_name or "",
        is_read1=alignment.is_read1,
        is_reverse=alignment.is_reverse,
        base_at_site=base_at_site,
        bq_at_site=bq_at_site,
        query_position=query_pos,
        query_length=alignment.query_length or 0,
        mean_bq=mean_bq,
        n_low_bq=n_low_bq,
        mapq=alignment.mapping_quality,
        edit_distance=edit_distance,
        softclip_len=_softclip_len(alignment),
        is_proper_pair=alignment.is_proper_pair,
        reference_start=alignment.reference_start,
        reference_end=alignment.reference_end or alignment.reference_start,
        template_length=alignment.template_length,
    )


def extract_fragments_at_site(
    alignment_file: pysam.AlignmentFile,
    site: CompendiumSite,
    config: FilterConfig = DEFAULT_FILTERS,
) -> Iterator[FragmentFeatures]:
    """Yield one FragmentFeatures per read pair covering a compendium site.

    Reads covering the site are grouped by name; mates that both overlap the site
    (the R1/R2 overlap region) are collapsed into a single fragment. Filters are
    applied inline (MAPQ/BQ at the pileup, structural flags in `_observe_mate`,
    insert size after assembly) so failing reads never materialize.

    Args:
        alignment_file: Open pysam AlignmentFile (CRAM or BAM).
        site: The compendium SNV site to query.
        config: Fragment-stage filter thresholds.

    Yields:
        FragmentFeatures for each fragment passing the filters at the site.
    """
    by_name: dict[str, list[MateObservation]] = {}
    # pileup returns, for each base in the reference, the reads that map that position.
    # ignore_overlaps=False is required: with the default (True), pysam collapses
    # overlapping mate pairs by summing their base qualities and zeroing one mate,
    # which corrupts our per-mate base-quality features and pre-empts the R1/R2
    # concordance we compute ourselves in build_fragment.
    for pileup_column in alignment_file.pileup(
        contig=site.chrom,
        start=site.pos,
        stop=site.pos + 1,
        truncate=True,
        min_mapping_quality=config.min_mapq,
        min_base_quality=config.min_bq,
        ignore_overlaps=False,
        stepper="nofilter",
    ):
        for pileup_read in pileup_column.pileups:
            mate = _observe_mate(site, pileup_read)
            if mate is not None:
                by_name.setdefault(mate.read_name, []).append(mate)

    for mates in by_name.values():
        fragment = build_fragment(site, mates)
        if fragment_passes(fragment, config):
            yield fragment


def iter_fragments(
    alignment_file_path: str | Path,
    compendium_vcf_path: str | Path,
    reference_fasta_path: str | Path | None = None,
    config: FilterConfig = DEFAULT_FILTERS,
) -> Iterator[tuple[CompendiumSite, list[FragmentFeatures]]]:
    """Stream (site, fragments) one compendium site at a time.

    The alignment file is opened once and sites are processed lazily, so only one
    site's fragments are held at a time. This is the streaming core that keeps
    genome-wide (agnostic) extraction from materializing the whole candidate set;
    `extract_all_fragments` wraps it into a dict for the small tumor-informed case.

    Args:
        alignment_file_path: Path to plasma CRAM/BAM file.
        compendium_vcf_path: Path to patient SNV compendium VCF.
        reference_fasta_path: Path to reference FASTA. Required for CRAM (reads are
            reference-compressed); optional for BAM, which stores the read sequence.
        config: Fragment-stage filter thresholds.

    Yields:
        (site, fragments) tuples; `fragments` may be empty for a covered-but-filtered site.

    Raises:
        ValueError: If a CRAM file is given without a reference FASTA.
    """
    path_str = str(alignment_file_path)
    if path_str.endswith(".cram") and reference_fasta_path is None:
        raise ValueError("A reference FASTA is required to read CRAM input")

    reference_filename = str(reference_fasta_path) if reference_fasta_path is not None else None

    with pysam.AlignmentFile(path_str, reference_filename=reference_filename) as aln:
        for site in load_compendium(compendium_vcf_path):
            yield site, list(extract_fragments_at_site(aln, site, config=config))


def extract_all_fragments(
    alignment_file_path: str | Path,
    compendium_vcf_path: str | Path,
    reference_fasta_path: str | Path | None = None,
    config: FilterConfig = DEFAULT_FILTERS,
) -> dict[CompendiumSite, list[FragmentFeatures]]:
    """Extract candidate fragments at all compendium sites into a dict.

    Thin wrapper over `iter_fragments` that materializes every site; suitable for the
    tumor-informed case (hundreds of sites). For genome-wide agnostic scanning consume
    `iter_fragments` directly instead of building this dict.

    Args:
        alignment_file_path: Path to plasma CRAM/BAM file.
        compendium_vcf_path: Path to patient SNV compendium VCF.
        reference_fasta_path: Path to reference FASTA (required for CRAM).
        config: Fragment-stage filter thresholds.

    Returns:
        Dict mapping each CompendiumSite to its list of FragmentFeatures.
    """
    return dict(
        iter_fragments(
            alignment_file_path=alignment_file_path,
            compendium_vcf_path=compendium_vcf_path,
            reference_fasta_path=reference_fasta_path,
            config=config,
        )
    )


def count_mapped_reads(
    alignment_file_path: str | Path,
    reference_fasta_path: str | Path | None = None,
) -> int:
    """Total mapped reads in the file — the PPM depth normalizer.

    Read from the index statistics (fast; no full pass). Requires a .bai/.csi index.

    Args:
        alignment_file_path: Path to the CRAM/BAM file.
        reference_fasta_path: Reference FASTA (required for CRAM).

    Returns:
        Number of mapped reads.
    """
    reference_filename = str(reference_fasta_path) if reference_fasta_path is not None else None
    with pysam.AlignmentFile(str(alignment_file_path), reference_filename=reference_filename) as aln:
        return aln.mapped
