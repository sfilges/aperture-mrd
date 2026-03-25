"""Extract candidate SNV-supporting reads from plasma cfDNA at compendium loci.

Implements the MRDetect candidate extraction step: for each site in the
patient-specific SNV compendium, extract every read covering that position
from the plasma CRAM/BAM, along with per-read alignment features needed
for downstream SVM classification.

References:
    Zviran et al. Nature Medicine (2020), Methods: "Plasma cfDNA SNV identification"
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path
import pysam


@dataclasses.dataclass(frozen=True, slots=True)
class CandidateRead:
    """A single read overlapping a compendium SNV site with extracted features."""

    chrom: str
    pos: int  # 0-based
    ref: str
    alt: str
    read_name: str
    read_base: str
    vbq: int  # variant base quality
    mrbq: float  # mean read base quality
    pir: float  # position in read (fraction, 0–1)
    mapq: int  # mapping quality
    is_read1: bool
    mate_cigar: str | None  # for paired-end concordance
    is_concordant: bool | None  # True if mate also supports alt at this position
    is_proper_pair: bool
    insert_size: int


@dataclasses.dataclass(frozen=True, slots=True)
class CompendiumSite:
    """A single SNV site from the patient compendium."""

    chrom: str
    pos: int  # 0-based
    ref: str
    alt: str


def load_compendium(vcf_path: str | Path) -> list[CompendiumSite]:
    """Load SNV sites from a compendium VCF.

    Reads all PASS (or unfiltered) biallelic SNV records from the VCF.
    Indels and multiallelic sites are skipped.

    Args:
        vcf_path: Path to the compendium VCF (.vcf.gz with .tbi index).

    Returns:
        List of CompendiumSite objects.
    """
    sites: list[CompendiumSite] = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for rec in vcf:
            # Skip non-SNVs
            if len(rec.ref) != 1:
                continue
            for alt in rec.alts or []:
                if len(alt) != 1:
                    continue
                sites.append(
                    CompendiumSite(
                        chrom=rec.chrom,
                        pos=rec.pos,  # pysam VariantFile uses 0-based
                        ref=rec.ref,
                        alt=alt,
                    )
                )
    return sites


def _position_in_read(aligned_pairs: list, query_pos: int, read_length: int) -> float:
    """Calculate the fractional position of a variant within the read (0–1)."""
    if read_length == 0:
        return 0.5
    return query_pos / read_length


def extract_candidates_at_site(
    alignment_file: pysam.AlignmentFile,
    site: CompendiumSite,
    min_mapq: int = 0,
    min_bq: int = 0,
) -> Iterator[CandidateRead]:
    """Extract all reads overlapping a compendium site with per-read features.

    No VAF filtering is applied — even a single supporting read is informative
    at the tumor fractions relevant to MRD (10⁻⁵). No soft-clipping or masking
    at the variant position.

    Args:
        alignment_file: Open pysam AlignmentFile (CRAM or BAM).
        site: The compendium SNV site to query.
        min_mapq: Minimum mapping quality (default 0, filtering deferred to SVM).
        min_bq: Minimum base quality (default 0, filtering deferred to SVM).

    Yields:
        CandidateRead for each read covering the site.
    """
    for pileup_column in alignment_file.pileup(
        contig=site.chrom,
        start=site.pos,
        stop=site.pos + 1,
        truncate=True,
        min_mapping_quality=min_mapq,
        min_base_quality=min_bq,
        stepper="nofilter",
    ):
        for pileup_read in pileup_column.pileups:
            if pileup_read.is_del or pileup_read.is_refskip:
                continue

            alignment = pileup_read.alignment
            query_pos = pileup_read.query_position
            if query_pos is None:
                continue

            read_base = alignment.query_sequence[query_pos]
            qualities = alignment.query_qualities
            vbq = int(qualities[query_pos]) if qualities is not None else 0
            mrbq = float(np.mean(qualities)) if qualities is not None else 0.0
            pir = _position_in_read(
                alignment.get_aligned_pairs(),
                query_pos,
                alignment.query_length or 0,
            )

            yield CandidateRead(
                chrom=site.chrom,
                pos=site.pos,
                ref=site.ref,
                alt=site.alt,
                read_name=alignment.query_name or "",
                read_base=read_base,
                vbq=vbq,
                mrbq=mrbq,
                pir=pir,
                mapq=alignment.mapping_quality,
                is_read1=alignment.is_read1,
                mate_cigar=alignment.get_tag("MC") if alignment.has_tag("MC") else None,
                is_concordant=None,  # resolved in paired-end concordance step
                is_proper_pair=alignment.is_proper_pair,
                insert_size=abs(alignment.template_length),
            )


def extract_all_candidates(
    cram_path: str | Path,
    compendium_path: str | Path,
    reference_path: str | Path,
    min_mapq: int = 0,
    min_bq: int = 0,
) -> dict[CompendiumSite, list[CandidateRead]]:
    """Extract candidate reads at all compendium sites from a plasma CRAM.

    Args:
        cram_path: Path to plasma CRAM/BAM file.
        compendium_path: Path to patient SNV compendium VCF.
        reference_path: Path to reference FASTA (required for CRAM).
        min_mapq: Minimum mapping quality filter.
        min_bq: Minimum base quality filter.

    Returns:
        Dict mapping each CompendiumSite to its list of CandidateReads.
    """
    sites = load_compendium(compendium_path)
    results: dict[CompendiumSite, list[CandidateRead]] = {}

    with pysam.AlignmentFile(
        str(cram_path), reference_filename=str(reference_path)
    ) as aln:
        for site in sites:
            candidates = list(
                extract_candidates_at_site(aln, site, min_mapq=min_mapq, min_bq=min_bq)
            )
            results[site] = candidates

    return results


def resolve_paired_end_concordance(
    candidates: dict[CompendiumSite, list[CandidateRead]],
) -> dict[CompendiumSite, list[CandidateRead]]:
    """Resolve R1/R2 paired-end concordance for candidate reads.

    For cfDNA fragments with ~165 bp insert size and 150 bp PE reads, R1 and R2
    overlap substantially. If both reads of a pair support the same variant at
    the same position, this provides strong evidence against sequencing error
    (which would affect only one read).

    Args:
        candidates: Output from extract_all_candidates.

    Returns:
        Updated dict with is_concordant field resolved.
    """
    resolved: dict[CompendiumSite, list[CandidateRead]] = {}

    for site, reads in candidates.items():
        # Group by read name to find pairs
        by_name: dict[str, list[CandidateRead]] = {}
        for read in reads:
            by_name.setdefault(read.read_name, []).append(read)

        updated_reads: list[CandidateRead] = []
        for _name, pair in by_name.items():
            if len(pair) == 2:
                # Both mates cover the site
                both_support = all(r.read_base == site.alt for r in pair)
                both_ref = all(r.read_base == site.ref for r in pair)
                concordant = both_support or both_ref
                for r in pair:
                    updated_reads.append(
                        dataclasses.replace(r, is_concordant=concordant)
                    )
            else:
                # Only one mate covers the site — concordance unknown
                for r in pair:
                    updated_reads.append(
                        dataclasses.replace(r, is_concordant=None)
                    )

        resolved[site] = updated_reads

    return resolved
