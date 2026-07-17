"""Tests for aperture_snv.extract — fragment assembly logic.

These cover the pure-Python fragment collapse (concordance, consensus base,
3'-distance, fragment length) without requiring BAM/CRAM files. Integration tests
that exercise the pysam pileup path require fixtures and are run separately.
"""

import pytest

from aperture_snv.extract import (
    CompendiumSite,
    MateObservation,
    _dist_from_3prime,
    _fragment_length,
    build_fragment,
    load_compendium,
)

SITE = CompendiumSite(chrom="chr1", pos=1000, ref="C", alt="T")

_VCF_HEADER = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class TestLoadCompendium:
    def test_dedups_and_skips_mnv(self, tmp_path):
        vcf = tmp_path / "c.vcf"
        vcf.write_text(
            _VCF_HEADER
            + "chr1\t100\t.\tC\tT\t.\tPASS\t.\n"  # SNV
            + "chr1\t100\t.\tC\tT\t.\tPASS\t.\n"  # exact duplicate (annotation) -> collapsed
            + "chr1\t200\t.\tG\tA\t.\tPASS\t.\n"  # distinct SNV
            + "chr1\t300\t.\tCC\tGG\t.\tPASS\t.\n"  # MNV -> skipped
        )
        sites = load_compendium(vcf)
        assert len(sites) == 2
        assert len(set(sites)) == 2
        assert CompendiumSite("chr1", 99, "C", "T") in sites  # pysam pos is 0-based


def _mate(
    read_name: str = "frag_1",
    base_at_site: str = "T",
    is_read1: bool = True,
    is_reverse: bool = False,
    bq_at_site: int = 35,
    query_position: int = 45,
    query_length: int = 150,
    mean_bq: float = 32.0,
    n_low_bq: int = 2,
    mapq: int = 60,
    edit_distance: int = 0,
    softclip_len: int = 0,
    is_proper_pair: bool = True,
    reference_start: int = 900,
    reference_end: int = 1050,
    template_length: int = 165,
) -> MateObservation:
    return MateObservation(
        read_name=read_name,
        is_read1=is_read1,
        is_reverse=is_reverse,
        base_at_site=base_at_site,
        bq_at_site=bq_at_site,
        query_position=query_position,
        query_length=query_length,
        mean_bq=mean_bq,
        n_low_bq=n_low_bq,
        mapq=mapq,
        edit_distance=edit_distance,
        softclip_len=softclip_len,
        is_proper_pair=is_proper_pair,
        reference_start=reference_start,
        reference_end=reference_end,
        template_length=template_length,
    )


class TestConcordance:
    def test_both_support_alt_is_concordant(self):
        frag = build_fragment(SITE, [_mate(base_at_site="T"), _mate(base_at_site="T")])
        assert frag.concordance == "concordant"
        assert frag.supports_alt is True
        assert frag.n_mates_at_site == 2

    def test_both_support_ref_is_concordant(self):
        frag = build_fragment(SITE, [_mate(base_at_site="C"), _mate(base_at_site="C")])
        assert frag.concordance == "concordant"
        assert frag.supports_alt is False

    def test_disagreeing_mates_are_discordant(self):
        frag = build_fragment(SITE, [_mate(base_at_site="T"), _mate(base_at_site="C")])
        assert frag.concordance == "discordant"

    def test_single_mate(self):
        frag = build_fragment(SITE, [_mate(base_at_site="T")])
        assert frag.concordance == "single"
        assert frag.n_mates_at_site == 1


class TestConsensusBase:
    def test_higher_quality_mate_wins_discordant_call(self):
        alt = _mate(base_at_site="T", bq_at_site=38)
        ref = _mate(base_at_site="C", bq_at_site=12)
        frag = build_fragment(SITE, [ref, alt])
        assert frag.read_base == "T"
        assert frag.supports_alt is True
        # vbq is the best base-quality evidence at the site.
        assert frag.vbq == 38


class TestEditDistanceAggregation:
    def test_takes_worst_mate(self):
        frag = build_fragment(SITE, [_mate(edit_distance=0), _mate(edit_distance=3)])
        assert frag.edit_distance == 3


class TestDistFrom3Prime:
    def test_forward_read(self):
        # forward: 3' end is the high query coordinate
        m = _mate(is_reverse=False, query_position=45, query_length=150)
        assert _dist_from_3prime(m) == 150 - 1 - 45

    def test_reverse_read(self):
        # reverse: query position 0 is the 3' end
        m = _mate(is_reverse=True, query_position=45, query_length=150)
        assert _dist_from_3prime(m) == 45

    def test_min_across_mates_used(self):
        fwd = _mate(is_reverse=False, query_position=140, query_length=150)  # 9 from 3'
        rev = _mate(is_reverse=True, query_position=45, query_length=150)  # 45 from 3'
        frag = build_fragment(SITE, [fwd, rev])
        assert frag.dist_3prime == 9


class TestFragmentLength:
    def test_prefers_template_length(self):
        assert _fragment_length([_mate(template_length=165)]) == 165

    def test_falls_back_to_reference_span(self):
        m = _mate(template_length=0, reference_start=900, reference_end=1060)
        assert _fragment_length([m]) == 160


def test_build_fragment_requires_a_mate():
    with pytest.raises(ValueError, match="at least one mate"):
        build_fragment(SITE, [])
