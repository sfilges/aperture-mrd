"""Tests for aperture_snv.extract — candidate extraction logic.

These tests cover the pure-Python logic (dataclasses, concordance resolution)
without requiring actual BAM/CRAM files. Integration tests with pysam
require test fixtures and are run separately.
"""

import pytest

from aperture_snv.extract import (
    CandidateRead,
    CompendiumSite,
    _position_in_read,
    resolve_paired_end_concordance,
)


def _make_read(
    read_name: str = "read_1",
    read_base: str = "T",
    is_read1: bool = True,
    **kwargs,
) -> CandidateRead:
    defaults = dict(
        chrom="chr1",
        pos=1000,
        ref="C",
        alt="T",
        vbq=30,
        mrbq=28.5,
        pir=0.3,
        mapq=60,
        mate_cigar=None,
        is_concordant=None,
        is_proper_pair=True,
        insert_size=165,
    )
    defaults.update(kwargs)
    return CandidateRead(
        read_name=read_name,
        read_base=read_base,
        is_read1=is_read1,
        **defaults,
    )


SITE = CompendiumSite(chrom="chr1", pos=1000, ref="C", alt="T")


class TestPositionInRead:
    def test_middle(self):
        assert _position_in_read([], 75, 150) == pytest.approx(0.5)

    def test_start(self):
        assert _position_in_read([], 0, 150) == pytest.approx(0.0)

    def test_zero_length(self):
        assert _position_in_read([], 0, 0) == 0.5


class TestResolvepairedEndConcordance:
    def test_both_support_alt(self):
        """Both mates support alt → concordant."""
        r1 = _make_read(read_name="pair1", read_base="T", is_read1=True)
        r2 = _make_read(read_name="pair1", read_base="T", is_read1=False)
        result = resolve_paired_end_concordance({SITE: [r1, r2]})
        assert all(r.is_concordant is True for r in result[SITE])

    def test_both_support_ref(self):
        """Both mates support ref → concordant (both agree, just not variant)."""
        r1 = _make_read(read_name="pair1", read_base="C", is_read1=True)
        r2 = _make_read(read_name="pair1", read_base="C", is_read1=False)
        result = resolve_paired_end_concordance({SITE: [r1, r2]})
        assert all(r.is_concordant is True for r in result[SITE])

    def test_discordant(self):
        """One mate supports alt, other supports ref → discordant."""
        r1 = _make_read(read_name="pair1", read_base="T", is_read1=True)
        r2 = _make_read(read_name="pair1", read_base="C", is_read1=False)
        result = resolve_paired_end_concordance({SITE: [r1, r2]})
        assert all(r.is_concordant is False for r in result[SITE])

    def test_single_read_unknown(self):
        """Single read (no mate at site) → concordance unknown."""
        r1 = _make_read(read_name="solo", read_base="T", is_read1=True)
        result = resolve_paired_end_concordance({SITE: [r1]})
        assert result[SITE][0].is_concordant is None

    def test_multiple_pairs(self):
        """Multiple read pairs are resolved independently."""
        reads = [
            _make_read(read_name="pair1", read_base="T", is_read1=True),
            _make_read(read_name="pair1", read_base="T", is_read1=False),
            _make_read(read_name="pair2", read_base="T", is_read1=True),
            _make_read(read_name="pair2", read_base="C", is_read1=False),
        ]
        result = resolve_paired_end_concordance({SITE: reads})
        by_name = {}
        for r in result[SITE]:
            by_name.setdefault(r.read_name, []).append(r)

        # pair1: both alt → concordant
        assert all(r.is_concordant is True for r in by_name["pair1"])
        # pair2: discordant
        assert all(r.is_concordant is False for r in by_name["pair2"])

    def test_empty_input(self):
        result = resolve_paired_end_concordance({})
        assert result == {}

    def test_no_reads_at_site(self):
        result = resolve_paired_end_concordance({SITE: []})
        assert result[SITE] == []
