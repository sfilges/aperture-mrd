"""Tests for aperture_snv.filters — pure fragment-stage filter predicates."""

from types import SimpleNamespace

from aperture_snv.extract import CompendiumSite, build_fragment
from aperture_snv.filters import (
    DEFAULT_FILTERS,
    FilterConfig,
    fragment_passes,
    read_flags_ok,
)

from .test_extract import _mate  # reuse the MateObservation builder

SITE = CompendiumSite(chrom="chr1", pos=1000, ref="C", alt="T")


def _aln(**flags):
    defaults = dict(
        is_unmapped=False,
        is_secondary=False,
        is_supplementary=False,
        is_duplicate=False,
        is_qcfail=False,
    )
    defaults.update(flags)
    return SimpleNamespace(**defaults)


class TestReadFlagsOk:
    def test_primary_mapped_read_passes(self):
        assert read_flags_ok(_aln()) is True

    def test_rejects_each_bad_flag(self):
        for flag in ("is_unmapped", "is_secondary", "is_supplementary", "is_duplicate", "is_qcfail"):
            assert read_flags_ok(_aln(**{flag: True})) is False, flag


class TestFragmentPasses:
    def _fragment(self, insert):
        return build_fragment(SITE, [_mate(template_length=insert)])

    def test_within_range_passes(self):
        assert fragment_passes(self._fragment(165), DEFAULT_FILTERS) is True

    def test_too_short_fails(self):
        assert fragment_passes(self._fragment(30), DEFAULT_FILTERS) is False

    def test_too_long_fails(self):
        assert fragment_passes(self._fragment(500), DEFAULT_FILTERS) is False

    def test_boundaries_inclusive(self):
        cfg = FilterConfig(min_insert=40, max_insert=240)
        assert fragment_passes(self._fragment(40), cfg) is True
        assert fragment_passes(self._fragment(240), cfg) is True
