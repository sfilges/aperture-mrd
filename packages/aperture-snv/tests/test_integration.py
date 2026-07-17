"""Integration tests exercising the real pysam pileup path against a tiny BAM.

The fixture (`tests/test_data/mini.bam`) is a slice of a real cfDNA WGS sample (~10x,
ichorCNA TF ~45%) restricted to windows around the compendium SNV sites, with headers
anonymized. Unlike the pure `build_fragment` unit tests, these run the actual pileup,
so they catch bugs in the pysam layer — e.g. pileup overlap handling silently mangling
base qualities.
"""

from pathlib import Path

import pytest

from aperture_snv.extract import extract_all_fragments

DATA = Path(__file__).parent / "test_data"
BAM = DATA / "mini.bam"
COMPENDIUM = DATA / "mini_compendium.vcf.gz"

pytestmark = pytest.mark.skipif(
    not (BAM.exists() and COMPENDIUM.exists()),
    reason="integration fixture (tests/test_data/mini.bam) not present",
)


@pytest.fixture(scope="module")
def fragments():
    result = extract_all_fragments(alignment_file_path=BAM, compendium_vcf_path=COMPENDIUM)
    return [f for lst in result.values() for f in lst]


def test_extraction_produces_fragments(fragments):
    assert len(fragments) > 0


def test_base_quality_is_not_overlap_summed(fragments):
    """Regression guard: pileup ignore_overlaps must be off.

    With pysam's default ignore_overlaps=True, overlapping mate pairs are collapsed
    by summing base qualities (e.g. 40+40=80), producing impossible Phred values.
    """
    assert max(f.vbq for f in fragments) <= 41


def test_paired_mates_are_grouped_and_concordant(fragments):
    """Overlapping R1/R2 at a site should collapse into one fragment with 2 mates."""
    paired = [f for f in fragments if f.n_mates_at_site == 2]
    assert paired, "expected some fragments covered by both mates"
    assert all(f.concordance in {"concordant", "discordant"} for f in paired)
    assert any(f.concordance == "concordant" for f in paired)


def test_alt_supporting_signal_present(fragments):
    """TF ~45% tumor-informed sites should yield alt-supporting fragments."""
    alt = [f for f in fragments if f.supports_alt]
    assert alt, "expected alt-supporting fragments at TF ~45%"
    assert all(f.read_base == f.alt for f in alt)


def test_concordance_values_are_valid(fragments):
    assert all(f.concordance in {"concordant", "discordant", "single"} for f in fragments)
    assert all(f.n_mates_at_site == (1 if f.concordance == "single" else 2) for f in fragments)
