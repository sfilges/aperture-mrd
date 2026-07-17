"""Tests for aperture_snv.score — PPM mutational load and detection z-scoring."""

import pytest

from aperture_snv.extract import build_fragment
from aperture_snv.score import (
    compute_ppm,
    compute_z_score,
    score_sample,
    sum_supporting_reads,
)

from .test_extract import SITE, _mate


def _alt_fragment(base="T"):
    return build_fragment(SITE, [_mate(base_at_site=base)])


class TestSumSupportingReads:
    def test_counts_only_alt_supporting(self):
        frags = [_alt_fragment("T"), _alt_fragment("C"), _alt_fragment("T")]
        assert sum_supporting_reads(frags) == 2.0

    def test_default_weight_is_one_per_fragment(self):
        assert sum_supporting_reads([_alt_fragment("T")]) == 1.0

    def test_weight_fn_applied_to_alt_only(self):
        frags = [_alt_fragment("T"), _alt_fragment("C")]  # one alt, one ref
        assert sum_supporting_reads(frags, weight_fn=lambda _f: 0.25) == 0.25

    def test_empty(self):
        assert sum_supporting_reads([]) == 0.0


class TestComputePPM:
    def test_basic(self):
        # 10 supporting reads out of 1e6 aligned → 10 PPM
        assert compute_ppm(10, 1_000_000) == pytest.approx(10.0)

    def test_fractional_supporting_reads(self):
        assert compute_ppm(2.5, 1_000_000) == pytest.approx(2.5)

    def test_zero_aligned_reads(self):
        assert compute_ppm(10, 0) == 0.0

    def test_scales_with_depth(self):
        assert compute_ppm(10, 2_000_000) == pytest.approx(5.0)


class TestComputeZScore:
    def test_at_mean_is_zero(self):
        assert compute_z_score(1.0, noise_mean=1.0, noise_std=0.5) == pytest.approx(0.0)

    def test_above_mean_positive(self):
        assert compute_z_score(5.0, noise_mean=1.0, noise_std=1.0) == pytest.approx(4.0)

    def test_zero_std_above_mean(self):
        assert compute_z_score(5.0, noise_mean=1.0, noise_std=0.0) == float("inf")

    def test_zero_std_at_mean(self):
        assert compute_z_score(1.0, noise_mean=1.0, noise_std=0.0) == 0.0


class TestScoreSample:
    def test_detected_case(self):
        result = score_sample(
            supporting_reads=500,
            total_aligned_reads=1_000_000,
            noise_mean=1.0,
            noise_std=1.0,
            z_threshold=1.2,
        )
        assert result.ppm == pytest.approx(500.0)
        assert result.is_detected is True
        assert result.z_score > 1.2

    def test_noise_level_not_detected(self):
        result = score_sample(
            supporting_reads=1,
            total_aligned_reads=1_000_000,
            noise_mean=1.0,
            noise_std=0.5,
            z_threshold=1.2,
        )
        assert result.is_detected is False

    def test_fields_populated(self):
        result = score_sample(
            supporting_reads=100,
            total_aligned_reads=500_000,
            noise_mean=1.0,
            noise_std=0.5,
        )
        assert result.supporting_reads == 100
        assert result.total_aligned_reads == 500_000
        assert isinstance(result.ppm, float)
        assert isinstance(result.z_score, float)
        assert isinstance(result.is_detected, bool)
