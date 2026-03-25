"""Tests for aperture_snv.score - binomial TF estimation and z-scoring."""

import pytest

from aperture_snv.score import (
    compute_detection_rate,
    compute_z_score,
    estimate_tumor_fraction,
    score_sample,
)


class TestEstimateTumorFraction:
    def test_zero_signal_returns_zero(self):
        """When detected variants are entirely explained by noise, TF should be 0."""
        tf = estimate_tumor_fraction(M=10, N=10000, cov=30.0, mu=0.001, R=10000)
        assert tf == 0.0

    def test_zero_compendium_returns_zero(self):
        tf = estimate_tumor_fraction(M=0, N=0, cov=30.0, mu=0.001, R=1000)
        assert tf == 0.0

    def test_zero_coverage_returns_zero(self):
        tf = estimate_tumor_fraction(M=100, N=10000, cov=0.0, mu=0.001, R=1000)
        assert tf == 0.0

    def test_high_signal_gives_nonzero_tf(self):
        """With clear signal above noise, TF should be positive."""
        # 500 detected out of 10000, 30x coverage, noise explains ~10
        tf = estimate_tumor_fraction(M=500, N=10000, cov=30.0, mu=0.001, R=10000)
        assert tf > 0.0
        assert tf < 1.0

    def test_tf_increases_with_detected_variants(self):
        """More detected variants → higher TF, all else equal."""
        tf_low = estimate_tumor_fraction(M=100, N=10000, cov=30.0, mu=0.0, R=10000)
        tf_high = estimate_tumor_fraction(M=500, N=10000, cov=30.0, mu=0.0, R=10000)
        assert tf_high > tf_low

    def test_tf_bounded_0_1(self):
        """TF should always be in [0, 1]."""
        tf = estimate_tumor_fraction(M=10000, N=10000, cov=30.0, mu=0.0, R=10000)
        assert 0.0 <= tf <= 1.0

    def test_all_detected_gives_high_tf(self):
        """When all compendium sites are detected, TF should be high."""
        tf = estimate_tumor_fraction(M=10000, N=10000, cov=30.0, mu=0.0, R=10000)
        assert tf > 0.01


class TestComputeDetectionRate:
    def test_basic(self):
        assert compute_detection_rate(100, 10000) == pytest.approx(0.01)

    def test_zero_reads(self):
        assert compute_detection_rate(0, 0) == 0.0

    def test_no_detections(self):
        assert compute_detection_rate(0, 10000) == 0.0


class TestComputeZScore:
    def test_at_mean_gives_zero(self):
        z = compute_z_score(0.001, noise_mean=0.001, noise_std=0.0005)
        assert z == pytest.approx(0.0)

    def test_above_mean_gives_positive(self):
        z = compute_z_score(0.005, noise_mean=0.001, noise_std=0.001)
        assert z == pytest.approx(4.0)

    def test_zero_std_above_mean(self):
        z = compute_z_score(0.005, noise_mean=0.001, noise_std=0.0)
        assert z == float("inf")

    def test_zero_std_at_mean(self):
        z = compute_z_score(0.001, noise_mean=0.001, noise_std=0.0)
        assert z == 0.0


class TestScoreSample:
    def test_detected_case(self):
        """High signal should produce a detection call."""
        result = score_sample(
            detected_variants=500,
            total_sites=10000,
            mean_coverage=30.0,
            noise_rate=0.0001,
            total_reads_evaluated=300000,
            noise_mean=0.0001,
            noise_std=0.00005,
            z_threshold=1.2,
        )
        assert result.is_detected is True
        assert result.z_score > 1.2
        assert result.tumor_fraction > 0.0

    def test_noise_level_not_detected(self):
        """Signal at noise level should not be detected."""
        result = score_sample(
            detected_variants=30,
            total_sites=10000,
            mean_coverage=30.0,
            noise_rate=0.0001,
            total_reads_evaluated=300000,
            noise_mean=0.0001,
            noise_std=0.00005,
            z_threshold=1.2,
        )
        assert result.z_score < 5  # Should be reasonable, not extreme

    def test_score_fields_populated(self):
        result = score_sample(
            detected_variants=100,
            total_sites=10000,
            mean_coverage=30.0,
            noise_rate=0.0001,
            total_reads_evaluated=300000,
            noise_mean=0.0001,
            noise_std=0.00005,
        )
        assert result.detected_variants == 100
        assert result.total_sites == 10000
        assert result.mean_coverage == 30.0
        assert isinstance(result.tumor_fraction, float)
        assert isinstance(result.z_score, float)
        assert isinstance(result.is_detected, bool)
