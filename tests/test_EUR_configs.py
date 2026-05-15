#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for EUR simulated data config files.

The EUR simulation data was generated with heritability ~0.4.
We test that >50% of estimates fall within a reasonable range (0.2 to 0.8).
"""

import pytest
import pandas as pd
import numpy as np
from Estimate.data_input.parser import read_flags
from Estimate.estimators.all_estimators import h2Estimation

# Tolerance for heritability estimates
H2_MIN = 0.2
H2_MAX = 0.8
MIN_PASS_RATIO = 0.5  # At least 50% of phenotypes should pass


@pytest.fixture(scope="module")
def config_dir():
    return "tests/test_data"


@pytest.fixture
def config_file():
    """Default config file for testing."""
    return "config_AdjHE.json"


@pytest.fixture
def args(config_file, config_dir):
    """Load config file and return args dictionary."""
    return read_flags({"argfile": f"{config_dir}/{config_file}"})


@pytest.fixture
def estimator(args):
    """Create h2Estimation object from args."""
    return h2Estimation(args=args)


def run_config(config_file, config_dir="tests/test_data"):
    """Helper to run MASH with a config file and return results."""
    args = read_flags({"argfile": f"{config_dir}/{config_file}"})
    est = h2Estimation(args=args)
    results = est.estimate()
    return results


def check_results_numeric(results, config_name):
    """Check that h2 estimates are valid numeric values (no NaN)."""
    total = len(results)
    if total == 0:
        pytest.fail(f"{config_name}: No results produced")
    
    valid = sum(1 for _, row in results.iterrows()
                if not pd.isna(row["h2"]) and not np.isinf(row["h2"]))
    
    print(f"\n{config_name}: {valid}/{total} valid h2 estimates")
    print(f"  h2 range: [{results['h2'].min():.4f}, {results['h2'].max():.4f}]")
    print(f"  Mean h2: {results['h2'].mean():.4f}")
    
    assert valid > 0, \
        f"{config_name}: No valid h2 estimates produced (all NaN or Inf)"


def check_results(results, config_name):
    """Check that >50% of h2 estimates are within [H2_MIN, H2_MAX]."""
    total = len(results)
    if total == 0:
        pytest.fail(f"{config_name}: No results produced")
    
    within_range = sum(1 for _, row in results.iterrows() 
                     if H2_MIN <= row["h2"] <= H2_MAX)
    ratio = within_range / total
    
    # Print summary for debugging
    print(f"\n{config_name}: {within_range}/{total} ({ratio:.1%}) phenotypes within [{H2_MIN}, {H2_MAX}]")
    print(f"  h2 range: [{results['h2'].min():.4f}, {results['h2'].max():.4f}]")
    print(f"  Mean h2: {results['h2'].mean():.4f}")
    
    assert ratio >= MIN_PASS_RATIO, \
        f"{config_name}: Only {ratio:.1%} ({within_range}/{total}) phenotypes within [{H2_MIN}, {H2_MAX}], need {MIN_PASS_RATIO:.0%}"


class TestAdjHE:
    """Tests for AdjHE method configs."""

    @pytest.mark.AdjHE
    def test_AdjHE_basic(self):
        """Test basic AdjHE with covariates."""
        results = run_config("config_AdjHE.json")
        check_results(results, "AdjHE basic")

    @pytest.mark.AdjHE
    def test_AdjHE_RV(self):
        """Test AdjHE with random_groups (site) and qcovar/covar_discrete."""
        results = run_config("config_AdjHE_RV.json")
        check_results(results, "AdjHE_RV")

    @pytest.mark.AdjHE
    def test_AdjHE_complex(self):
        """Test AdjHE with multiple phenotype files, #FID eigenvec, and participant_id IID."""
        results = run_config("config_complexAdjHE.json")
        check_results_numeric(results, "AdjHE complex")


class TestGCTA:
    """Tests for GCTA method configs."""

    @pytest.mark.GCTA
    def test_GCTA_basic(self):
        """Test GCTA with basic covariates."""
        results = run_config("config_GCTA.json")
        check_results(results, "GCTA basic")

    @pytest.mark.GCTA
    def test_GCTA_named(self):
        """Test GCTA with named covariates."""
        results = run_config("config_GCTA_named.json")
        check_results(results, "GCTA named")

    @pytest.mark.GCTA
    def test_GCTA_csv(self):
        """Test GCTA with CSV phenotype file."""
        results = run_config("config_GCTA_csv.json")
        check_results(results, "GCTA csv")

    @pytest.mark.GCTA
    def test_GCTA_discrete(self):
        """Test GCTA with discrete covariates."""
        results = run_config("config_GCTA_discrete.json")
        check_results(results, "GCTA discrete")

    @pytest.mark.GCTA
    def test_GCTA_fixed_mixed(self):
        """Test GCTA with fixed mixed effects."""
        results = run_config("config_GCTA_fixed_mixed.json")
        check_results(results, "GCTA fixed_mixed")

    @pytest.mark.GCTA
    def test_GCTA_mixed(self):
        """Test GCTA with mixed model."""
        results = run_config("config_GCTA_mixed.json")
        check_results(results, "GCTA mixed")

    @pytest.mark.GCTA
    def test_GCTA_mixed_both(self):
        """Test GCTA with mixed both covariate file."""
        results = run_config("config_GCTA_mixed_both.json")
        check_results(results, "GCTA mixed_both")

    @pytest.mark.GCTA
    def test_GCTA_multi_covar(self):
        """Test GCTA with multiple covariate files."""
        results = run_config("config_GCTA_multi_covar.json")
        check_results(results, "GCTA multi_covar")


class TestSWD:
    """Tests for SWD method configs."""

    @pytest.mark.SWD
    def test_SWD_basic(self):
        """Test SWD with site as random variable."""
        results = run_config("config_SWD.json")
        check_results(results, "SWD")


class TestCovbat:
    """Tests for Covbat (Combat) method configs."""

    @pytest.mark.Combat
    def test_Covbat_basic(self):
        """Test Covbat with site as random variable."""
        results = run_config("config_COMBAT.json")
        check_results(results, "Covbat")


class TestAllConfigs:
    """Meta-test: run all configs and summarize results."""

    def test_all_configs_run(self):
        """Test that all configs can run without errors."""
        configs = [
            "config_AdjHE.json",
            "config_AdjHE_RV.json",
            "config_GCTA.json",
            "config_GCTA_named.json",
            "config_GCTA_csv.json",
            "config_GCTA_discrete.json",
            "config_GCTA_fixed_mixed.json",
            "config_GCTA_mixed.json",
            "config_GCTA_mixed_both.json",
            "config_GCTA_multi_covar.json",
            "config_complexAdjHE.json",
            "config_COMBAT.json",
        ]
        for config in configs:
            try:
                results = run_config(config)
                assert results is not None and len(results) > 0, \
                    f"{config} produced no results"
            except Exception as e:
                pytest.fail(f"{config} failed to run: {e}")
