"""
End-to-end test: create binary phenotypes from existing test data,
estimate h2 with AdjHE using liability-scale transformation,
and verify estimates are reasonable.
"""
import sys
sys.path.insert(0, 'src')

import numpy as np
import pandas as pd
import tempfile
import os

from Estimate.estimators.all_estimators import h2Estimation, load_n_estimate
from Estimate.data_input.load_data import load_everything


def test_binary_pheno_h2_estimation():
    """
    Load existing EUR2 test data, create binary phenotype from
    continuous pheno_1, estimate h2 with AdjHE using liability-scale
    transformation, and verify the estimate is in a reasonable range.
    
    This test verifies the end-to-end binary phenotype pipeline:
    1. Create binary phenotype from continuous source
    2. Load data and apply liability-scale transformation
    3. Estimate h2 using AdjHE
    4. Verify estimate is positive and within [0, 1]
    """
    test_dir = "tests/test_data"
    pheno_file = os.path.join(test_dir, "EUR_simulation2.pheno")
    grm_prefix = os.path.join(test_dir, "EUR2")
    
    # Read original phenotype file
    pheno_df = pd.read_csv(pheno_file, sep="\t")
    
    # Create binary phenotype from pheno_1 (top ~10% as cases)
    source_pheno = "pheno_1"
    cont_values = pheno_df[source_pheno].values
    threshold = np.percentile(cont_values, 90)
    binary_values = (cont_values > threshold).astype(int)
    binary_pheno_name = f"{source_pheno}_binary"
    
    # Create temp phenotype file with binary phenotype
    tmpdir = tempfile.mkdtemp()
    binary_pheno_file = os.path.join(tmpdir, "binary_pheno.pheno")
    pheno_df[binary_pheno_name] = binary_values
    pheno_df[["FID", "IID", binary_pheno_name]].to_csv(
        binary_pheno_file, sep="\t", index=False
    )
    
    # Load data with binary phenotype
    config = {
        "prefix": grm_prefix,
        "pheno": binary_pheno_file,
        "npc": [0],
        "continuousPhenos": [],
        "binPhenos": [binary_pheno_name],
        "prevalence": {binary_pheno_name: float(np.mean(binary_values))},
        "preprocess": "None",
        "Method": "AdjHE",
        "k": 0,
        "std": False,
        "qcovar": None,
        "covar_discrete": None,
        "random_groups": "None",
        "Naive": False,
        "loop_covars": False,
    }
    
    df, GRM, phenotypes, ids = load_everything(args=config)
    
    # Create h2Estimation object
    ests = h2Estimation(args=config)
    ests.df = df
    ests.GRM = GRM
    ests.continuousPhenos = []
    ests.binPhenos = [binary_pheno_name]
    ests.phenotypes = [binary_pheno_name]
    
    # Apply liability-scale transformation to binary phenotypes
    ests._apply_liability_scale(ests.binPhenos)
    
    # Estimate h2 for binary phenotype
    r = load_n_estimate(
        df=ests.df, nnpc=0, mp=binary_pheno_name, GRM=ests.GRM,
        std=True, Method="AdjHE", random_groups=None,
        homo=True, PC_effect="mixed", qcovar=None,
        covar_discrete=None, all_cols=None
    )
    
    h2_est = r["h2"].values[0]
    
    # Cleanup
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)
    
    # Verify: h2 should be positive and within [0, 1]
    assert h2_est >= 0, f"Negative h2 estimate: {h2_est}"
    assert h2_est <= 1.0, f"h2 > 1: {h2_est}"
    assert h2_est > 0.05, f"h2 estimate too low: {h2_est}"
