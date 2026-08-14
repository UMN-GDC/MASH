#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for data loading functionality, specifically PC/eigenvec file handling.
"""

import pytest
import pandas as pd
import numpy as np
import tempfile
import os

from Estimate.data_input.load_data import load_tables
from Estimate.data_input.parser import read_flags


def test_loading_eigenvec_hash_fid_uppercase_pcs():
    """Test loading eigenvec file with #FID header and uppercase PC columns."""
    # Create test eigenvec file with #FID and uppercase PC columns
    eigenvec_content = """#FID	IID	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
0	SAMPLE1	0.1	0.2	0.3	0.4	0.5	0.6	0.7	0.8	0.9	1.0
1	SAMPLE2	0.11	0.12	0.13	0.14	0.15	0.16	0.17	0.18	0.19	0.20
2	SAMPLE3	0.21	0.22	0.23	0.24	0.25	0.26	0.27	0.28	0.29	0.30
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    try:
        # Create minimal args for loading
        args = {
            "PC": eigenvec_file,
            "pheno": None,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID"
        }
        
        # Create dummy IDs dataframe (normally comes from GRM file)
        ids = pd.DataFrame({
            "FID": ["0", "0", "0"],
            "IID": ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        # Load tables
        df = load_tables(ids, args)
        
        # Check that FID and IID are present
        assert "FID" in df.columns
        assert "IID" in df.columns
        
        # Check that PC columns were lowercased
        expected_pc_cols = [f"pc{i}" for i in range(1, 11)]
        for pc_col in expected_pc_cols:
            assert pc_col in df.columns, f"Missing column {pc_col}"
            
        # Check that original uppercase PC columns are not present (should be lowercased)
        for i in range(1, 11):
            assert f"PC{i}" not in df.columns, f"Uppercase column PC{i} should have been lowercased"
            
        # Check values are correct
        assert df.loc[df["IID"] == "SAMPLE1", "pc1"].iloc[0] == 0.1
        assert df.loc[df["IID"] == "SAMPLE2", "pc5"].iloc[0] == 0.15
        assert df.loc[df["IID"] == "SAMPLE3", "pc10"].iloc[0] == 0.30
        
    finally:
        os.unlink(eigenvec_file)


def test_loading_eigenvec_hash_fid_with_n_a():
    """Test loading eigenvec file with #FID header and 'n/a' values treated as NA."""
    # Create test eigenvec file with n/a values
    eigenvec_content = """#FID	IID	PC1	PC2	PC3
0	SAMPLE1	0.1	n/a	0.3
1	SAMPLE2	n/a	0.2	0.4
2	SAMPLE3	0.3	0.4	n/a
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    try:
        args = {
            "PC": eigenvec_file,
            "pheno": None,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID"
        }
        
        ids = pd.DataFrame({
            "FID": ["0", "0", "0"],
            "IID": ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        df = load_tables(ids, args)
        
        # Check that n/a values were converted to NaN
        assert pd.isna(df.loc[df["IID"] == "SAMPLE1", "pc2"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "SAMPLE2", "pc1"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "SAMPLE3", "pc3"].iloc[0])
        
        # Check that valid values are preserved
        assert df.loc[df["IID"] == "SAMPLE1", "pc1"].iloc[0] == 0.1
        assert df.loc[df["IID"] == "SAMPLE2", "pc2"].iloc[0] == 0.2
        assert df.loc[df["IID"] == "SAMPLE3", "pc1"].iloc[0] == 0.3
        
    finally:
        os.unlink(eigenvec_file)


def test_loading_eigenvec_participant_id_custom_iid():
    """Test loading eigenvec file with participant_id column and custom iid_col."""
    # Create test eigenvec file with participant_id instead of IID
    eigenvec_content = """#FID	participant_id	PC1	PC2	PC3
0	SAMPLE1	0.1	0.2	0.3
1	SAMPLE2	0.4	0.5	0.6
2	SAMPLE3	0.7	0.8	0.9
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    try:
        args = {
            "PC": eigenvec_file,
            "pheno": None,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID",  # We want IID as the standard column name
        }
        
        ids = pd.DataFrame({
            "FID": ["0", "0", "0"],
            "IID": ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        df = load_tables(ids, args)
        
        # Check that participant_id was used to populate IID column
        assert "IID" in df.columns
        assert set(df["IID"].tolist()) == {"SAMPLE1", "SAMPLE2", "SAMPLE3"}
        
        # Check PC columns were lowercased
        assert "pc1" in df.columns
        assert "pc2" in df.columns
        assert "pc3" in df.columns
        
        # Check values are correct
        assert df.loc[df["IID"] == "SAMPLE1", "pc1"].iloc[0] == 0.1
        assert df.loc[df["IID"] == "SAMPLE2", "pc3"].iloc[0] == 0.6
        assert df.loc[df["IID"] == "SAMPLE3", "pc2"].iloc[0] == 0.8
        
    finally:
        os.unlink(eigenvec_file)


def test_loading_eigenvec_no_header_fallback():
    """Test that files without proper headers fall back to the original logic."""
    # Create test file WITHOUT header (just data)
    eigenvec_content = """0	SAMPLE1	0.1	0.2	0.3
1	SAMPLE2	0.4	0.5	0.6
2	SAMPLE3	0.7	0.8	0.9
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    try:
        args = {
            "PC": eigenvec_file,
            "pheno": None,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID"
        }
        
        ids = pd.DataFrame({
            "FID": ["0", "0", "0"],
            "IID": ["SAMPLE1", "SAMPLE2", "SAMPLE3"]
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        df = load_tables(ids, args)
        
        # Should have FID, IID, and PC columns (pc1, pc2, pc3)
        assert "FID" in df.columns
        assert "IID" in df.columns
        assert "pc1" in df.columns
        assert "pc2" in df.columns
        assert "pc3" in df.columns
        
        # Check values - first data row should be used as data, not header
        assert df.loc[df["IID"] == "SAMPLE1", "pc1"].iloc[0] == 0.1
        assert df.loc[df["IID"] == "SAMPLE2", "pc2"].iloc[0] == 0.5
        assert df.loc[df["IID"] == "SAMPLE3", "pc3"].iloc[0] == 0.9
        
    finally:
        os.unlink(eigenvec_file)


def test_covar_filter_discrete():
    """Test that covar_filter filters correctly on a discrete covariate column."""
    # Create eigenvec file with 6 samples
    eigenvec_content = "#FID	IID	PC1	PC2\n0	S1	0.1	0.2\n0	S2	0.3	0.4\n0	S3	0.5	0.6\n0	S4	0.7	0.8\n0	S5	0.9	1.0\n0	S6	1.1	1.2\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    # Create covariate file with dummy_sex column (M/F)
    cov_content = "participant_id\tdummy_sex\nS1\tM\nS2\tF\nS3\tM\nS4\tM\nS5\tF\nS6\tM\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write(cov_content)
        cov_file = f.name
    
    try:
        args = {
            "PC": eigenvec_file,
            "pheno": None,
            "covar": [cov_file],
            "covar_filter": "dummy_sex==M",
            "fid_col": "FID",
            "iid_col": "IID",
        }
        
        ids = pd.DataFrame({
            "FID": ["0"] * 6,
            "IID": ["S1", "S2", "S3", "S4", "S5", "S6"],
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        df = load_tables(ids, args)
        
        assert "dummy_sex" in df.columns, "dummy_sex column missing after load_tables"
        present = df["dummy_sex"].dropna()
        assert (present == "M").all(), f"Expected all present dummy_sex values to be 'M', got {present.unique()}"
        # S2 and S5 were F, so they should be NaN (filtered then left-joined)
        assert pd.isna(df.loc[df["IID"] == "S2", "dummy_sex"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "S5", "dummy_sex"].iloc[0])
        
    finally:
        os.unlink(eigenvec_file)
        os.unlink(cov_file)


def test_pheno_filter_numeric():
    """Test that pheno_filter filters correctly on a numeric column."""
    eigenvec_content = "#FID	IID	PC1	PC2\n0	S1	0.1	0.2\n0	S2	0.3	0.4\n0	S3	0.5	0.6\n0	S4	0.7	0.8\n0	S5	0.9	1.0\n0	S6	1.1	1.2\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.eigenvec', delete=False) as f:
        f.write(eigenvec_content)
        eigenvec_file = f.name
    
    # Create phenotype file with both FID and IID
    pheno_content = "FID\tIID\tpheno_1\n0\tS1\t1\n0\tS2\t5\n0\tS3\t2\n0\tS4\t3\n0\tS5\t-1\n0\tS6\t0\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write(pheno_content)
        pheno_file = f.name
    
    try:
        args = {
            "PC": eigenvec_file,
            "pheno": pheno_file,
            "covar": None,
            "pheno_filter": "pheno_1<2",
            "fid_col": "FID",
            "iid_col": "IID",
        }
        
        ids = pd.DataFrame({
            "FID": ["0"] * 6,
            "IID": ["S1", "S2", "S3", "S4", "S5", "S6"],
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)
        
        df = load_tables(ids, args)
        
        assert "pheno_1" in df.columns, "pheno_1 column missing after load_tables"
        present = df["pheno_1"].dropna()
        # Only S1 (1), S5 (-1), S6 (0) satisfy pheno_1<2
        assert set(present.tolist()) == {1, -1, 0}, f"Unexpected pheno_1 values: {present.tolist()}"
        # S2 (5), S3 (2), S4 (3) were filtered out → NaN after left join
        assert pd.isna(df.loc[df["IID"] == "S2", "pheno_1"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "S3", "pheno_1"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "S4", "pheno_1"].iloc[0])
        
    finally:
        os.unlink(eigenvec_file)
        os.unlink(pheno_file)


def test_loading_custom_na_values_from_config():
    """Test that custom na_values from config (list) are recognized as missing."""
    pheno_content = """FID\tIID\tpheno_1\tpheno_2
0\tS1\t-777\t5
0\tS2\t3\t-888
0\tS3\t4\t6
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write(pheno_content)
        pheno_file = f.name
    try:
        args = {
            "PC": None,
            "pheno": pheno_file,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID",
            "na_values": [-777, -888],
        }
        ids = pd.DataFrame({
            "FID": ["0", "0", "0"],
            "IID": ["S1", "S2", "S3"],
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)

        df = load_tables(ids, args)

        # Sentinel values should be converted to NaN
        assert pd.isna(df.loc[df["IID"] == "S1", "pheno_1"].iloc[0])
        assert pd.isna(df.loc[df["IID"] == "S2", "pheno_2"].iloc[0])
        # Valid values should be preserved
        assert df.loc[df["IID"] == "S1", "pheno_2"].iloc[0] == 5
        assert df.loc[df["IID"] == "S2", "pheno_1"].iloc[0] == 3
        assert df.loc[df["IID"] == "S3", "pheno_1"].iloc[0] == 4
    finally:
        os.unlink(pheno_file)


def test_loading_without_na_values_keeps_sentinels():
    """Test that without na_values config, sentinel values are kept as-is."""
    pheno_content = """FID\tIID\tpheno_1
0\tS1\t-777
0\tS2\t3
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write(pheno_content)
        pheno_file = f.name
    try:
        args = {
            "PC": None,
            "pheno": pheno_file,
            "covar": None,
            "fid_col": "FID",
            "iid_col": "IID",
        }
        ids = pd.DataFrame({
            "FID": ["0", "0"],
            "IID": ["S1", "S2"],
        })
        ids["FID"] = ids["FID"].astype(str)
        ids["IID"] = ids["IID"].astype(str)

        df = load_tables(ids, args)

        # Without na_values, -777 is a real value, not NaN
        assert df.loc[df["IID"] == "S1", "pheno_1"].iloc[0] == -777
        assert df.loc[df["IID"] == "S2", "pheno_1"].iloc[0] == 3
    finally:
        os.unlink(pheno_file)


def test_read_flags_na_values_from_config():
    """Test that na_values is exposed from the config (JSON argfile) as a list."""
    import json
    config = {
        "na_values": [-777, -888],
        "PC": None,
        "pheno": None,
        "covar": None,
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(config, f)
        config_file = f.name
    try:
        args = read_flags({"argfile": config_file})
        assert args["na_values"] == [-777, -888]
    finally:
        os.unlink(config_file)


def test_read_flags_na_values_normalization():
    """Test that na_values accepts a single value and coerces numeric strings."""
    from Estimate.data_input.parser import normalize_na_values
    assert normalize_na_values(None) is None
    assert normalize_na_values(-777) == [-777]
    assert normalize_na_values(["-777", "-888"]) == [-777, -888]
    assert normalize_na_values(["NA", "-888"]) == ["NA", -888]


if __name__ == "__main__":
    # Run tests manually if executed directly
    test_loading_eigenvec_hash_fid_uppercase_pcs()
    print("✓ test_loading_eigenvec_hash_fid_uppercase_pcs passed")
    
    test_loading_eigenvec_hash_fid_with_n_a()
    print("✓ test_loading_eigenvec_hash_fid_with_n_a passed")
    
    test_loading_eigenvec_participant_id_custom_iid()
    print("✓ test_loading_eigenvec_participant_id_custom_iid passed")
    
    test_loading_eigenvec_no_header_fallback()
    print("✓ test_loading_eigenvec_no_header_fallback passed")
    
    test_covar_filter_discrete()
    print("✓ test_covar_filter_discrete passed")
    
    test_pheno_filter_numeric()
    print("✓ test_pheno_filter_numeric passed")
    
    test_loading_custom_na_values_from_config()
    print("✓ test_loading_custom_na_values_from_config passed")
    
    test_loading_without_na_values_keeps_sentinels()
    print("✓ test_loading_without_na_values_keeps_sentinels passed")
    
    test_read_flags_na_values_from_config()
    print("✓ test_read_flags_na_values_from_config passed")
    
    test_read_flags_na_values_normalization()
    print("✓ test_read_flags_na_values_normalization passed")
    
    print("\nAll tests passed!")