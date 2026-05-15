# MASH Config Reference

## Overview

MASH uses JSON configuration files to specify all parameters for heritability estimation. This document describes the current config format (as of April 2026) which uses `qcovar`, `covar_discrete`, and `npc` instead of the deprecated `fixed_effects` parameter.

## Basic Structure

```json
{
  "PC": "path/to/pcas.eigenvec",
  "covar": "path/to/covariates.covar",
  "prefix": "path/to/GRM_prefix",
  "pheno": "path/to/phenotypes.pheno",
  "out": "path/to/output_results",
  "npc": 2,
  "mpheno": [1],
  "Method": "GCTA",
  "qcovar": ["age", "bmi"],
  "covar_discrete": ["sex", "site"]
}
```

## Parameter Reference

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `prefix` | string | Prefix of GRM files in GCTA binary format (`prefix.grm.bin`, `prefix.grm.N.bin`, `prefix.grm.id`) |
| `pheno` | string or array | Path(s) to phenotype file(s). Supports a single file path, comma-separated string, or array of multiple files (e.g., `["file1.tsv", "file2.tsv"]`). Multiple files are merged by FID/IID. Supports `.parquet`, `.csv`, `.tsv`/`.tab`, or whitespace-delimited (PLINK format) |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `PC` | string | null | Path to principal components file. Supports header (with FID/IID columns, including `#FID` notation) or no-header (PLINK format). PC column names (`PC1, PC2, ...`) are normalized to lowercase |
| `npc` | integer or array | all PCs | Number of PCs to use. If integer, use top N PCs. If array (e.g., `[3, 5, 10]`), run estimates for each value |
| `covar` | string or array | null | Path(s) to covariate file(s). Supports comma-separated string or array of paths. Multiple files merged on FID/IID |
| `qcovar` | array of strings | null | Column names to treat as quantitative covariates (GCTA only). If null and `covar` provided, auto-detected |
| `covar_discrete` | array of strings | null | Column names to treat as discrete/categorical covariates (GCTA only). If null and `covar` provided, auto-detected |
| `mpheno` | integer or array | 1 | Phenotype column(s) to analyze. 1 = third column (after FID, IID). Accepts list for multiple phenotypes |
| `Method` | string | "AdjHE" | Estimation method: `AdjHE`, `GCTA`, `PredLMM`, `SWD`, `Combat`, or `Covbat` |
| `RV` | string | null | Column name to use as random effect (for methods that support it, e.g., AdjHE with site effects) |
| `out` | string | "MASH_results" | Output file prefix (`.csv` and `.log` will be appended) |
| `k` | integer | all | Rows to process at once when restoring GCTA GRM binary file. Lower values save memory for large datasets |
| `std` | boolean | false | Run SAdj-HE (standardized) instead of UAdj-HE (unstandardized). Note: standardized version may have bugs |
| `loop_covs` | boolean | false | Loop over covariates, including all previous covariates in each iteration |
| `pheno_filter` | string | null | Filter phenotype data (e.g., `"age>30"`, `"sex==1"`). Supports `==`, `!=`, `>`, `<`, `>=`, `<=` |
| `covar_filter` | string | null | Filter covariate data (e.g., `"site==1"`, `"age>=18"`). Same operators as `pheno_filter` |
| `iid_col` | string | `"IID"` | Custom name for the IID column in input files. Set to `"participant_id"` when files use that column name. When `participant_id` is found without a separate FID column, FID is set equal to IID automatically |
| `fid_col` | string | `"FID"` | Custom name for the FID column in input files. Typically kept as `"FID"` |

## Covariate Specification

### New Approach (Recommended)

Use `qcovar` for quantitative covariates and `covar_discrete` for categorical covariates:

```json
{
  "covar": "data/covariates.covar",
  "qcovar": ["pc1", "pc2", "age", "bmi"],
  "covar_discrete": ["sex", "site"]
}
```

### Auto-Detection (GCTA Only)

If `qcovar` and `covar_discrete` are both null, MASH auto-detects covariate types:

- Columns with <35 unique values (and >1 unique value) → discrete/categorical
- Columns with ≥35 unique values → quantitative

```json
{
  "covar": "data/covariates.covar",
  "qcovar": null,
  "covar_discrete": null
}
```

### PC Columns

PC columns (named `pc1`, `pc2`, etc.) are handled separately via `npc` parameter:

```json
{
  "PC": "data/pcas.eigenvec",
  "npc": 2
}
```

PC columns are automatically filtered from `qcovar`/`covar_discrete` to avoid duplication.

## GCTA-Specific Notes

For GCTA method:
- `qcovar` columns are passed via `--qcovar` flag to GCTA
- `covar_discrete` columns are passed via `--covar` flag to GCTA
- PC columns specified by `npc` are included in `qcovar`
- The `RV` (random variable) parameter is NOT passed to GCTA (GCTA doesn't support random effects)

## Example Configs

### GCTA with Auto-Detected Covariates

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_GCTA_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "GCTA",
  "qcovar": null,
  "covar_discrete": null
}
```

### GCTA with Explicit Covariate Specification

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_GCTA_named_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "GCTA",
  "qcovar": ["pc1", "pc2", "age"],
  "covar_discrete": ["sex"]
}
```

### AdjHE with Multiple Phenotype Files and Custom IID

```json
{
  "PC": "tests/test_data/eigen_dummy.eigenvec",
  "covar": ["tests/test_data/age_dummy.tsv", "tests/test_data/sex_dummy.tsv"],
  "prefix": "tests/test_data/EUR3",
  "pheno": ["tests/test_data/nihtb_dummy.tsv", "tests/test_data/smri_dummy.tsv"],
  "out": "tests/results/EUR_AdjHE_results_complex",
  "npc": [2],
  "mpheno": ["phenotype_1", "phenotype_2", "phenotype_3",
    "phenotype_4", "phenotype_5", "phenotype_6"],
  "Method": "AdjHE",
  "qcovar": ["ab_p_demo_age"],
  "covar_discrete": ["ab_g_stc__cohort_sex"],
  "iid_col": "participant_id",
  "fid_col": "FID"
}
```
*This example uses multiple phenotype files merged by FID/IID, `#FID` notation in the PC file header, `participant_id` as the IID column name, and automatic `session_id` handling.*

### AdjHE with Site as Random Effect

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_AdjHE_RV_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "AdjHE",
  "RV": "site",
  "qcovar": ["age", "sex"]
}
```

### AdjHE Basic (No Random Effect)

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_AdjHE_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "AdjHE",
  "qcovar": ["age", "sex"]
}
```

### SWD Method

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_SWD_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "SWD",
  "qcovar": ["age", "sex"]
}
```

### Combat/Covbat Method

```json
{
  "PC": "tests/test_data/EUR2.eigenvec",
  "covar": "tests/test_data/EUR_covariates.cov",
  "prefix": "tests/test_data/EUR2",
  "pheno": "tests/test_data/EUR_simulation2.pheno",
  "out": "tests/results/EUR_COMBAT_results",
  "npc": 2,
  "mpheno": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "Method": "Combat",
  "qcovar": ["age", "sex"],
  "covar_discrete": ["site"]
}
```

## Migration from Old Format

If you have configs using the deprecated `fixed_effects` parameter:

**Old format:**
```json
{
  "fixed_effects": ["pc1", "pc2", "age", "sex"]
}
```

**New format:**
```json
{
  "npc": 2,
  "qcovar": ["age", "sex"]
}
```

**Mapping:**
- PC column names in `fixed_effects` → specify via `npc` (integer) and omit from `qcovar`/`covar_discrete`
- Non-PC columns in `fixed_effects` → add to `qcovar` (if quantitative) or `covar_discrete` (if categorical)

## File Format Support

All input files (phenotype, covariate, PC) support multiple formats:

| Extension | Delimiter | Header | Notes |
|-----------|-----------|--------|-------|
| `.parquet` | N/A | N/A | Apache Parquet format |
| `.csv` | comma (`,`) | recommended | Comma-separated values |
| `.tsv` or `.tab` | tab (`\t`) | recommended | Tab-separated values |
| `.txt`, `.phe`, `.phen`, `.covar` | whitespace | optional | PLINK format (default) |

All files should have the first two columns as FID (Family ID) and IID (Individual ID).

**#FID support:** PC files can use `#FID` as the first column header (e.g., `#FID IID PC1 PC2 ...`). This is automatically detected and normalized to `FID`.

**Custom IID column names:** Use `"iid_col"` to specify a custom IID name (e.g., `"participant_id"`). When `participant_id` is found without a separate FID column, FID is automatically set to match IID.

**session_id handling:** If any covariate or phenotype file contains a `session_id` column, it is automatically dropped during loading to prevent duplicate column conflicts when merging multiple files.

## Output Format

MASH produces a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `h2` | Heritability estimate |
| `var(h2)` | Variance of h2 estimate |
| `G` | Genetic variance component |
| `E` | Environmental variance component |
| `pval` | P-value (if available) |
| `pheno` | Phenotype name |
| `PCs` | Number of PCs used |
| `Covariates` | Covariate names separated by "+" |
| `time` | Computational time (seconds) |
| `mem` | Peak memory usage (MB) |
| `N` | Sample size |

## Testing Your Config

Use the test suite to verify your config works:

```bash
# Activate conda environment
conda activate MASH

# Run all config tests
python -m pytest tests/test_EUR_configs.py -v

# Run specific method tests
python -m pytest tests/test_EUR_configs.py -v -m GCTA
python -m pytest tests/test_EUR_configs.py -v -m AdjHE
```
