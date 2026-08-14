![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Tests](https://img.shields.io/badge/tests-14%20passing-brightgreen?style=for-the-badge)



# adjustedHE

Adj-HE (Adjusted HE) is a computationally efficient method to estimate [Single Nucleotide Polymorphism (SNP)](https://www.cancer.gov/publications/dictionaries/genetics-dictionary/def/single-nucleotide-polymorphism)-heritability in presence of population substructure for biobank-scale data. It is a simplification of the [Haseman- Elston regression (HE)](https://pubmed.ncbi.nlm.nih.gov/4157472/). For details of this statistical method, please refer/cite:

## AdjHE
Lin, Z., Seal, S., & Basu, S. (2020). Estimating SNP heritability in presence of population substructure in large biobank-scale data. bioRxiv. https://doi.org/10.1101/2020.08.05.236901

## AdjHE-RE (Site effects extension)
Coffman C, Feczko E, Larsen B, Tervo-Clemmens B, Conan G, Lundquist JT, Houghton A, Moore LA, Weldon K, McCollum R, Perrone AJ, Fayzullobekova B, Madison TJ, Earl E, Dominguez OM, Fair DA, Basu S. Heritability estimation of subcortical volumes in a multi-ethnic multi-site cohort study. bioRxiv [Preprint]. 2024 Jan 12:2024.01.11.575231. doi: 10.1101/2024.01.11.575231. PMID: 38260520; PMCID: PMC10802572.

# Installation
In desired virtual environment install as
```bash
pip install git+https://github.com/UMN-GDC/MASH.git
```

which is then loadable within Python as `import gdcMASH` or callable from command line as 

```bash
MASH --argfile <path to argfile>
```

If you prefer to have an image so the dependent technologies are handeled for you you can install it with apptainer (singularity) as follows 

```bash
curl -o container.def https://raw.githubusercontent.com/UMN-GDC/MASH/refs/heads/master/container.def
apptainer build MASH.sif container.def
```

which is then callable with


```bash
apptainer run MASH.sif MASH --argfile <path to argfile>
```

## Testing

MASH includes a test suite that validates configuration and estimation methods. All 14 tests currently pass.

### Running Tests
Activate the MASH conda environment first, then run:
```bash
python -m pytest tests/test_EUR_configs.py -v
```
Or run without activating the environment:
```bash
conda run -n Mash python -m pytest tests/test_EUR_configs.py -v
```

### Test Organization
Tests are grouped into classes with pytest marks for easy selection:
| Test Class | Marker | Description |
|------------|--------|-------------|
| `TestAdjHE` | `AdjHE` | Tests AdjHE method |
| `TestGCTA` | `GCTA` | Tests GCTA method |
| `TestSWD` | `SWD` | Tests SWD method |
| `TestCovbat` | `Covbat` | Tests Combat/Covbat method |

Run a specific test class with:
```bash
python -m pytest tests/test_EUR_configs.py -m GCTA -v
```

## Adjusted-HE with closed form

```MASH``` (or ```Estimate.py```) estimates SNP-heritability via closed form formula with single [Genetic Relatedness Matrix (GRM)](https://ibg.colorado.edu/cdrom2020/medland/tuesday1/Tuesday1.pdf) as input. It is suggested to use this version on a server with sufficient memory when sample size is less than 100k. In our paper, analyzing a 45k sample took less than 2 minutes and about 40 GB memory.

Please check the input description with ```MASH --help```.

## Arguments

It is reccomended that users define a `.json` file containing all of the arguments for analysis. This will help both with organization and with reproducibility. This means that the **argfile would be the only argument**. Users can also define all filepaths and variable selections manually using command line flags if desired.

| Input | Description |
|:-------|-------------------|
| --argfile ARGFILE.json | COND REQUIRED. ARGFILE.json, *string*, is the filename to be passed containing all information for PC's, covariates, phenotypes, and GRM. This takes priority over all other arguments. [See the example arfile included under the Example directory.](https://github.com/UMN-GDC/MASH/blob/master/Example/Argfile.json) |
| --prefix PREFIX|  REQUIRED. *string* PREFIX is the prefix of GRM file with GCTA binary GRM format. (`PREFIX.grm.bin`, `PREFIX.grm.N.bin` and `PREFIX.grm.id`)|
| --pheno PHENO |  REQUIRED. PHENO, *string or array*, is the phenotype file. Supports a single file path, a comma-separated string, or an array of multiple files (e.g., `["file1.tsv", "file2.tsv"]`). Multiple files are merged by FID/IID. Supports `.parquet`, `.csv`, `.tsv`/`.tab`, or whitespace-delimited formats. First two columns must be FID and IID, followed by phenotype columns. |
| --mpheno m| OPTIONAL. *list of integers, integer, or "ALL"*, Default=1. If you have multiple phenotypes in the file, you can specify by `--mpheno m`. Otherwise, the first phenotype will be used. Note that 1 refers to the third column of the file since we skip over the FID and IID columns. If passed a list, estimates will be computed for every phenotype specified. Use `"ALL"` (case-insensitive) to run across **all** phenotypes in the phenotype file. |
| --PC PC | OPTIONAL. PC, *string*, is the name of PCs file. Supports header (with FID/IID columns, including `#FID` notation) or no-header (PLINK format). Column names can be `PC1, PC2` or `pc_1, pc_2`. |
| --npc n | OPTIONAL. *integer*, Default = all PCs in the PC file will be used. You can specify top n PCs to be adjusted by `--npc n`.|
| --covar COVAR | OPTIONAL. COVAR, *string or array*, is the name of covariate file(s). Supports multiple files as comma-separated string or array of paths. Multiple files are merged on FID/IID. For GCTA, if no explicit qcovar or covar_discrete specified, columns will be auto-classified as discrete (<35 unique values) or quantitative. |
| --qcovar QCOVARS | OPTIONAL. List of quantitative covariate column names (e.g., `--qcovar age bmi`). Setting `qcovar: null` uses **all** covariate columns; setting `qcovar: []` uses **no** covariates. For GCTA, if not specified, columns are auto-classified as quantitative or discrete based on number of unique values. |
| --covar_discrete DISCRETE_COVARS | OPTIONAL. List of discrete/categorical covariate column names (e.g., `--covar_discrete sex site`). Setting `covar_discrete: null` uses **all** covariate columns; setting `covar_discrete: []` uses **no** covariates. For GCTA, if not specified, columns are auto-classified based on number of unique values. |
| --Method METHOD | OPTIONAL. Specify estimation method: `AdjHE` (default), `GCTA`, `PredLMM`, `SWD`, `Combat`, or `Covbat`. |
| --RV RANDOM_VAR | OPTIONAL. Specify a column name to use as a random effect (for methods that support it, e.g., AdjHE with site effects). |
| --iid_col COL | OPTIONAL. Name of the IID column in input files. Default: `"IID"`. Supports custom names like `"participant_id"`. When `participant_id` is found without a separate FID column, FID is set equal to IID automatically. |
| --k k| OPTIONAL. *integer*. You can specify the number of rows in restoring the GCTA GRM binary file into matrix each time. If not provide, it will process the whole GRM at one time. When you have a relative large sample size, specifying `--k k` can speed up the computation and save the memory. |
| --std | OPTIONAL. Run SAdj-HE by specifying `--std`. Otherwise, UAdj-HE will be computed.  (There are potential bugs with the standardized version, so it is recommended to use unstandardized for now).|
| --loop_covs| OPTIONAL: Default= False. If True, loop over the ORDERED set of covariates including all previous covariates in each iteration. **Note: The order in which the covariates are controlled for is based upon the researcher's best judgements.**|
| --pheno_filter FILTER | OPTIONAL. Filter phenotype data by a condition (e.g., `"age>30"`, `"sex==1"`). Supports: `==`, `!=`, `>`, `<`, `>=`, `<=`. Column names must exist in the phenotype file. |
| --covar-filter FILTER | OPTIONAL. Filter by covariate values (e.g., `"site==1"`, `"age>=18"`). Same operators as pheno_filter. Column names must exist in the covariate file. |
| --na-values NA_VALUES | OPTIONAL. Additional values to recognize as NA/missing in delimited input files (e.g., `--na-values -777 -888`). Accepts a list of values. Sentinel codes like `-777`/`-888` (common in ABCD data) will then be treated as missing instead of real numbers. |

## Description of Input File Formats

MASH supports multiple file formats with automatic format detection based on file extension:

| Extension | Delimiter | Header | Notes |
|-----------|-----------|--------|-------|
| `.csv` | comma (`,`) | recommended | Comma-separated values |
| `.tsv` or `.tab` | tab (`\t`) | recommended | Tab-separated values |
| `.parquet` | N/A | N/A | Apache Parquet format |
| `.txt`, `.phe`, `.phen`, `.covar` | whitespace | not required | PLINK format (default) |

All files should have the first two columns as FID (Family ID) and IID (Individual ID). Column headers are recommended and will be preserved for use in `qcovar`, `covar_discrete`, and `npc`. PC column names are normalized to lowercase (`PC1` → `pc1`) for consistency.

**#FID support:** PC files can use `#FID` (hashed FID) as the first column header instead of `FID`. This is automatically detected and normalized. Files with `#FID` and `IID` as first two columns followed by capital `PC1, PC2, ...` columns are supported.

**Custom IID column names:** Use `"iid_col"` to specify a custom IID column name (e.g., `"participant_id"`). When `participant_id` is found without a separate FID column, FID is automatically set equal to IID. This is common in neuroimaging datasets.

**session_id handling:** If any covariate or phenotype file contains a `session_id` column, it is automatically dropped during loading to prevent duplicate column conflicts when merging multiple files.

### Data Filtering
You can filter your phenotype or covariate data at load time using `--pheno_filter` and `--covar_filter`:

**Supported operators:**
- `==` (equal)
- `!=` (not equal)
- `>` (greater than)
- `<` (less than)
- `>=` (greater or equal)
- `<=` (less or equal)

**Examples:**
```json
{
  "pheno_filter": "age>30",
  "covar_filter": "sex==1"
}
```
```bash
# Command line usage
MASH --pheno data/phenotypes.tsv --covar data/covariates.tsv --pheno-filter "age>=18" --covar-filter "site!=2"
```

Filters are applied before merging, so `--pheno-filter` affects phenotype rows and `--covar-filter` affects covariate rows. Only subjects passing both filters (and present in the GRM) will be included in analysis.

### Missing Values (`na_values`)
By default, MASH recognizes the standard missing markers (empty values, `NA`, `NaN`, `n/a`, etc.). Some datasets (e.g., ABCD) encode missing data with sentinel values such as `-777` or `-888`. Use `na_values` to declare additional values that should be treated as missing:

```json
{
  "na_values": [-777, -888]
}
```
```bash
# Command line usage
MASH --argfile config.json --na-values -777 -888
```

Custom NA values are *added* to the built-in list, so you don't need to re-specify standard markers. `na_values` applies to all delimited text files (phenotype, covariate, and PC files).

### Covariate Selection (`qcovar` / `covar_discrete`)
Control which covariates are used in the model:

| Setting | Behavior |
|---------|----------|
| `null` (default) | Use **all** covariate columns |
| `[]` | Use **no** covariates |
| `["age", "sex"]` | Use exactly those columns |

```json
{
  "qcovar": null,
  "covar_discrete": null
}
```
```json
{
  "qcovar": [],
  "covar_discrete": []
}
```

For the `GCTA` method, `null` triggers automatic classification of covariate columns into quantitative vs. discrete based on the number of unique values.

### Running All Phenotypes (`mpheno: "ALL"`)
Set `mpheno` to `"ALL"` (case-insensitive) to estimate heritability across **every** phenotype in the phenotype file:

```json
{
  "mpheno": "ALL"
}
```
```bash
# Command line usage
MASH --argfile config.json --mpheno ALL
```

### Example Input Files

**Phenotype file (.csv with header):**
```csv
FID,IID,pheno_1,pheno_2
0, EUR20_EUR20, 1.2, 3.4
0, EUR60_EUR60, 2.1, 1.8
```

**Covariate file (.tsv with header):**
```tsv
FID	IID	site	age	sex
0	EUR20_EUR20	1	45	M
0	EUR60_EUR20	2	32	F
```

**PC file (PLINK format without header):**
```
0	EUR20_EUR20	-0.017	-0.006
0	EUR60_EUR60	-0.076	0.060
```

### Example Configuration Files

**Basic GCTA with auto-detected covariates:**
```json
{
  "PC": "data/pcas.eigenvec",
  "covar": "data/covariates.covar",
  "prefix": "data/my_GRM",
  "pheno": "data/phenotypes.pheno",
  "out": "results/my_results",
  "npc": [3, 5, 10],
  "mpheno": [1, 2, 3],
  "Method": "GCTA",
  "qcovar": null,
  "covar_discrete": null
}
```
*When qcovar/covar_discrete are null, MASH auto-detects: columns with <35 unique values → discrete, others → quantitative.*

**GCTA with explicit covariate specification:**
```json
{
  "PC": "data/pcas.eigenvec",
  "covar": "data/covariates.covar",
  "prefix": "data/my_GRM",
  "pheno": "data/phenotypes.pheno",
  "out": "results/my_results",
  "npc": [5],
  "mpheno": [1],
  "Method": "GCTA",
  "qcovar": ["PC1", "PC2", "age", "bmi"],
  "covar_discrete": ["sex", "site"]
}
```

**AdjHE with site as random effect:**
```json
{
  "PC": "data/pcas.eigenvec",
  "covar": "data/covariates.csv",
  "prefix": "data/my_GRM",
  "pheno": "data/phenotypes.pheno",
  "out": "results/my_results",
  "npc": [5],
  "mpheno": [1],
  "Method": "AdjHE",
  "qcovar": ["age", "sex"],
  "random_groups": "site"
}
```

**AdjHE with filtering (e.g., include only adults from specific site):**
```json
{
  "PC": "data/pcas.eigenvec",
  "covar": "data/covariates.csv",
  "prefix": "data/my_GRM",
  "pheno": "data/phenotypes.pheno",
  "out": "results/my_results",
  "npc": [5],
  "mpheno": [1],
  "Method": "AdjHE",
  "pheno_filter": "age>=18",
  "covar_filter": "site==1"
}
```

**AdjHE with multiple phenotype files and custom IID column:**
```json
{
  "PC": "data/eigen.eigenvec",
  "covar": ["data/age.tsv", "data/sex.tsv"],
  "prefix": "data/my_GRM",
  "pheno": ["data/nihtb_phenotypes.tsv", "data/mri_phenotypes.tsv"],
  "out": "results/complex_results",
  "npc": [2],
  "mpheno": ["phenotype_1", "phenotype_2", "phenotype_3",
    "phenotype_4", "phenotype_5", "phenotype_6"],
  "Method": "AdjHE",
  "qcovar": ["age"],
  "covar_discrete": ["sex"],
  "iid_col": "participant_id",
  "fid_col": "FID"
}
```
*Multiple phenotype files are merged by FID/IID. The `iid_col` parameter tells MASH to look for `participant_id` instead of `IID`. Files can use `#FID` as the first column header. The `session_id` column is automatically handled.*

### Example of computing GRM and eigenvectors from .bed files
In your own data analysis, you may need to computed the grm and eigenvectors yourself so here are the steps to do that so that you will just be able to run the above examples afterward. Great resrources are [here](https://yanglab.westlake.edu.cn/software/gcta/#PCA) and [here](https://ibg.colorado.edu/cdrom2021/Day04-yengo/Day4_practical_Boulder2021_v4.pdf).
```

# Starting with .bed files
# 1. Calculate the GRM (note you can filter for alleles with unacceptably low MAF's at this step if not already done in Plink)
gcta64 --bfile Example/geno --maf 0.05 --make-grm-bin --out Example/grm

# 2. Calculate eigenvectors
gcta64  --grm Example/grm --pca 20  --out Example/pcas
```



## Output
A .csv with heritability estimate (h2), standard error (SE), phenotype used (Pheno), number of prinicpal components controlled for (PCs), list of covariates controled for separated by a "+" (Covariates), computational time (Time for analysis(s)), and peak memory (Memory Usage (Mb)) are also provided. See the [results](https://github.com/UMN-GDC/MASH/blob/master/Example/results.csv) included in the Example folder.

# Example
Included in this repo is an example dataset for users to practice with and check results to the included results file.

## Data description
Example data is included from [this paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010151) which simulated 10,000 SNP's for 5,000 individuals with the following files

| File Name | Brief description | File contents|
|-----------|-------------------|--------------|
| pheno.phen| phenotype file | First two columns contain FID, IID, one phenotype column (though there can be more) |
| covar.txt | covariate file | First two columns contain FID, IID, one covariate column (though there can be more) |
| geno (.bed, .bim, .fam) | suite of files containing genotypes in PLINK format | For more info see [PLINK](https://www.cog-genomics.org/plink/) |
| grm (.grm.id, .grm.bin, .grm.N.bin) | suite of files GRM's.| For more info see [Making a GRM](https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM) |
| pcas (.eigenval, .eigenvec, .log) | suite of files containing PC's. | For more info see [GCTA PCA](https://yanglab.westlake.edu.cn/software/gcta/#PCA) |
|AdjHE_results.csv | results from example listed below | |
| Argfile.json | File containing all arguments to run AdjHE exampel | |
| Batch_Arg_file.txt | File containing information how to create sets of batch files| |


It is to be carefully noted that the order of the individuals in all the files (phenotype, covariate, GRM) have to be the same.
The example dataset includes simulated phenotypes from the linked paper with a true heritability of 0.8. Also note that the phenotypes in the .fam file are not used.


## Running with included example data 
First, change into the "Basu_herit" directory then run the following commands that specify the arguments
```
module load python
prefix='EXAMPLE/grm' #We partitioned the GRM into 200 parts.
id='/EXAMPLE/ukbiobank.grm.id'
pheno='Example/phen.pheno'
PC='Example/pcas.eigenvec'
covar='Example/covar.csv'
out='Example/results.csv'

python Estimate.py --prefix ${prefix} --PC ${PC} --npc 10  --covar ${covar} --pheno ${pheno} --mpheno 1 --out ${out}
```
This should result in estimates for heritability stored in a .csv with the estimated heritability. For this dataset, the simulated heritability was 80%. Notice that the estimate is sensitive to the number of Prinicipal components included in the model since the data was simulated to have population stratification. The covariates don't have much of an influence on the estimates since they were not included in the simulation of this dataset. Compare your results with the [results included in the Example folder](https://github.com/coffm049/Basu_herit/blob/master/Example/results.csv.

## Running Example with argfile
```
python Estimate.py --argfile Example/Arg_file.json
```
**Note that this is running the same example method as the previous example, only in this case, all of the arguments are contained within the Arg_file.txt file. This helps with reproducibility and creating batch scripts.**

# Creating Batch scripts (Coming soon)

### SLURM script example 
The SLURM script contains two portions: the first request resources from the supercomputer, the second contains the code you actually want to run. In order to make sure your code runs properly, make sure you are changing to the proper working direcotry with the "cd" command. Then simply call the function you want to run. [Here](https://github.com/coffm049/Basu_herit/blob/master/Example/Estimate.SLURM) is an example SLURM script where we are simply running the previous example utilizing the "argfile". For information on the anatomy of SLURM scripts please see information from [MSI](https://www.msi.umn.edu/content/job-submission-and-scheduling-slurm).


# FAQ's
## Sensitivity to number of PC's
See simulation 2 in Section 3 "Results" in [the paper](https://doi.org/10.1101/2020.08.05.236901). the efficiency of this method of estimation allows for users to do some amount of model selection nicluding with the number of PC's which could drastically impact heritability estimates.

<!--![Alt](https://repobeats.axiom.co/api/embed/11759d6c6f5bb629dd90af840be633628d5d6add.svg "Repobeats analytics image")
-->

