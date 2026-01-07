![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)


# adjustedHE

Adj-HE (Adjusted HE) is a computationally efficient method to estimate [Single Nucleotide Polymorphism (SNP)](https://www.cancer.gov/publications/dictionaries/genetics-dictionary/def/single-nucleotide-polymorphism)-heritability in presence of population substructure for biobank-scale data. It is a simplification of the [Haseman- Elston regression (HE)](https://pubmed.ncbi.nlm.nih.gov/4157472/). For details of this statistical method, please refer/cite:
 
Lin, Z., Seal, S., & Basu, S. (2020). Estimating SNP heritability in presence of population substructure in large biobank-scale data. bioRxiv. https://doi.org/10.1101/2020.08.05.236901

## Adjusted-HE with closed form

```Estimate.py```  estimates SNP-heritability via closed form formula with single [Genetic Relatedness Matrix (GRM)](https://ibg.colorado.edu/cdrom2020/medland/tuesday1/Tuesday1.pdf) as input. It is suggested to use this version on a server with sufficient memory when sample size is less than 100k. In our paper, analyzing a 45k sample took less than 2 minutes and about 40 GB memory.

Please check the input description with ```./Estimate.py --help```.

## Arguments
It is reccomended that users define a .json file containing all of the arguments for analysis. This will help both with organization and with reproducibility. This means that the **argfile would be the only argument**. Users can also define all filepaths and variable selections manually using command line flags if desired. 

|&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Input &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description
:-------|-------------------
| --argfile ARGFILE.json | COND REQUIRED. ARGFILE.json, *string*, is the filename to be passed containing all information for PC's, covariates, phenotypes, and GRM. This takes priority over all other arguments. [See the example arfile included under the Example directory.](https://github.com/coffm049/Basu_herit/blob/master/Example/Argfile.json). |
| --prefix PREFIX|  REQUIRED. *string* PREFIX is the prefix of GRM file with GCTA binary GRM format. (```PREFIX.grm.bin```, ```PREFIX.grm.N.bin``` and ```PREFIX.grm.id```)|
| --pheno PHENO.phen |  REQUIRED. PHENO.phen, *string*, is the name of phenotype file, following GCTA phenotype file format (space delimited text file) but with column names). The first two columns are FID and IID and phenotypes start from the third column. |
| --mpheno m| OPTIONAL. *list of integers or integer*, Default=1. If you have multiple phenotypes in the file, you can specify by ```--mpheno m```. Otherwise, the first phenotype will be used. Note that 1 refers to the third column of the file since we skip over the FID and IID columns. If passed a list, estimates will be computed for every phenotype specified. |
| --PC PC | OPTIONAL. PC, *string*, is the name of PCs file, following GCTA (space delimited, no column names) ```--pca``` file (same as plink ```--pca```). The third column is the first PC, the forth column is the second PC...|
| --npc n | OPTIONAL. *integer*, Default = all PCs in the PC file will be used. You can specify top n PCs to be adjusted by ```--npc n```.|
| --covar COVAR | OPTIONAL. COVAR, *string*, is the name of covariate file, following GCTA ```--qcovar``` file format or .csv file format. It may contain sex, age, etc. *Note that this file does not include principal components, which need to be include seperately by ```--PC PC```*.|
|--covars COVARS| OPTIONAL. COVARS is the *list of integers* specifying which covariates to control for from the covariate file. column numbering does not include the FID and IID columns (therefore 1 refers to the third column of the file). **Note that this is an ordered list if used in conjunction with the ```--loop_covs``` flag.**|
| --k k| OPTIONAL. *integer*. You can specify the number of rows in restoring the GCTA GRM binary file into matrix each time. If not provide, it will process the whole GRM at one time. When you have a relative large sample size, specifying ```--k k``` can speed up the computation and save the memory. |
| --std | OPTIONAL. Run SAdj-HE by specifying ```--std```. Otherwise, UAdj-HE will be computed.  (There are potential bugs with the standardized version, so it is reccommended to use unstandardized for now).|
| --loop_covs| OPTIONAL: Default= False. If True, loop over the ORDERED set of user defined covariates including all previous covariates in each iteration. **Note: The order in which the covrariates are controlled for is based upon the researchers best judjements. In other words, include the most likely **|

## Descrtiption of Inputs
Here are illustrative examples of what files might look like. The phenotyp, covariates, and principal componet files have the first two columns that are the Family ID (FID) and the Individuals ID (IID). They are then followed by values specific to each file type (phenotypes for the phenotype file, covariates for the covaraiates file, and PC loadings for the PC file (An exmaple explaining how to compute PC's is given in the section "Example of computing GRM and eigenvectors from .bed files" below). **Note: Both the phenotype and covariates files should have column headers, whereas the PC file should not.** See the example [pheno](https://github.com/coffm049/Basu_herit/blob/master/Example/pheno.phen), [covariate](https://github.com/coffm049/Basu_herit/blob/master/Example/covar.txt), and [PC](https://github.com/coffm049/Basu_herit/blob/master/Example/pcas.eigenvec) file in the examples folder.

| Column | Column Contents |
|------------|----------|
| 1 | FID,*string*, Unique family identifier for each family in the study. |
| 2 | IID, *string*, Unique individual identifier for each individual in the study.|
| 3-infinite | *numeric*, measurements for PC loading, phenotype, or covariate measures depending on file type |

**Note: All files need the first two columns to be FID and IID, respectively. Also any missing values will remove the observation from the given analysis. **

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
A .csv with heritability estimate (h2), standard error (SE), phenotype used (Pheno), number of prinicpal components controlled for (PCs), list of covariates controled for separated by a "+" (Covariates), computational time (Time for analysis(s)), and peak memory (Memory Usage (Mb)) are also provided. See the [results](https://github.com/coffm049/Basu_herit/blob/master/Example/results.csv) included in the Example folder.

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

????? ALSO LIST THE GENERAL RESOURCES NEEDED TO STARTUP THIS ANALYSIS FOR OTHER TyPES OF CLUSTERS????

# UNDER CONSTRUCTION


## Adjusted-HE with regression version (Still being tested)

For large sample size (e.g. biobank size), it is suggested to use ```AdjHE_reg_s1.py``` and ```AdjHE_reg_s2.py``` to perform linear regression to get the heritability estimation. In our estimation of the UKB data, we partitioned the whole data set (n = 305,639) into 200 parts and allocated 15GB memory to each job. And the maximum running time for 200 jobs was less than 10 minutes (not include construcing GRM and computing PCs).




# Building Image  
```
sudo docker build -f Dockerfile -t adjhe:1.1 .
sudo docker save adjhe:1.1 > AdjHE_1.1.tar

module load singularity
cd ~/tools/Basu_herit
singularity build Estimate.sif docker-archive://AdjHE_1.1.tar

```


Start an interactive job on mesabi that you will use to build the container. The size of the /tmp directory allocated with the –tmp flag will be dependent on the size of the container you are trying to build. A rule of thumb that has seemed to work okay is to request at least twice the size of the container you are trying to build. 

```
srun -N 1 --ntasks-per-node=1  --tmp=100g --mem-per-cpu=30g -t 3:00:00 -p interactive --pty bash
```



# FAQ's
## Sensitivity to number of PC's
See simulation 2 in Section 3 "Results" in [the paper](https://doi.org/10.1101/2020.08.05.236901). the efficiency of this method of estimation allows for users to do some amount of model selection nicluding with the number of PC's which could drastically impact heritability estimates.

![Alt](https://repobeats.axiom.co/api/embed/11759d6c6f5bb629dd90af840be633628d5d6add.svg "Repobeats analytics image")

