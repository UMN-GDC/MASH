# adjustedHE

Adj-HE (Adjusted HE) is a computational efficient method to estimate [Single Nucleotide Polymorphism (SNP)](https://www.cancer.gov/publications/dictionaries/genetics-dictionary/def/single-nucleotide-polymorphism)-heritability in presence of population substructure for biobank-scale data. It is a simplification of the [Haseman- Elston regression (HE)](https://pubmed.ncbi.nlm.nih.gov/4157472/). For details of this statistical method, please refer/cite:
 
Lin, Z., Seal, S., & Basu, S. (2020). Estimating SNP heritability in presence of population substructure in large biobank-scale data. bioRxiv. https://doi.org/10.1101/2020.08.05.236901

## Adjusted-HE with closed form formula

```AdjHE.py```  estimates SNP-heritability via closed form formula with single [Genetic Relatedness Matrix (GRM)](https://ibg.colorado.edu/cdrom2020/medland/tuesday1/Tuesday1.pdf) as input. It is suggested to use this version on a server with sufficient memory when sample size is less than 100k. In our paper, analyzing a 45k sample took less than 2 minutes and about 40 GB memory.

Please check the input description with ```./AdjHE.py --help```.

## Arguments
It is reccomended that users define a .json file containing all of the arguments for analysis. This will help both with organization and with reproducibility. This means that the **argfile would be the only argument**. Users can also define all filepaths and variable selections manually using command line flags if desired. 

|&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Input &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description
:-------|-------------------
| --argfile ARGFILE.json | COND REQUIRED. ARGFILE.json, *string*, is the filename to be passed containing all information for PC's, covariates, phenotypes, and GRM. This takes priority over all other arguments. [See the example arfile included under the Example directory.](https://github.com/coffm049/Basu_herit/blob/master/Example/Argfile.json). |
| --prefix PREFIX|  REQUIRED. *string* PREFIX is the prefix of GRM file with GCTA binary GRM format. (```PREFIX.grm.bin```, ```PREFIX.grm.N.bin``` and ```PREFIX.grm.id```)|
| --pheno PHENO.phen |  REQUIRED. PHENO.phen, *string*, is the name of phenotype file, following GCTA phenotype file format (space delimited text file) but with column names). The first two columns are FID and IID and phenotypes start from the third column. |
| --mpheno m| OPTIONAL. *integer*, Default=1. If you have multiple phenotypes in the file, you can specify by ```--mpheno m```. Otherwise, the first phenotype will be used.|
| --PC PC | OPTIONAL. PC, *string*, is the name of PCs file, following GCTA (space delimited, no column names) ```--pca``` file (same as plink ```--pca```). The third column is the first PC, the forth column is the second PC...|
| --npc n | OPTIONAL. *integer*, Default = all PCs in the PC file will be used. You can specify top n PCs to be adjusted by ```--npc n```.|
| --covar COVAR | OPTIONAL. COVAR, *string*, is the name of covariate file, following GCTA ```--qcovar``` file format or .csv file format. It may contain sex, age, etc. *Note that this file does not include principal components, which need to be include seperately by ```--PC PC```*.|
|--covars COVARS| OPTIONAL. COVARS is the *list of integers* specifying which covariates to control for from the covariate file. column numbering does not include the FID and IID columns. **Note that this is an ordered list if used in conjunction with the ```--loop_covs``` flag.**|
| --k k| OPTIONAL. *integer*. You can specify the number of rows in restoring the GCTA GRM binary file into matrix each time. If not provide, it will process the whole GRM at one time. When you have a relative large sample size, specifying ```--k k``` can speed up the computation and save the memory. |
| --std | OPTIONAL. Run SAdj-HE by specifying ```--std```. Otherwise, UAdj-HE will be computed.  (There are potential bugs with the standardized version, so it is reccommended to use unstandardized for now).|
| --loop_covs| OPTIONAL: Default= False. If True, loop over the ORDERED set of user defined covariates including all previous covariates in each iteration. **Note: The order in which the covrariates are controlled for is based upon the researchers best judjements. In other words, include the most likely **|

## Descrtiption of Inputs
Here are illustrative examples of what files might look like. All three file types containing data for analysis have the first two columns that are the Family ID (FID) and the Individuals ID (IID). They are then followed by values specific to each file type (phenotypes for the phenotype file, covariates for the covaraiates file, and PC loadings for the PC file. **Note: Both the phenotype and covariates files should have column headers, whereas the PC file should not.** See the example [pheno](https://github.com/coffm049/Basu_herit/blob/master/Example/pheno.phen), [covariate](https://github.com/coffm049/Basu_herit/blob/master/Example/covar.txt), and [PC](https://github.com/coffm049/Basu_herit/blob/master/Example/pcas.eigenvec) file in the examples folder.

| Column | Column Contents |
|------------|----------|
| 1 | FID,*string*, Unique family identifier for each family in the study. |
| 2 | IID, *string*, Unique individual identifier for each individual in the study.|
| 3-infinite | *numeric*, measurements for PC loading, phenotype, or covariate measures depending on file type |

**Note: All files need the first two columns to be FID and IID, respectively. Also any missing values will remove the observation from the given analysis. **


## Output
A .csv with heritability estimate (h2), standard error (SE), phenotype used (Pheno), number of prinicpal components controlled for (PCs), list of covariates controled for separated by a "+" (Covariates), computational time (Time for analysis(s)), and peak memory (Memory Usage (Mb)) are also provided. See the [results](https://github.com/coffm049/Basu_herit/blob/master/Example/results.csv) included in the Example folder.

# Example
Included in this repo is an example dataset for users to practice with and check results to the included results file.

## Data description
Example data is included from [this paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010151) with the following files

    a phenotype file: pheno.phen
    a covariate file: covar.txt
    PLINK Binary files: geno (.bed, .bim, .fam)
    GCTA GRM files: grm (.grm.id, .grm.bin, .grm.N.bin)
    eigenvector files: pcas (.eigenval, .eigenvec, .log)

There are 5000 individuals and 10,000 SNPs. The first two columns of the phenotype and the covariate files have the FID and individual ID IID of each individual. The phenotype file has a single phenotype and the covariate file has a single covariate. With the binary files, the GRM files have been computed using GCTA. It is to be craefully noted that the order of the individuals in all the files (phenotype, covariate, GRM) have to be the same.
The example dataset includes simulated phenotypes from the linked paper with a true heritability of 0.8.


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

python AdjHE.py --prefix ${prefix} --PC ${PC} --npc 10  --covar ${covar} --pheno ${pheno} --mpheno 1 --out ${out}
```
This should result in estimates for heritability stored in a .csv with the estimated heritability. For this dataset, the simulated heritability was 70%. Notice that the estimate is sensitive to the number of Prinicipal components included in the model since the data was simulated to have population stratification. The covariates don't have much of an influence on the estimates since they were not included in the simulation of this dataset. Compare your results with the [results included in the Example folder](https://github.com/coffm049/Basu_herit/blob/master/Example/results.csv.

## Running Example with argfile
```
python AdjHE.py --argfile Example/Arg_file.txt
```
**Note that this is running the same example method as the previous example, only in this case, all of the arguments are contained within the Arg_file.txt file. This helps with reproducibility and creating batch scripts.**

# Creating Batch scripts (Coming soon)

### SLURM script example 
The SLURM script contains two portions: the first request resources from the supercomputer, the second contains the code you actually want to run. In order to make sure your code runs properly, make sure you are changing to the proper working direcotry with the "cd" command. Then simply call the function you want to run. In this case, we are simply running the previous example utilizing the "argfile". For information on the anatomy of SLURM scripts please see information from [MSI](https://www.msi.umn.edu/content/job-submission-and-scheduling-slurm).
```
#!/bin/bash -l        
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=sample_email@umn.edu 

cd ~/PATH/TO/DIR/
module load python 
python AdjHE.py --argfile Example/Arg_file.txt
```

### Example of computing GRM and eigenvectors from .bed files
In your own data analysis, you may need to computed the grm and eigenvectors yourself so here are the steps to do that so that you will just be able to run the above examples afterward. Great resrources are [here](https://yanglab.westlake.edu.cn/software/gcta/#PCA) and [here](https://ibg.colorado.edu/cdrom2021/Day04-yengo/Day4_practical_Boulder2021_v4.pdf).
```

# Starting with .bed files
# 1. Calculate the GRM (note you can filter for alleles with unacceptably low MAF's at this step if not already done in Plink)
gcta64 --bfile Example/geno --maf 0.05 --make-grm-bin --out Example/grm

# 2. Calculate eigenvectors
gcta64  --grm Example/grm --pca 20  --out Example/pcas
```
Then follow the examples above directly.



# UNDER CONSTRUCTION


## Adjusted-HE with regression version (Still being tested)

For large sample size (e.g. biobank size), it is suggested to use ```AdjHE_reg_s1.py``` and ```AdjHE_reg_s2.py``` to perform linear regression to get the heritability estimation. In our estimation of the UKB data, we partitioned the whole data set (n = 305,639) into 200 parts and allocated 15GB memory to each job. And the maximum running time for 200 jobs was less than 10 minutes (not include construcing GRM and computing PCs).



#### ```AdjHE_reg_s1.py```

*You should FIRST call ```AdjHE_reg_s1.py``` before ```AdjHE_reg_s2.py```.*
Please check the input description with ```./AdjHE_reg_s1.py --help```.

Arguments

|&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Input &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description
:-------|-------------------
| --prefix PREFIX|  REQUIRED. The prefix of partitioned GRM file with GCTA binary GRM format. (e.g. ```--prefix ukbiobank.part_200_```) Please refer to ```--make-grm-part``` in https://cnsgenomics.com/software/gcta/#MakingaGRM.|
| --pheno PHENO|  REQUIRED. The name of phenotype file, following GCTA phenotype file format. The first two columns are FID and IID and phenotypes start from the third column. |
| --Npart N| REQUIRED. The number of parts you specify in ```--make-grm-part N i``` when constructing GRM by parts.|
| --id ID| REQUIRED. The file including the ids of *ALL* samples. It follows GCTA ```PREFIX.grm.id``` format. You can create this file by ```cat ukbiobank.part_200_*.grm.id > ukbiobank.grm.id```.|
| --job i| REQUIRED. The current run for computing the i-th part of GRM. (i.e., It would take ukbiobank.part_200_i.* as input) You can easily parallel jobs using PBS Job Array by specifying ```--job ${PBS_ARRAYID}```. (See below)|
| --out DIRECTORY| REQUIRED. This is a existing directory to store all intermediate results and a Log file. After processing all parts of the GRM, it should generate ```Npart``` intermediate files (e.g., ```pheno1_1```,...,```pheno1_200```). The Log file will provide the computation time for each run. |
| --mpheno m| OPTIONAL. If you have multiple phenotypes in the file, you can specify by ```--mpheno m```. Otherwise, the first phenotype will be used.|
| --PC PC| OPTIONAL. The name of PCs file, following GCTA ```--pca``` file (same as plink ```--pca```). The third column is the first PC, the forth column is the second PC...|
| --npc n| OPTIONAL. You can specify top n PCs to be adjusted by ```--npc n```. Otherwise, all PCs in the PC file will be used.|
| --covar COVAR| OPTIONAL. The name of covariate file, following GCTA ```--qcovar``` file format. It may contain sex, age, etc. *Note that this file does not include principal components, which need to be include seperately by ```--PC PC```*.|
| --std | OPTIONAL. Run SAdj-HE by specifying ```--std```. Otherwise, UAdj-HE will be computed.|



#### ```AdjHE_reg_s2.py```

*You should only call ```AdjHE_reg_s2.py``` once after finishing ```AdjHE_reg_s1.py``` with N parts GRM successfully.* Please check the input description with ```./AdjHE_reg_s2.py --help```.

Arguments

|&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Input &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description
:-------|-------------------
| --pheno PHENO|  REQUIRED. This should be consistent with ```--pheno``` in ```AdjHE_reg_s1.py```. |
| --Npart N| REQUIRED. This should be consistent with ```--Npart``` in ```AdjHE_reg_s1.py```.|
| --id ID| REQUIRED. This should be consistent with ```--id``` in ```AdjHE_reg_s1.py```.|
| --out DIRECTORY| REQUIRED. This should be consistent with ```--out``` in ```AdjHE_reg_s1.py```. Make sure you have m intermediate result files in this directoory. The point estimate and standard error of heritability will be output to the Log file generated in ```AdjHE_reg_s1.py```.|
| --mpheno m| OPTIONAL. This should be consistent with ```--mpheno``` in ```AdjHE_reg_s1.py```.|
| --PC PC| OPTIONAL. This should be consistent with ```--PC``` in ```AdjHE_reg_s1.py```|
| --npc n| OPTIONAL. This should be consistent with ```--npc``` in ```AdjHE_reg_s1.py```.|
| --covar COVAR| OPTIONAL. This should be consistent with ```--covar``` in ```AdjHE_reg_s1.py```.|
| --std | OPTIONAL. This should be consistent with ```--std``` in ```AdjHE_reg_s1.py```|

Portable Batch System (or simply PBS) is the name of computer software that performs job scheduling. Its primary task is to allocate computational tasks, i.e., batch jobs, among the available computing resources. It is often used in conjunction with UNIX cluster environments.

Here we give an example of running regression adjusted-HE for large sample size in parallel using PBS scripts.

The first script below ```reg_s1.pbs``` is to run ```AdjHE_reg_s1.py```.

```
#!/bin/bash -l
#PBS -l walltime=01:30:00,nodes=1:ppn=10,pmem=2500MB
#PBS -m a
#PBS -M user@email.com
#PBS -j oe
#PBS -o ./pbs
#PBS -e ./pbs
#PBS -t 1-200

module load python
prefix='/PATH/TO/PARTGRM/ukbiobank.part_200_' #We partitioned the GRM into 200 parts.
id='/PATH/TO/PARTGRM/ukbiobank.grm.id'
pheno='/PATH/TO/PHENO/ukbiobank.pheno'
PC='/PATH/TO/PC/ukbiobank.eigenvec'
covar='/PATH/TO/COVAR/ukbiobank_sex_age.txt'
out='/PATH/TO/RESULT/DIRECTORY'

python AdjHE_reg_s1.py --prefix ${prefix} --PC ${PC} --npc 10  --covar ${covar} --pheno ${pheno} --mpheno 1 --job ${PBS_ARRAYID} --Npart 200 --id ${id} --out ${out}
```

The second script below ```reg_s2.pbs``` is to run ```AdjHE_reg_s2.py```.
```
#!/bin/bash -l
#PBS -l walltime=1:00:00,nodes=1:ppn=10,pmem=2580MB
#PBS -m a
#PBS -M user@email.com
#PBS -j oe
#PBS -o ./pbs
#PBS -e ./pbs

module load python
id='/PATH/TO/PARTGRM/ukbiobank.grm.id'
pheno='/PATH/TO/PHENO/ukbiobank.pheno'
PC='/PATH/TO/PC/ukbiobank.eigenvec'
covar='/PATH/TO/COVAR/ukbiobank_sex_age.txt'
out='/PATH/TO/RESULT/DIRECTORY'

python AdjHE_reg_s2.py --PC ${PC} --npc 10  --covar ${covar} --pheno ${pheno} --mpheno 1 --job ${PBS_ARRAYID} --Npart 200 --id ${id} --out ${out}
```
The result will be saved in /PATH/TO/RESULT/DIRECTORY/pheno1.log .


# Building docker  
```
sudo docker build -f Dockerfile -t adjhe:1.1 .
sudo docker save adjhe:1.1 > AdjHE_1.1.tar
```


Start an interactive job on mesabi that you will use to build the container. The size of the /tmp directory allocated with the –tmp flag will be dependent on the size of the container you are trying to build. A rule of thumb that has seemed to work okay is to request at least twice the size of the container you are trying to build. 

```
srun -N 1 --ntasks-per-node=1  --tmp=100g --mem-per-cpu=30g -t 3:00:00 -p interactive --pty bash
```

Load singularity, go to the directory where you have the .tar file for the image, and build it as a singularity image:

```
module load singularity
cd ~/tools/Basu_herit
singularity build AdjHE.sif docker-archive://AdjHE_1.1.tar
```






