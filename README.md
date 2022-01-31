# adjustedHE
## Description

Adj-HE (Adjusted HE) is a computational efficient method to estimate SNP-heritability in presence of population substructure for biobank-scale data. For details of this statistical method, please refer/cite:
 
Lin, Z., Seal, S., & Basu, S. (2020). Estimating SNP heritability in presence of population substructure in large biobank-scale data. bioRxiv. https://doi.org/10.1101/2020.08.05.236901

## Adjusted-HE with closed form formula version

```AdjHE.py```  estimates SNP-heritability via closed form formula with single GRM as input. It is suggested to use this version on a server with sufficient memory when sample size is less than 100k. In our paper, analyzing a 45k sample only used less than 2 minutes and about 40 GB memory.

Please check the input description with ```./AdjHE.py --help```.

Arguments

|&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Input &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Description
:-------|-------------------
| --prefix PREFIX|  REQUIRED. The prefix of GRM file with GCTA binary GRM format. (```PREFIX.grm.bin```, ```PREFIX.grm.N.bin``` and ```PREFIX.grm.id```)|
| --pheno PHENO|  REQUIRED. The name of phenotype file, following GCTA phenotype file format (space delimited, no column names). The first two columns are FID and IID and phenotypes start from the third column. |
| --mpheno m| OPTIONAL. If you have multiple phenotypes in the file, you can specify by ```--mpheno m```. Otherwise, the first phenotype will be used.|
| --PC PC| OPTIONAL. The name of PCs file, following GCTA (space delimited, no column names) ```--pca``` file (same as plink ```--pca```). The third column is the first PC, the forth column is the second PC...|
| --npc n| OPTIONAL. You can specify top n PCs to be adjusted by ```--npc n```. Otherwise, all PCs in the PC file will be used.|
| --covar COVAR| OPTIONAL. The name of covariate file, following GCTA ```--qcovar``` file format or .csv file format. It may contain sex, age, etc. *Note that this file does not include principal components, which need to be include seperately by ```--PC PC```*.|
| --k k| OPTIONAL. You can specify the number of rows in restoring the GCTA GRM binary file into matrix each time. If not provide, it will process the whole GRM at one time. When you have a relative large sample size, specifying ```--k k``` can speed up the computation and save the memory. |
| --std | OPTIONAL. Run SAdj-HE by specifying ```--std```. Otherwise, UAdj-HE will be computed.  (There are potential bugs with the standardized version, so it is reccommended to use unstandardized for now).|



The output of ```AdjHE.py``` contains heritability estimation and its standard error in a .csv file. Computational time and peak memory are also provided.


```
module load python
prefix='EXAMPLE/grm' #We partitioned the GRM into 200 parts.
id='/EXAMPLE/ukbiobank.grm.id'
pheno='Example/phen.pheno'
PC='Example/pcas.eigenvec'
covar='Example/covar.csv'
out='/PATH/TO/RESULT/DIRECTORY'

python AdjHE.py --prefix ${prefix} --PC ${PC} --npc 10  --covar ${covar} --pheno ${pheno} --mpheno 1 --out ${out}
```

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

### PBS script example 
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
