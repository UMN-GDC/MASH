# Basu_herit

# adjustedHE

## Adjusted-HE with closed form formula version

```AdjHE_formula.py```  estimates SNP-heritability via closed form formula with single GRM as input.

Please check the input description with ```./AdjHE_formula.py --help```.

```AdjHE_formula.py``` can take a single GRM, a covariate file, a PC file and a phenotype file (with multiple phenotypes) as input files.
1. ```--prefix``` follows GCTA binary GRM format. (```PREFIX.grm.bin```, ```PREFIX.grm.N.bin``` and ```PREFIX.grm.id```)
2. ```--pheno``` follows GCTA phenotype file format. The first two columns are FID and IID and phenotypes start from the third column. If you have multiple phenotypes in the file, please specify by ```--mpheno``` (DEFAULT is 1, i.e., the third column of ```--pheno```)
3. ```--covar``` follows GCTA ```--qcovar``` file format. It may contain sex, age, etc. *Principal components should be included in a seperate file by ```--PC```*.
4. ```--PC``` follows GCTA ```--pca``` file (same as plink ```--pca```). The third column is the first PC, the forth column is the second PC... You can specify the number of PCs to be adjusted by ```--npc``` (DEFAULT is all PCs in the file). *Also make sure these PCs are computed based on the same set of individuals in GRM.*
5. You can use ```--std``` to compute Standardised-Adj-HE. Default is to compute Unstandardized-Adj-HE.
6. You can use ```--k``` to specify the number of rows in restoring the GRM (restore GCTA GRM format into matrix) each time. If not provide, it will process the whole GRM at one time. When you have a relative large sample size, controlling *k* can speed up the computation.

The output of ```AdjHE_formula.py``` contains heritability estimation and its standard error. Computational time and peak memory are also provided.

## Adjusted-HE with regression version

For large sample size (e.g. biobank size), you can use ```AdjHE_reg_s1.py``` and ```AdjHE_reg_s2.py``` to perform linear regression to get the heritability estimation.
This allows you to take part GRM as input and compute efficiently.


#### ```AdjHE_reg_s1.py```

You should FIRST call ```AdjHE_reg_s1.py``` before ```AdjHE_reg_s2.py```. Please check the input description with ```./AdjHE_reg_s1.py --help```.

1. ```--pheno```, ```--covar```, ```--PC``` follow the same format as mentioned in ```AdjHE_formula.py```.
2. ```--prefix``` follows GCTA binary partitioned GRM format. (e.g. ukbiobank.part_200_') Please refer to ```--make-grm-part``` in https://cnsgenomics.com/software/gcta/#MakingaGRM.
3. ```--Npart``` specifies the number of parts (m) you specify in ```--make-grm-part m i``` when constructing GRM by parts.
4. ```--id``` follows GCTA ```PREFIX.grm.id``` format. It contains ids of *ALL* samples. You can create this file by ```cat ukbiobank.part_200_*.grm.id > ukbiobank.grm.id```.
5. ```--job i``` computes the i-th part of GRM in the current run. (e.g., It would take ukbiobank.part_200_i.* as input) You can easily parallel jobs using PBS Job Array by specifying ```--job ${PBS_ARRAYID}```. In our estimation of the UKB data, we partitioned the whole data set (n = 305,639) into 200 parts and allocated 15GB memory to each job.
6. ```--out ``` specifies the directory to store all intermediate results and a Log file. After processing all parts of the GRM, it should generate m intermediate files (e.g., ```pheno1_1```,...,```pheno1_200```). The Log file will provide the computation time for each run.

#### ```AdjHE_reg_s2.py```

You should only call ```AdjHE_reg_s2.py``` once after finishing ```AdjHE_reg_s1.py``` with m part GRM successfully. Please check the input description with ```./AdjHE_reg_s2.py --help```.

1. ```--pheno```, ```--covar```, ```--PC``` follow the same format as mentioned in ```AdjHE_formula.py```. ```--Npart```, ```--id```  follow the same format as mentioned in ```AdjHE_reg_s1.py```
2. ```--out``` specifies the *SAME* directory as you use in ```AdjHE_reg_s1.py```. Make sure you have m intermediate result files in this directoory. The point estimate and standard error of heritability will be output to the Log file generated in ```AdjHE_reg_s1.py```.
# Basu_herit
