#!/bin/bash

# TODO: 
# x Make sure risk groups affect outcome
# x make confounders affect endophenotypes
#   x covary with ancestry
#   x covary with risk group
# - simulate outcome from the endos

conda activate datasci

python -m Simulate --argfile Simulate/test.json

# calculate fist 10 PCs using plink2
plink2 --bfile temp/simulation --pca 10 --out temp/simulation

plink2 --bfile temp/simulation --glm --out temp/simulation --covar temp/simulation.covar --pheno temp/simulation.pheno --pheno-name Y0,Y1,Y2,Y3,Y4 

for i in $(seq 0 4) ; 
do
    awk -F'\t' 'NR==1 || $7 == "ADD" {print}' temp/simulation.Y${i}.glm.linear > temp/simulation.Y${i}.glm.linear.filtered
    #awk '{$7=""; print $0}' temp/simulation.fam > temp/simulation.fam2
    #mv temp/simulation.fam2 temp/simulation.fam
    PRSice_linux --base temp/simulation.PHENO1.glm.linear.filtered --score avg --target temp/simulation,temp/simulation.fam --pheno temp/simulation.pheno \
     --out temp/simulation.Y${i} --A1 REF --stat BETA --snp ID --bp POS --pvalue P --keep-ambig --cov  temp/simulation.covar --cov-col subj_ancestries --cov-factor  subj_ancestries \
     --pheno-col Y${i}
done   


Rscript Estimate/estimators/ppmxPRS.R
