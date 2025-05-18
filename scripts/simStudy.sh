#!/bin/bash

for i in {1..100}   # 100 simulations
do
    conda activate datasci
    python -m Simulate --argfile Simulate/test.json
    plink2 --bfile temp/simulation --glm --out temp/simulation --covar temp/simulation.covar
    awk -F'\t' 'NR==1 || $7 == "ADD" {print}' temp/simulation.PHENO1.glm.linear > temp/simulation.PHENO1.glm.linear.filtered
    awk '{$7=""; print $0}' temp/simulation.fam > temp/simulation.fam2
    mv temp/simulation.fam2 temp/simulation.fam
    PRSice_linux --base temp/simulation.PHENO1.glm.linear.filtered --score avg --target temp/simulation,temp/simulation.fam \
     --out temp/simulation --A1 REF --stat BETA --snp ID --bp POS --pvalue P --keep-ambig --cov  temp/simulation.covar --cov-col subj_ancestries --cov-factor  subj_ancestries
    Rscript Estimate/estimators/ppmxPRS.R
done

