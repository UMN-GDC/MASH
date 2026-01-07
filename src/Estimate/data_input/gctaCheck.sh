#/bin/bash

PC=/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/Eigens/no_rels.eigenvec
covar=/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/02_Phenotypes/Covars.tsv
prefix=/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/no_rels/no_rels
pheno=/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/02_Phenotypes/asegs_no_names.phen
out=/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/height_AdjHERE.csv
mpheno=[anthro_height_calc_x]
fixed_effects=[age female household_income]

gcta64  --grm $prefix --pheno $pheno --mpheno 1 --reml --out outname
    --qcovar  /panfs/jay/groups/31/rando149/coffm049/ABCD/Results/02_Phenotypes/Quant_covars.tsv
    #--covar /panfs/jay/groups/31/rando149/coffm049/ABCD/Results/02_Phenotypes/Disc_covars.tsv