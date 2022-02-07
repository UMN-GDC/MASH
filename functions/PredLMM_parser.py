import timeit
import logging
import resource
import argparse
from argparse import RawTextHelpFormatter

#-----------------------Set up user inputs for parsing----------------------------------------
start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running PredLMM',description='In a general scenario, one should compute the GRM files first with the PLINK Binary files using GCTA software. Then, follow the jupyter notebook: PredLMM_notebook for estimating heribaility and variance of a phenotype adjusting the availble covariates.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--prefix', type=str, help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]',required=True)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',nargs="+",type=int, default=1,help='Specify which phenotype to use from phenotype file (Can be a list)')
parser.set_defaults(mpheno=1)
parser.add_argument('--out',type=str, help='Specify the output file name. [required]',required=True)
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE (i.e., with standardization)')

args = parser.parse_args()

prefix= args.prefix 
covar = args.covar
pheno = args.pheno
outprefix = args.out
mpheno = args.mpheno


