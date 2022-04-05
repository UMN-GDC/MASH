import argparse
from argparse import RawTextHelpFormatter

import timeit

start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running adjusted HE regression',description='This program gives estimation in formula fashion.\n Make sure you have enough memory to store GRM matrix in python.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files.\nIf --npc is not specified then all PCs in the file will be used.')
parser.add_argument('--npc', nargs="+", type=int, help='Specify the number of PCs to be adjusted')
parser.add_argument('--prefix', type=str, help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]',required=False)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=False)
parser.add_argument('--mpheno',nargs="+",type=int, default=1,help='Specify which phenotype to use from phenotype file (Can be a list)')
parser.add_argument('--k',type=int,help='Specify the number of rows in restoring the GRM each time.\n This could affect the computation time and memory especially when sample size is large. If not provide, it will process the whole GRM at one time.')
parser.add_argument('--out',type=str, help='Specify the output file name. [required]',required=False)
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE (i.e., with standardization)')
parser.add_argument('--covars',nargs="+",type=int, default=1,help='Specify which covariates to control for from the covariate file. Should be a list of the column numbers not including the FID and IID columns')
parser.set_defaults(PC="None", npc=0, covar="None", mpheno=1, k=0, prefix = "None", pheno = "None", out = "None")

# Or accept a file with all arguments
parser.add_argument("--argfile", default=None, type=str, help="Filename to be passed containing all information for PC's, covariates, phenotypes, and grm")
# Set defaults

args = vars(parser.parse_args())


#args['f'] = '/home/christian/Research/Stat_gen/tools/Basu_herit/Example/Arg_file.txt'
if args['argfile'] != None :
    d= {}
    with open(args['argfile']) as f:
        for line in f:
           (key, val) = line.split("=")
           # remove line break
           d[key] = val[:-1]
    args.update(**d)
    
# Ensure types 
args["k"] = int(args["k"])

try:
    # convert a string of arugments sto a list
    args['mpheno'] = eval(args['mpheno'])
except:
    # Convert a single integer value to a list
    args['mpheno'] = list(args['mpheno'])
  
try:
    # convert a string of arugments sto a list
    args['npc'] = eval(args['npc'])
except:
    # Convert a single integer value to a list
    args['npc'] = list(args['npc'])
    

    
## Do the same for specified covariates
try:
    # Convert to    list of integers of agrument is a list
    args["covars"] = eval(args["covars"])
except: 
    # convert single integer to integer list 
    args["covars"] = args["covars"]


# try:
#     args["npc"] = int(args["npc"])
# except:
#     args["npc"] = eval(args["npc"])

print(args)



