import argparse
from argparse import RawTextHelpFormatter
from functions import *


start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running adjusted HE regression',description='This program gives estimation in formula fashion.\n Make sure you have enough memory to store GRM matrix in python.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files.\nIf --npc is not specified then all PCs in the file will be used.')
parser.set_defaults(PC="NULL")
parser.add_argument('--npc', type=int, help='Specify the number of PCs to be adjusted')
parser.set_defaults(npc=0)
parser.add_argument('--prefix', type=str, help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]',required=True)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',type=int, default=1,help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
parser.add_argument('--k',type=int,help='Specify the number of rows in restoring the GRM each time.\n This could affect the computation time and memory especially when sample size is large. If not provide, it will process the whole GRM at one time.')
parser.set_defaults(k=0)
parser.add_argument('--out',type=str, help='Specify the output file name. [required]',required=True)
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE (i.e., with standardization)')

args = parser.parse_args()


npc = args.npc
outprefix = args.out
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename=outprefix+'.log',filemode='a')
for arg, value in sorted(vars(args).items()):
    logging.info("Argument %s: %r", arg, value)

prefix = args.prefix

G = ReadGRMBin(prefix)
ids = G['id']
n_phen_nona = ids.shape[0]


# Read data function that can load csv pheno and txt file types
def read_datas(file_path) :
 if(file_path.split(".")[-1] == "csv"):
  dat = pd.read_csv(file_path)
 elif(file_path.split(".")[-1] == "phen"):
  dat = pd.read_table(file_path, sep = " " )
 elif(file_path.split(".")[-1] == "txt"):
  dat = pd.read_table(file_path, sep = " " )
 # remove the unintentional columns that sometimes happen with phenotype and csv filetypes
 dat = dat[dat.columns.drop(list(dat.filter(regex='Unnamed')))]
 # only keep intersections of the individuals in the GRM and in the covariates and phenotypes
 intersection_indiv = np.intersect1d(ids.iloc[:,1], dat.IID)
 dat = dat[dat.IID.isin(intersection_indiv)]
 # store it as an array for efficiency and for certain linear alg functions to work
 dat = dat.drop(["FID", "IID"], axis = 1)
# dat = np.asarray(dat.drop(["FID", "IID"], axis = 1))
 return(dat)

# seed the covariates matrix with a column of 1's for the intercept
cov_selected = np.ones(n_phen_nona)

# load phenotypes and covariates
y = read_datas(args.pheno)

# read in covariates if nonnull
if (args.covar != "NULL"):
 X = read_datas(args.covar)

# stack the covariates onto the incercepts
cov_selected = np.column_stack((cov_selected,X))

# onlyt load pcs if non null
if (args.PC != "NULL"):
    PCs = pd.DataFrame(np.loadtxt(args.PC))
    PCs.index = PCs.iloc[:,0].astype("int32")
    if (args.npc == -9):
        npc = PCs.shape[1] - 2
    if (npc != 0):
        final_PC = PCs.loc[intersection_indiv]
        final_PC = final_PC.values[:,2:(2+npc)]
        cov_selected = np.column_stack((cov_selected,final_PC))

# y = y.iloc[:,args.mpheno+1]


# only regress out covariates if they are entered
res_y = regout(cov_selected, y)
# this is reshaping 1-D array to vector in numpy, this might cause problems for multivariate regression
res_y = np.reshape(np.asarray(y), -1)


start_read = timeit.default_timer()
nmarkers = G['N']
x = G['diag'].astype('float64')
n_phen_nona = G['diag'].size
GRM_array_nona = np.zeros((n_phen_nona,n_phen_nona))
temp_i = 0
temp  = 0
k = args.k

if(k == 0):
    k = n_phen_nona

l = list(range(k,n_phen_nona,k))
l.append(n_phen_nona)
for i in l:
    cor = multirange(range(temp_i,i))
    GRM_array_nona[cor['b'],cor['a']] = G['off'][temp:temp+len(cor['b'])]
    GRM_array_nona.T[cor['b'],cor['a']] = G['off'][temp:temp+len(cor['b'])]
    temp = temp + len(cor['b'])
    del(cor)
    temp_i = i

GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']
logging.info('GRM matrix restored done. It takes: '+str(timeit.default_timer() - start_read)+' seconds.')


trace_A = np.sum(x)
trace_A2 = 2*np.sum(G['off'].astype('float64')**2) + np.sum(x**2)
#print(trace_A)
#print(trace_A2)
del(G)

# print(res_y.head())

if (args.std == True):
    h2,se = myformula1(A=GRM_array_nona,y=res_y,trA=trace_A,trA2=trace_A2, npc = npc)
else:
    h2,se= myformula2(A=GRM_array_nona,y=res_y,trA=trace_A,trA2=trace_A2, npc =npc)
logging.info('h2: '+str(h2))
logging.info('Standard error: '+str(se))
logging.info('It takes: '+str(timeit.default_timer() - start_time)+' seconds.')
logging.info('Memory usage:'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
