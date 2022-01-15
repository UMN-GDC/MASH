#-----------------------Loading  the required Module------------------------------------------
import timeit
from functions.load_data import *
from functions.math_functions import *

#-----------------------Set up user inputs for parsing----------------------------------------
start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running PredLMM',description='In a general scenario, one should compute the GRM files first with the PLINK Binary files using GCTA software. Then, follow the jupyter notebook: PredLMM_notebook for estimating heribaility and variance of a phenotype adjusting the availble covariates.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--prefix', type=str, help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]',required=True)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',type=int, default=1,help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
parser.add_argument('--out',type=str, help='Specify the output file name. [required]',required=True)
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE (i.e., with standardization)')

args = parser.parse_args()


outprefix = args.out
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename=outprefix+'.log',filemode='a')
for arg, value in sorted(vars(args).items()):
    logging.info("Argument %s: %r", arg, value)

prefix = args.prefix


start_read = timeit.default_timer()

G = ReadGRMBin(prefix)
N = len(G['diag'])
GRM = csr_matrix((N, N));GRM_array = GRM.todense().A1.reshape(N, N)
idx = np.tril_indices(N,-1,N);idy = np.triu_indices(N,1,N);id_diag = np.diag_indices(N)
GRM_array[idx] = G['off'];GRM_array[id_diag] = G['diag'];GRM_array[idy] = GRM_array.T[idy]
GRM_array = np.float32(GRM_array)



#-----------------------convert the GRM to h5py format for faster loading------------------- 
#hf = h5py.File('Data/example_grm.h5', 'w')
#hf.create_dataset('dataset_1', data=GRM_array)
#hf.close()

#-----------------------loading GRM in h5py format------------------------------------------- 
#hf = h5py.File('Data/example_grm.h5', 'r')
#GRM_array= np.array(hf.get('GRM'),dtype="float32")




y = read_datas(args.pheno) 
X = read_datas(args.covar)



#----------------------Knot selection and selecting corresponding vectors----------------------------
subsample_size = 500;
sub_sample = sorted(np.random.choice(range(0,N),subsample_size,replace=False))
non_subsample = np.setdiff1d(range(0,N),sub_sample)
indices = np.hstack((sub_sample,non_subsample))
GRM_array = np.float32(GRM_array[np.ix_(indices,indices)].T)
y = y[indices]; X=X[indices]; X_T = X.T;

G_selected = GRM_array[range(0,subsample_size),:][:,range(0,subsample_size)]
y_sub = y[range(0,subsample_size)]; X_sub=X[range(0,subsample_size)]; X_subT=X_sub.T

logging.info('knot selection and corresponding vectors seleected: '+str(timeit.default_timer() - start_read)+' seconds.')

#------------------Fitting LMM using only the selected subsample (set of knots)-------------------------
A_selc = np.copy(G_selected)-Identity(subsample_size)
result_subsample = derivative_minim_sub(y_sub, X_sub, X_subT, G_selected, A_selc, subsample_size)
# print(result_subsample)


#------------------Running PredLMM----------------------------------------------------------------------
Ct =  np.copy(GRM_array[range(0,subsample_size),:],order='F')
C12 = Ct[:,range(subsample_size,N)]
id_diag = np.diag_indices(N)
diag_G_sub = GRM_array[id_diag]
G_inv = inv(G_selected).T
GRM_array[np.ix_(range(subsample_size,N),range(subsample_size,N))] = sgemm(alpha=1,a=C12.T,b=sgemm(alpha=1,a=G_inv,b=C12))
del G_inv, C12
add = copy(-GRM_array[id_diag] + diag_G_sub) ## diagonal adjustment
np.fill_diagonal(GRM_array, - 1 + diag_G_sub)


result_full = derivative_minim_full(y, X, X_T, Ct, id_diag, add, G_selected, GRM_array, N)
# print(result_full)


final_result = {'GREML (sub) =>':result_subsample,'PredLMM =>': result_full}
# print(final_result)

logging.info('output: '+str(final_result))