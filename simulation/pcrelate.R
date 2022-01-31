library(GWASTools)
library(GENESIS)
library(SNPRelate)
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
job = as.numeric(Sys.getenv("PBS_ARRAYID"))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied.")
##supply default values
}else{
    for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
        }
}
bed.fn <- paste0(prefix,'.bed')
fam.fn <- paste0(prefix,'.fam')
bim.fn <- paste0(prefix,'.bim')
pheno <- read.table(paste0(prefix,'.pheno'))

gdsfile <- paste0(prefix,'.gds')
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gdsfile, family=TRUE,cvt.chr="int", cvt.snpid="auto", verbose=FALSE)
gds <- GdsGenotypeReader(gdsfile)
gds_genoData <-  GenotypeData(gds)
genoIterator <- GenotypeBlockIterator(gds_genoData, snpBlock=100)
ptm <- proc.time()
mypcair <- pcair(gds)
proc.time() - ptm

ytPcy_func<-function(number_pc){
  summand = 0
  if(number_pc == 0){summand = 0}
  else for(i in 1:number_pc){
      PC_iy = crossprod(grm_svd$u[,i],y)
      summand = summand + (grm_svd$d[i]-1)*PC_iy^2-(grm_svd$d[i]-1) }
  return(summand)
}

eig2<-function(number_pc){
      summand = 0
        if(number_pc == 0){summand = 0}
          else for(i in 1:number_pc){
           summand = summand + (grm_svd$d[i]-1)^2 }
    return(summand)
}

naive_he <- function(y){
    tn = sum(y)**2/npeople1
    sigg = npeople1*ytAy - trace_A*yty
        sigg = sigg-ytAy+tn*trace_A # add 1's
        sige = trace_A2*yty - trace_A*ytAy
        sige = sige-tn*trace_A2 # add 1's
    sigg/(sigg+sige)
  }

for(i in 3:npc){
ptm <- proc.time()
mypcrel <- pcrelate(genoIterator, pcs = mypcair$vectors[,1:i,drop=FALSE],training.set = mypcair$unrels)
GRM <- pcrelateToMatrix(mypcrel)
npeople1=nrow(GRM)
IndivID=c(1:npeople1)
IndivID=as.character(IndivID)
GRM1=GRM[IndivID,IndivID]
GRM = as.matrix(GRM1)
res_y = lm(pheno[,3]~mypcair$vectors[,1:i])$residuals
npeople1 = length(res_y)

trace_A2 = sum((GRM*GRM))
trace_A = sum(diag(GRM))
#y = scale(res_y)
y = res_y
ytAy = t(y)%*%GRM%*%y
yty = sum(y^2)
h2 = naive_he(res_y)

print(proc.time() - ptm)
write.table(h2,paste0('pcrelate_',job),row.names=FALSE,col.names=FALSE,quote=F,append=T)

}
