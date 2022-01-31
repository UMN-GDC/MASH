#library(genetics)
library(data.table)
library(snpStats)
library(popkin)
library(reshape2)

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied.")
##supply default values
    nclus = 4
    npeople = 1000
    heritability = 0.8
    a1 = 0
    a2 = 0
    a3 = 0 
    a4 = 0
    F1=0.15
    F2=0.15
    F3=0.15
    F4=0.15
    nmarkers = 15000
    dir = getwd()
}else{
    for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
        }
}
### Simulation for homogeneous Pop####
job = as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(seed) #for GRM




freq1A=rep(0,nmarkers)
freq2A=rep(0,nmarkers)
freq3A=rep(0,nmarkers)
freq4A=rep(0,nmarkers)
marker1A=matrix(0,nrow=npeople,ncol=nmarkers)
marker2A=matrix(0,nrow=npeople,ncol=nmarkers)
marker3A=matrix(0,nrow=npeople,ncol=nmarkers)
marker4A=matrix(0,nrow=npeople,ncol=nmarkers)
Pop1=matrix(0,nrow=npeople,ncol=nmarkers)
Pop2=matrix(0,nrow=npeople,ncol=nmarkers)
Pop3=matrix(0,nrow=npeople,ncol=nmarkers)
Pop4=matrix(0,nrow=npeople,ncol=nmarkers)
p0=runif(nmarkers,0.1,0.9)
p1=rbeta(nmarkers,(1-F1)/F1*p0,(1-F1)/F1*(1-p0))
p2=rbeta(nmarkers,(1-F2)/F2*p1,(1-F2)/F2*(1-p1))
p3=rbeta(nmarkers,(1-F3)/F3*p2,(1-F3)/F3*(1-p2))
p4=rbeta(nmarkers,(1-F4)/F4*p3,(1-F4)/F4*(1-p3))
freq1A=p1
freq2A=p2
freq3A=p3
freq4A=p4
  for(i in 1:length(freq1A)){marker1A[,i]=rbinom(npeople,2,freq1A[i]); Pop1[,i]= (marker1A[,i] -2*freq1A[i])/(sqrt(2*freq1A[i]*(1-freq1A[i])))}
  for(i in 1:length(freq2A)){marker2A[,i]=rbinom(npeople,2,freq2A[i]);Pop2[,i]= (marker2A[,i] -2*freq2A[i])/(sqrt(2*freq2A[i]*(1-freq2A[i])))}
  for(i in 1:length(freq3A)){marker3A[,i]=rbinom(npeople,2,freq3A[i]);Pop3[,i]= (marker3A[,i] -2*freq3A[i])/(sqrt(2*freq3A[i]*(1-freq3A[i])))}
  for(i in 1:length(freq4A)){marker4A[,i]=rbinom(npeople,2,freq4A[i]);Pop4[,i]= (marker4A[,i] -2*freq4A[i])/(sqrt(2*freq4A[i]*(1-freq4A[i])))}
  
  genotype =rbind(marker1A,marker2A,marker3A,marker4A) 
  pgenotype<-apply(genotype,2,function(x){
        sum(x)/(2*length(x)) 
  })
  out <- which(pgenotype<0.05 | pgenotype>0.95)
  remained <- setdiff(1:nmarkers,out) 

    set.seed(job) #for phenotype
  ncausalsnps =   pcausalsnps * length(remained)
  samplecausalsnps=sample(remained,ncausalsnps)
#  heritability= as.numeric(Args[8])
  
  CausalGeno1=marker1A[,samplecausalsnps]
  CausalGeno2=marker2A[,samplecausalsnps]
  CausalGeno3=marker3A[,samplecausalsnps]
  CausalGeno4=marker4A[,samplecausalsnps]
  
  CausalFreq1=freq1A[samplecausalsnps]
  CausalFreq2=freq2A[samplecausalsnps]
  CausalFreq3=freq3A[samplecausalsnps]
  CausalFreq4=freq4A[samplecausalsnps]
  
  beta = rnorm(length(samplecausalsnps),0,sqrt(heritability/(length(samplecausalsnps)))) 
  
  sigg1 = var(CausalGeno1%*%beta)
  sigg2 = var(CausalGeno2%*%beta)
  sigg3 = var(CausalGeno3%*%beta)
  sigg4 = var(CausalGeno4%*%beta)
  sige1 = sigg1/heritability-sigg1
  sige2 = sigg2/heritability-sigg2
  sige3 = sigg3/heritability-sigg3
  sige4 = sigg4/heritability-sigg4
  c(sige1,sige2,sige3,sige4)
  c(sigg1,sigg2,sigg3,sigg4)

#  a1 = 0
#  a2 = 0.1
#  a3 = 0.2
#  a4 = 0.3
  phenotype1=CausalGeno1%*%beta + rnorm(npeople,0,sqrt(sige1)) + a1
  phenotype2=CausalGeno2%*%beta + rnorm(npeople,0,sqrt(sige2)) + a2
  phenotype3=CausalGeno3%*%beta + rnorm(npeople,0,sqrt(sige3)) + a3
  phenotype4=CausalGeno4%*%beta + rnorm(npeople,0,sqrt(sige4)) + a4
 
  c(var(phenotype1),var(phenotype2),var(phenotype3),var(phenotype4))
if(nclus==4){
  phenotype=c(phenotype1,phenotype2,phenotype3,phenotype4)
  genotype=rbind(marker1A[,remained],marker2A[,remained],marker3A[,remained],marker4A[,remained])
  if(job==101){
snpmatrix=new('SnpMatrix',genotype)
stratum=rep(1:nclus,each=npeople)
i1=1:npeople
i2=(npeople+1):(2*npeople)
i3=(2*npeople+1):(3*npeople)
i4=(3*npeople+1):(4*npeople)
f12=Fst(snpmatrix[c(i1,i2),],stratum[c(i1,i2)])
f13=Fst(snpmatrix[c(i1,i3),],stratum[c(i1,i3)])
f14=Fst(snpmatrix[c(i1,i4),],stratum[c(i1,i4)])
f32=Fst(snpmatrix[c(i3,i2),],stratum[c(i3,i2)])
f42=Fst(snpmatrix[c(i4,i2),],stratum[c(i4,i2)])
f34=Fst(snpmatrix[c(i3,i4),],stratum[c(i3,i4)])
print(weighted.mean(f12$Fst,f12$weight,na.rm=T))
print(weighted.mean(f13$Fst,f13$weight,na.rm=T))
print(weighted.mean(f14$Fst,f14$weight,na.rm=T))
print(weighted.mean(f32$Fst,f32$weight,na.rm=T))
print(weighted.mean(f42$Fst,f42$weight,na.rm=T))
print(weighted.mean(f34$Fst,f34$weight,na.rm=T))
  }
}
if(nclus==2){
  phenotype=c(phenotype1,phenotype2)
genotype=rbind(marker1A,marker2A)
#snpmatrix=new('SnpMatrix',genotype)
#stratum=rep(1:nclus,each=npeople)
#f1=Fst(snpmatrix,stratum)
#print(weighted.mean(f1$Fst,f1$weight,na.rm=T))
}

# c(var(phenotype1),var(phenotype2),var(phenotype3),var(phenotype4))/var(phenotype)
c(sigg1/var(phenotype1),sigg2/var(phenotype2),sigg3/var(phenotype3),sigg4/var(phenotype4))
#  std_y = scale(phenotype)
  
  ###create ped and map file for plink
  npeople1=npeople*nclus
  IndivID=c(paste(1:npeople1,sep=""))
  FamID=c(paste(1:npeople1,sep=""))
  FatherID=rep(0,npeople1)
  MotherID=rep(0,npeople1)
  Sex=c(rep(1,npeople1/2),rep(2,npeople1/2))
  #phenotype=rep(-9,npeople1)
  #plink.ped=data.frame(FamID,IndivID, FatherID,MotherID, Sex,phenotype, Genotype_data)
  plink.fam=data.frame(FamID,IndivID, FatherID,MotherID, Sex,phenotype)
  gp.fam=plink.fam
  colnames(gp.fam) <- c('FID','IID','dad','mom','sex','pheno')
  gp.fam$population = rep(1:nclus,each=npeople)
  #box=data.frame(phenotype,pop=c(rep('pop1',npeople),rep('pop2',npeople),rep('pop3',npeople),rep('pop4',npeople)))
  #yp<-summary(aov(box$phenotype~box$pop))[[1]]$`Pr(>F)`[1]
  #write.table(yp,'yp.txt',row.names=F,col.names=F,append=T)
  snpnames=c(paste("snp",1:ncol(genotype),sep=""))
  chrom=rep(1,ncol(genotype))
  bpdis=sort.default(sample(1:902606,ncol(genotype)))
  plink.map=data.frame(chrom,snpnames,bpdis)
  
  # transpose the allele dosages
  # use as.matrix to ensure that we (efficiently) transpose the allele dosage numbers and nothing else
  tgeno = t(as.matrix(genotype))
destfile=paste('fst',seed,'.pdf',sep='')  
  if (job==1){
    stratum <- rep(1:nclus,each=npeople)
    kinship <- popkin(tgeno,stratum)
    pairwise_fst <- pwfst(kinship) # estimated matrix
    pdf(destfile)
leg_title <- expression(paste('Pairwise ', F[ST])) # fancy legend label
plot_popkin(pairwise_fst, labs = stratum, leg_title = leg_title)
dev.off()
}
  # build the dosage file itself
  geno.dosage = data.table(cbind(paste0("snp", 1:ncol(genotype)), "A", "T", tgeno))
  nocausal.dosage = geno.dosage[-samplecausalsnps,]
  fwrite(geno.dosage, file = paste0("sim",job,".dosage"), quote = F, sep = "\t", col.names=F)
  #fwrite(nocausal.dosage, file = "nocausal.dosage", quote = F, sep = "\t", col.names=F)
  # could add header, this one should work with --id-delim "-" in PLINK call
  #colnames(ceu.dosage) = c("SNP", "A1", "A2", paste0("id", 1:nrow(ceu), "-id", 1:nrow(ceu)))
  
  # write the dosage file to disk
  # set col.names=T if you want a header on the file  
  write.table(plink.fam[,c(1,2,6)],paste0('sim',job,'.pheno'),quote=F,col.names=F,row.names=F)   
  write.table(plink.fam,paste0('sim',job,'.fam'),quote=F,col.names = F,row.names = F)
  write.table(plink.map,paste0('sim',job,'.map'),quote=F,col.names = F,row.names = F)
  write.table(t(c(var(phenotype1),var(phenotype2),var(phenotype3),var(phenotype4))/var(phenotype)),paste('varYratio',seed,job,sep='_'),append=T,row.names=F,col.names=F,quote=F)
  write.table(t(c(var(phenotype1),var(phenotype2),var(phenotype3),var(phenotype4))),paste('varY',seed,job,sep='_'),append=T,row.names=F,col.names=F,quote=F)

      
