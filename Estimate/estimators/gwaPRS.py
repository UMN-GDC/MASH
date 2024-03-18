import subprocess


def fitPlink2GWAS(prefix):

    command = f"plink2 --bfile {prefix} --glm --pheno {prefix}.pheno --covar {prefix}.covar"

    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    
    pass


