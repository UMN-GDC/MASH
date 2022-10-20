#!/bin/bash -l
#SBATCH --time=8:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH -o full_analysis.out
#SBATCH -e full_analysis.err

# 2 hour wall time for random effect

module load python

cd ~/tools/Basu_herit

python AdjHE.py --argfile simulations/Analysis_jsons/ABCD/basic.json
python AdjHE.py --argfile simulations/Analysis_jsons/ABCD/basic_sites.json
python AdjHE.py --argfile simulations/Analysis_jsons/ABCD/site_rv.json  
python AdjHE.py --argfile simulations/Analysis_jsons/ABCD/MOM.json

