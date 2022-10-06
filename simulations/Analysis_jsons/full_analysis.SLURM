#!/bin/bash -l
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH -o full_analysis.out
#SBATCH -e full_analysis.err

module load python

cd ~/tools/Basu_herit

python AdjHE.py --argfile simulations/Analysis_jsons/first_full.json 
