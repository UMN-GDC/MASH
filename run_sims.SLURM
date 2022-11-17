#!/bin/bash -l        
#SBATCH --time=8:30:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu

cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit


module load python3
python run_sims.py
