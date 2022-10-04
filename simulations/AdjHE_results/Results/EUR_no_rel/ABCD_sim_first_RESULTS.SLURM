#!/bin/bash -l        
#SBATCH --time=0:40:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu 

cd ~/tools/Basu_herit
module load python 
python AdjHE.py --argfile simulations/AdjHE_results/Results/EUR_no_rel/Synth_sims_Argfile.json
