#!/bin/bash -l        
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-10 




cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit


module load python3

python Simulator.py --argfile Simulations/Het/Het.json --out Het/Results/5000/$SLURM_ARRAY_TASK_ID




