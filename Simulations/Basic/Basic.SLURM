#!/bin/bash -l        
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-10 




cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit


module load python3

python Simulator.py --argfile Simulations/Basic/Basic.json --out Basic/Results/2000/$SLURM_ARRAY_TASK_ID




