#!/bin/bash -l        
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-4 




cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu2


module load python3

python Simulator.py --argfile Simulations/EQUAL_Sim.json --out EQUAL_sims/$SLURM_ARRAY_TASK_ID


