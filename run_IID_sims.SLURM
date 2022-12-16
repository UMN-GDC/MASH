#!/bin/bash -l        
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-4 




cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu2


module load python3

python Simulator.py --argfile Simulations/IID_Sim.json --out IID_sims/2000/$SLURM_ARRAY_TASK_ID


