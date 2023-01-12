#!/bin/bash -l        
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-10 



cd /home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Sim_params 
# cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit/Simulations/Sim_params


#module load python3

# Make json files for parallel simulations

python ../../simulation_parallelizer.py --argfile Add_clusts_sites.json


# loop over each file created
sim_list=$(ls Temp)
=$(echo $sim_list | cut -d" " -f$SLURM_ARRAY_TASK_ID)

rm Temp/*


