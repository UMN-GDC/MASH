#!/bin/bash -l        
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH --array=1-50 



#cd /home/christian/Research/Stat_gen/tools/Basu_herit 
cd /panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit


module load python3

# loop over each file created
sim_list=$(ls Simulations/Sim_params/Temp)
SIM=$(echo $sim_list | cut -d" " -f$SLURM_ARRAY_TASK_ID)


# Run estimation
python Simulate.py --argfile Simulations/Sim_params/Temp/${SIM} --out Simulations/Results/1k25S_5C_Fixed_site_${SLURM_ARRAY_TASK_ID}

# remove that file
rm Simulations/Sim_params/Temp/${SIM}


