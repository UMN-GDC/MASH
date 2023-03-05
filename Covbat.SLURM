#!/bin/bash -l 
#SBATCH --time=1:30:00
#SBATCH --ntasks=1
#SBATCH --mem=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=coffm049@umn.edu


module load python3
cd ~/tools/Basu_herit
python AdjHE.py --argfile Covbat.json
