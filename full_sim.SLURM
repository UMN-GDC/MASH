#!/bin/bash -l
#SBATCH --time=6:00:00
#SBATCH --ntasks=4
#SBATCH --mem=22g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH -o full_sim.out
#SBATCH -e full_sim.err

module load python

cd ~/tools/Basu_herit

python ~/tools/Basu_herit/full_sim.py 
