#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=coffm049@umn.edu
#SBATCH -o sim_fit.out
#SBATCH -e sim_fit.err

module load python3

cd ~/tools/Basu_herit

python ~/tools/Basu_herit/Simulate_genos.py 
