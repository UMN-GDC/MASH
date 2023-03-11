#!/bin/bash -l        
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=sample_email@umn.edu 

cd ~/PATH/TO/DIR/
module load python 
python AdjHE.py --argfile Example/Arg_file.txt
