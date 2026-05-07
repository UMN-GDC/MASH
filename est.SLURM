#!/bin/bash
#SBATCH --job-name=subcorts
#SBATCH --output=logs/subcors.log
#SBATCH --mem=32G
#SBATCH --time=1:00:00  # Give the master enough time to finish the whole pipe

source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate MASH

MASH --argfile argfile2.json
