#!/bin/bash -l        

module load python3

# Make json files for parallel simulations

python ../../simulation_parallelizer.py --argfile Add_clusts_sites.json


