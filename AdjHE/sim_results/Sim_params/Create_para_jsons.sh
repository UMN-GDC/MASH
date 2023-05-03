#!/bin/bash -l        

module load python3

# Make json files for parallel simulations

python ../../simulation_parallelizer.py --argfile Covariates_non_random_sites.json 


