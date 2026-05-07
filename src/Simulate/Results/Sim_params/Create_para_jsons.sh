#!/bin/bash -l        

module load python3

# Make json files for parallel simulations

python ../../hpc_utilities/simulation_parallelizer.py --argfile Covariates_random_sites.json 


