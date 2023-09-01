# make sure your working directory is MASH
from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest


# The dataframe needs to have the subject ancestries contained in a column called "subj_ancestries"

rng = np.random.default_rng(123)
# Instantiate a simulator If the simulation is comprised of multiple ancestries, put a filepath containing the FID, IID in the same order as the genotpes, and a column with the ancestries
sim = pheno_simulator(rng = rng, plink_prefix="Estimate/examples/Genotype_files/geno", subjAncestries= None)
# Detemrine how many sites, genetic clusters you want
sim.full_sim(nclusts=2, nsites=2)

