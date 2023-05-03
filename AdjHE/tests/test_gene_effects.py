import os
import numpy as np
import pytest
import pandas as pd
from pandas_plink import read_plink
from AdjHE.simulation_helpers.sim_gene_effects import sim_gen_effects
from AdjHE.estimation.all_estimators import Basu_estimation

rng = np.random.default_rng(12345)

def test_sim_gen_effects() : 
    # read in 1kg dataset from examples (chr1)
    (bim, fam, bed) = read_plink("AdjHE/examples/Input_files/1kg_data/chr1")
    df = pd.read_table("AdjHE/examples/Input_files/1kg_data/igsr_samples.tsv")[["Sample name", "Population code"]]

    Gene_contrib = sim_gen_effects(rng = rng, genotypes = bed.T, df = df, causals= None, prop_causal = 0.1, site_dep = False, variance_propto_frequency = False)


    


