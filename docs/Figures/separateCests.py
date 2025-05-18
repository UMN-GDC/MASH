from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#%% Basic testing for single site and cluster
rng = np.random.default_rng(123)

def simNEstSeparate(rng = 123, h2Hom = 0.5, h2Het = [0.0, 0.0], N=500): 
    sim = pheno_simulator(rng = rng, nsubjects= N)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = h2Hom, h2Het= h2Het, alpha = 0)
    est = h2Estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    result = est.estimate(mpheno = ["Y0"], npc = [1], Method = "GCTA", fixed_effects= ["Xc"]) 
    simvar= est.df[["homo_contrib", "errors"]].var()
    simest = simvar["homo_contrib"] / np.sum(simvar)

    df1= est.df.query("subj_ancestries==0")
    G1 = est.GRM[df1.index, :][:, df1.index]
    df2 = est.df.query("subj_ancestries==1")
    G2 = est.GRM[df2.index, :][:, df2.index]
    
    est1 = h2Estimation()
    est1.GRM = G1
    est1.df = df1
    simvar1 = est.df[["homo_contrib", "errors"]].var()
    simest1 = simvar1["homo_contrib"] / np.sum(simvar1)
    
    est2 = h2Estimation()
    est2.GRM = G2
    est2.df = df2
    simvar2 = est.df[["homo_contrib", "errors"]].var()
    simest2 = simvar2["homo_contrib"] / np.sum(simvar2)

    result1 = est1.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc"])
    result2 = est2.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc"])

    return pd.DataFrame({"EstCombined" : result["h2"], "Est1": result1["h2"], "Est2": result2["h2"],
                         "SimCombined" : simest, "Sim1" : simest1, "Sim2" : simest2})

reps =50 
df = pd.DataFrame(np.zeros((reps, 6)), columns = ["EstCombined", "Est1", "Est2", "SimCombined", "Sim1", "Sim2"])

for i in range(reps) :
    df.iloc[i,:] = simNEstSeparate(rng = i, h2Hom = 0.25, h2Het = [0.25, 0.25])

# Save the dataframe to csv
df.to_csv("SimProfilesDiffering.csv", index= False, header= True)

#%% Plotting the distributions of the simulations
df2 = pd.melt(frame = df,
        var_name = "Method",
        value_name = "Value")
df2["Sim"] = df2.Method.str.startswith("Sim")


sns.histplot(data = df2.query("Sim == False"), x= "Value", hue = "Method", alpha =0.25)
plt.xlabel("Simulated heritability")
plt.savefig("SimProfilesDiffering.png", dpi = 300)
