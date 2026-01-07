import numpy as np
import pandas as pd

# Make a function that simulates the outcome "Z" as a linear combination of all columns starting with "Y" from the dataframe df

def sim_outcome(df, beta):
    # Extract the columns starting with "Y"
    Y_cols = [col for col in df.columns if col.startswith("Y")]
    
    # Simulate the outcome
    df["Z"] = np.dot(df[Y_cols], beta) + np.random.normal(0, 1, df.shape[0])
    
    return df