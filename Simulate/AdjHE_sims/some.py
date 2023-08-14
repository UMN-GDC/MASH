import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("2k_SWD_COMBAT_0.csv")

# Seaborn boxplot for columns Combat and SWD accumulating all of the reps
# faceted by Sites and Distribution
sns.boxplot(y="SWD", hue="Sites", data=df, palette="Set3")
plt.show()


#%
df.query("Sites == 'S25'")
