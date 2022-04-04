import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
sns.set(context = "talk", style = "ticks")

# load data (control and df where genes with cor>0.6 removed)
df1 = pd.read_csv("../../data/bifurcation_timecourse/bifurcation_timecourse_thres1.csv")
df2 = pd.read_csv("../../data/bifurcation_timecourse/bifurcation_timecourse_thres0.6.csv")

labels = ["original","correlation filter"]

df1["filter"] = "original"
df2["filter"] = "correlation filter"
df_list = [df1,df2]

# make tidy, combine and select first 6 PCs
df_list = [df.melt(id_vars = ["time", "cell", "filter"]) for df in df_list]
df3 = pd.concat(df_list)
keep_pcs = ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]
df3 = df3.loc[df3.variable.isin(keep_pcs),:]

# compute differences vs Th0
df_th0 = df3.loc[df3["cell"] == "Th0"]
df_th0 = df_th0.rename(columns = {"value" : "control"}).drop(columns = ["cell"])
df3 = pd.merge(df3, df_th0, how = "inner")
df3["delta_PC"] = df3["value"] - df3["control"]
df3 = df3.loc[df3["cell"] != "Th0"]

# plotting
mypal = ["tab:blue", "tab:red", "purple"]
g = sns.FacetGrid(df3, col="variable",  row= "filter", hue = "cell",
                  hue_order= ["Th1", "Th2", "Th1/2"], palette= mypal,
                  sharey = True, margin_titles = True)

g.map(sns.scatterplot,"time", "delta_PC")
g.map(sns.lineplot,"time", "delta_PC")

sns.despine(right = False, top = False)
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(xlabel= "time (h)", ylabel = r"$\Delta_{Th0}$ PC score")
g.set(xticks = [0,24,48,72,96,120])
plt.show()

