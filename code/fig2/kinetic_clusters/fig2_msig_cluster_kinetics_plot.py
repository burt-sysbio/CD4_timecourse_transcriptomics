#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:22:36 2021

@author: burt
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 16:10:41 2021

@author: burt
"""
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

sns.set(context = "poster", style = "ticks", palette = "dark")
df = pd.read_csv("../../data/peine_logfc_tidy.csv")

nclust = 4
fname = "peine_kinetic_deg_" + str(nclust) + "_clusters.csv"
df_clust = pd.read_csv(fname)

# keep only one cluster
clust = "c1" 
df_clust = df_clust[["gene", clust]]

# keep only genes found in cluster annot
df = df.loc[df.gene.isin(df_clust.gene),:]

# merge dataframes
df = pd.merge(df, df_clust, how = "left", on = "gene")
# get number of genes in cluster

df["c1"] = df["c1"].apply(str)


df2 = df[["gene", clust]].drop_duplicates()
counts = df2.groupby([clust]).count()

mypal = ["tab:blue", "tab:grey", "tab:red", "tab:purple"]

g = sns.relplot(data = df, x = "time", y = "value", hue = "c1", col = "celltype", kind = "line", aspect = 1.2)
    
g.set(xlabel = "time (h)", ylabel = "expr. norm.", xlim = (0, df.time.max()))
sns.despine(top = False, right = False)
plt.show()

# only keep one celltype
df = df.loc[df["celltype"] == "Th2", :]
myblack = str(0.2)
mypal = "Greys"
g = sns.relplot(data = df, x = "time", y = "value", hue = "c1", kind = "line", aspect = 1.2, ci = None,
                linewidth = 4, height = 4, legend = False, palette = mypal)
    
g.set(xlabel = "time (h)", ylabel = "expr. norm.", xlim = (0, df.time.max()))
sns.despine(top = False, right = False)
plt.show()

g.savefig("../../figures/fig2/kinetic_genes_clusters.pdf")
g.savefig("../../figures/fig2/kinetic_genes_clusters.svg")

# plot barplot
counts = counts.reset_index()
g = sns.catplot(data = counts, x = "c1", y = "gene", kind = "bar", palette= mypal)
g.set(xlabel = "cluster", ylabel = "genes in cluster")
plt.show()
g.savefig("../../figures/fig2/kinetic_genes_barplot.svg")