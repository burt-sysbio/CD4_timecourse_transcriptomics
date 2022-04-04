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

sns.set(context = "poster", style = "ticks")
df = pd.read_csv("../../data/peine_logfc_tidy.csv")

nclust = 8
fname = "peine_th1_th2_deg_" + str(nclust) + "_clusters.csv"
df_clust = pd.read_csv(fname)

# keep only one cluster
clust = "c1" 
df_clust = df_clust[["gene", clust]]

# keep only genes found in cluster annot
df = df.loc[df.gene.isin(df_clust.gene),:]

# merge dataframes
df = pd.merge(df, df_clust, how = "left", on = "gene")
# get number of genes in cluster

df2 = df[["gene", clust]].drop_duplicates()
counts = df2.groupby([clust]).count()

mypal = ["tab:blue", "tab:purple", "tab:red", "grey"]
g = sns.relplot(data = df, x = "time", y = "value", col = "c1", hue = "celltype", kind = "line",
                palette = mypal, hue_order = ["Th1", "ThMix", "Th2", "Th0"],
                facet_kws = {"sharey" : False}, col_wrap = 3, aspect = 1.2)

clusts = [str(i+1) for i in range(nclust)]
for ax, n, c in zip(g.axes.flat, counts.values, clusts):
    title = r"c" + c + ", " + str(n[0]) + " genes"
    ax.set_title(title)
    
g.set(xlabel = "time (h)", ylabel = "expr. norm.", xlim = (0, df.time.max()))
sns.despine(top = False, right = False)
plt.show()

sdir = "../../figures/fig2/msigpro_timecourse_" + str(nclust) + "_clusters"
g.savefig(sdir + ".svg")
g.savefig(sdir + ".pdf")
