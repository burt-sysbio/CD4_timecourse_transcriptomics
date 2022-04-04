import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style = "ticks", context = "talk")

dir1 = "../../data/correlation/hclust_output/hclust_corr_thm.csv"
dir2 = "../../data/correlation/hclust_output/hclust_corr_th0.csv"

baseline = "degs_Th1_Th2_cor0.3"
dir3 = "../../data/linear_model/linear_model_processed_" + baseline + ".csv"

df_clust_thm = pd.read_csv(dir1)
df_clust_th0 = pd.read_csv(dir2)

df_clust_thm["model"] = "Th1/2"
df_clust_th0["model"] = "Th0"

df_clust = pd.concat([df_clust_th0, df_clust_thm])

df_linmo = pd.read_csv(dir3)
#
# # from linear model only select relevant variables
# df_linmo = df_linmo[["gene", "category"]]
# df_clust = df_clust[["gene", "category"]]
#
#
df_linmo["type"] = "linear model"
df_clust["type"] = "hclust"
#
df = pd.concat([df_clust, df_linmo])

df.loc[df["category"]=="Ind","category"] = "Independent"

df2 = df.groupby(["model", "type"])["category"].value_counts(normalize = True).mul(100).rename('percent').reset_index()

for x in ["Th0", "Th1/2"]:
    df_bar = df2.loc[df2["model"]== x]
    g = sns.catplot(data = df_bar, x = "category", y = "percent", hue = "type",
                    kind = "bar",
                    palette= ["0.2", "0.8"], aspect = 1.0,
                    legend=False)

    g.set(xlabel = "", ylabel = "genes (%)")
    g.set_xticklabels(rotation=90)
    sns.despine(top = False, right = False)

    plt.legend(loc = "upper right")
    plt.show()

    if x == "Th1/2":
        x = "ThMix"
    g.savefig("../../figures/fig4/barplots/barplot_" + x + baseline + ".svg")
    g.savefig("../../figures/fig4/barplots/barplot_" + x + baseline + ".pdf")


# plot individually hclust and linear model
sns.set(context = "poster", style = "ticks")
for x in ["hclust", "linear model"]:
    df_bar = df2.loc[df2["type"] == x]
    df_bar = df_bar.loc[df_bar["model"] == "Th1/2"]

    g = sns.catplot(data = df_bar, x = "category", y = "percent",
                    kind = "bar", palette= ["0.1"], aspect = 0.8,
                    legend=False)

    g.set(xlabel = "", ylabel = "genes (%)")
    g.set_xticklabels(rotation=90)
    sns.despine(top = False, right = False)

    plt.show()

    g.savefig("../../figures/fig4/barplots/barplot_th12_" +x + baseline + ".svg")
    g.savefig("../../figures/fig4/barplots/barplot_th12_" +x+ baseline + ".pdf")