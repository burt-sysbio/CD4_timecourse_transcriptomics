import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style = "ticks", context = "paper")

df_tidy = pd.read_csv("../../data/peine_logfc_tidy.csv")

nclust = str(3)

#df_clust = pd.read_csv("../../data/clustered_genes/kinetic_genes_20210714_"+nclust+ "clust.csv")
# new cluster order for three clusters used throughout paper
df_clust = pd.read_csv("../../data/clustered_genes/kinetic_genes_clusters_renamed_20211214.csv")
### provide manual cluster annotation
df_clust["annot_new"] = df_clust["annot"]

if int(nclust)==4:
    cell = "Th2"
    # clust 3 goes to clust 4 and vice versa
    df_clust.loc[(df_clust["celltype"] == cell) & (df_clust["annot"] == 3), "annot_new"] = 4
    df_clust.loc[(df_clust["celltype"] == cell) & (df_clust["annot"] == 4), "annot_new"] = 3

    # reassign annot new to annot
    df_clust["annot"] = df_clust["annot_new"]

if int(nclust)==5:
    #
    df_clust.loc[(df_clust["celltype"] == "Th0") & (df_clust["annot"] == 2), "annot_new"] = 4
    df_clust.loc[(df_clust["celltype"] == "Th0") & (df_clust["annot"] == 4), "annot_new"] = 2

    df_clust.loc[(df_clust["celltype"] == "Th1") & (df_clust["annot"] == 3), "annot_new"] = 4
    df_clust.loc[(df_clust["celltype"] == "Th1") & (df_clust["annot"] == 4), "annot_new"] = 3

    df_clust.loc[(df_clust["celltype"] == "Th2") & (df_clust["annot"] == 5), "annot_new"] = 4
    df_clust.loc[(df_clust["celltype"] == "Th2") & (df_clust["annot"] == 4), "annot_new"] = 5

    # reassign annot new to annot
    df_clust["annot"] = df_clust["annot_new"]

if int(nclust)==6:
    #
    df_clust.loc[(df_clust["celltype"] == "Th0") & (df_clust["annot"] == 2), "annot_new"] = 5
    df_clust.loc[(df_clust["celltype"] == "Th0") & (df_clust["annot"] == 5), "annot_new"] = 2

    df_clust.loc[(df_clust["celltype"] == "Th1") & (df_clust["annot"] == 3), "annot_new"] = 5
    df_clust.loc[(df_clust["celltype"] == "Th1") & (df_clust["annot"] == 5), "annot_new"] = 3

    df_clust.loc[(df_clust["celltype"] == "Th1/2") & (df_clust["annot"] == 5), "annot_new"] = 6
    df_clust.loc[(df_clust["celltype"] == "Th1/2") & (df_clust["annot"] == 6), "annot_new"] = 5

    # reassign annot new to annot
    df_clust["annot"] = df_clust["annot_new"]


df = pd.merge(df_tidy, df_clust, how = "inner", on = ["gene", "celltype"])

# plot kinetic genes in clusters --> columns:clusters, hue: celltype
xlabel = "time (h)"
ylabel = "log2FC(day0)"
xticks = [0,30,60,90,120]
yticks = [-2,-1,0,1,2]

g = sns.relplot(data = df, x = "time", col = "annot", hue = "celltype",
                y = "value", kind = "line",
                palette = ["tab:grey", "tab:blue", "tab:red", "purple"],
                hue_order = ["Th0", "Th1", "Th2", "Th1/2"],
                aspect = 1.0,
                legend = False)
g.set(xlabel = xlabel, ylabel = ylabel, xlim = (0,120),
      xticks = xticks, yticks = yticks)
sns.despine(top = False, right = False)
plt.show()

#g.savefig("../../figures/fig2/kinetic_clusters/clustering_timecourse1_" + nclust + "clust.svg")
#g.savefig("../../figures/fig2/kinetic_clusters/clustering_timecourse1_" + nclust + "clust.pdf")

if int(nclust) == 3:
    # get mean and sem for data frame
    df2 = df.groupby(["annot", "time"])["value"].agg(["sem", "std", "mean"])
    df2 = df2.reset_index()

    g = sns.FacetGrid(df2, hue="annot",
                      palette=["0.7", "0.4", "k"],
                      ylim=[-2, 2], xlim=[0, 120],
                      despine=False,
                      height=2, aspect = 1.2
                      )

    g.map_dataframe(sns.lineplot, x="time", y="mean")
    g.map_dataframe(plt.errorbar, x="time", y="mean", yerr="sem", capsize=3)
    g.set(xticks=[0, 24, 48, 72, 96, 120], xlabel="time (h)", ylabel="log2FC(day0)",
          yticks=[-2, -1, 0, 1, 2])
    #g.savefig("../../figures/fig2/kinetic_clusters/clustering_timecourse2_" + nclust + "clust.svg")
    #g.savefig("../../figures/fig2/kinetic_clusters/clustering_timecourse2_" + nclust + "clust.pdf")


# plot kinetic genes in clusters --> style:clusters
#df_clust = df_clust.loc[df_clust["celltype"] != "Th0"]
df_clust[["gene2", "iso1", "iso2"]] = df["gene"].str.split(r"_|\.", expand = True)
df_clust.drop(columns = ["gene"], inplace = True)
df_clust.rename(columns = {"gene2" : "gene"}, inplace = True)

def get_gene_no(df):
    mygenes = df["gene"].str.split(r"_|\.", expand = True)
    mygenes = mygenes.iloc[:,0].drop_duplicates()
    return(len(mygenes))

test = df.groupby(["celltype", "annot"]).apply(get_gene_no).reset_index()
test.columns = ["celltype", "annot", "ngenes"]

g = sns.catplot(data = test, x = "annot", y = "ngenes", hue = "celltype", kind = "bar",
                palette=["tab:grey", "tab:blue", "tab:red", "purple"],
                hue_order=["Th0", "Th1", "Th2", "Th1/2"], legend = False,
                aspect = 1.2, height = 2)
g.set(xlabel = "cluster", ylabel = "genes in cluster")
sns.despine(right = False, top = False)
plt.show()
g.savefig("../../figures/fig2/kinetic_clusters/clustering_barplot_" + nclust +"clust.svg")
g.savefig("../../figures/fig2/kinetic_clusters/clustering_barplot_" + nclust +"clust.pdf")

