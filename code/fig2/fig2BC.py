import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style = "ticks", context = "paper")

df_tidy = pd.read_csv("../../data/peine_logfc_tidy.csv")

#df_clust = pd.read_csv("../../data/clustered_genes/kinetic_genes_20210714_"+nclust+ "clust.csv")
# load cluster data
df_clust = pd.read_csv("../../data/clustered_genes/kinetic_genes_clusters_renamed_20211214.csv")
#
df = pd.merge(df_tidy, df_clust, how = "inner", on = ["gene", "celltype"])

# plot kinetic genes in clusters --> columns:clusters, hue: celltype
xlabel = "time (h)"
ylabel = "log2FC(day0)"
xticks = [0,30,60,90,120]
yticks = [-2,-1,0,1,2]

# get mean and sem for data frame
df2 = df.groupby(["annot", "time"])["value"].agg(["sem", "std", "mean"])
df2 = df2.reset_index()

g = sns.FacetGrid(df2, hue="annot", palette=["0.7", "0.4", "k"], ylim=[-2, 2], xlim=[0, 120])
g.map_dataframe(sns.lineplot, x="time", y="mean")
g.map_dataframe(plt.errorbar, x="time", y="mean", yerr="sem", capsize=3)
g.set(xticks=xticks, xlabel=xlabel, ylabel=ylabel, yticks=yticks)
g.savefig("../../figures/fig2B.pdf")


# plot kinetic genes in clusters --> style:clusters
#df_clust = df_clust.loc[df_clust["celltype"] != "Th0"]
df_clust[["gene2", "iso1", "iso2"]] = df["gene"].str.split(r"_|\.", expand = True)
df_clust.drop(columns = ["gene"], inplace = True)
df_clust.rename(columns = {"gene2" : "gene"}, inplace = True)

def get_gene_no(df):
    """
    get number of genes for each cell type and cluster
    """
    mygenes = df["gene"].str.split(r"_|\.", expand = True)
    # remove isoforms
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
g.savefig("../../figures/fig2C.pdf")

