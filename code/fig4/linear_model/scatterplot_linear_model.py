import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_scatterplot(df, model, baseline):

    models_available = df["model"].drop_duplicates().values
    assert model in models_available

    # only focus on this model
    df = df.loc[df["model"] == model, :]

    if model == "Th0":
        sname = "Th0"
    else:
        sname = "ThMix"

    fig, ax = plt.subplots(figsize=(4.5, 4))
    palette = ["tab:blue", "tab:red", "k", "grey"]

    palette = ["#D95F02", "#1B9E77", "#E6AB02", "#666666"]
    hue_order = ["Ind", "Superpos", "Th1like", "Th2like"]

    ptsize = 20

    sns.scatterplot(data=df, x="beta_th1", y="beta_th2", s=ptsize, hue="category", alpha=1,
                    palette = palette, hue_order = hue_order, legend=False)

    # add symmetric lines
    lw = 1
    #ax.axline((0, 0), slope=np.tan(np.pi / 8), ls="-", c="k", lw=lw)
    #ax.axline((0, 0), slope=np.tan(3 * np.pi / 8), ls="-", c="k", lw=lw)
    #ax.axline((0, 0), slope=-np.tan(np.pi / 8), ls="-", c="k", lw=lw)
    #ax.axline((0, 0), slope=-np.tan(3 * np.pi / 8), ls="-", c="k", lw=lw)
    ax.axhline(0, color = "k", lw  = lw, alpha = 1, ls = "--")
    ax.axvline(0, color = "k", lw  = lw, alpha = 1, ls = "--")

    xticks = [-1, 0, 1]
    yticks = [-1, 0, 1]

    #ax.set_xlim([-1,1])
    #ax.set_ylim([-1,1])

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    xlabel = r"$\beta_\mathrm{Th1}$"
    ylabel = r"$\beta_\mathrm{Th2}$"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.show()
    fig.savefig(sdir + "linear_model_coeff_" + sname + baseline + ".svg")
    fig.savefig(sdir + "linear_model_coeff_" + sname + baseline +".pdf")

# plotting params
sdir = "../../../figures/fig4/linear_model/"
sns.set(style = "ticks", context = "talk")
palette = ["tab:grey", "tab:blue", "tab:red", "tab:purple"]

# read data
baseline = "degs_Th1_Th2_cor0.3"
baseline2 = "kinetic_th12"
df = pd.read_csv("../../../data/linear_model/linear_model_processed_"+ baseline + ".csv")

def norm_coeff(df):
    df["beta_th1"] = df["beta_th1"] / df["beta_th1"].max()
    df["beta_th2"] = df["beta_th2"] / df["beta_th2"].max()
    return df
# normalize coefficients?
df = norm_coeff(df)
#df2 = pd.read_csv("../../../data/linear_model/linear_model_processed_"+ baseline2 + ".csv")

plot_scatterplot(df, "Th1/2", baseline)
plot_scatterplot(df, "Th0", baseline)
#plot_scatterplot(df2, "Th1/2", baseline2)
#plot_scatterplot(df2, "Th0", baseline2)

# add angle information for categories
df.loc[:,["angle"]] = np.arctan2(df.loc[:,"beta_th2"], df.loc[:,"beta_th1"])
angle_superpos1 = (np.pi/8)
angle_superpos2 = (3*np.pi/8)

# assign new clusters
df["clust_new"] = "Superpos"
df.loc[(df["angle"] > (np.pi/8)) & (df["angle"] < (3*np.pi/8)), "clust_new"] = "Superpos"
df.loc[(df["angle"] > (-7*np.pi/8)) & (df["angle"] < (-5*np.pi/8)), "clust_new"] = "Superpos"
df.loc[(df["angle"] < (-3*np.pi/8)) & (df["angle"] > (-5*np.pi/8)), "clust_new"] = "Th2like"
df.loc[(df["angle"] > (3*np.pi/8)) & (df["angle"] < (5*np.pi/8)), "clust_new"] = "Th2like"
df.loc[(df["angle"] < (np.pi/8)) & (df["angle"] > (-1*np.pi/8)), "clust_new"] = "Th1like"
df.loc[(df["angle"] > (7*np.pi/8)) | (df["angle"] < (-7*np.pi/8)), "clust_new"] = "Th1like"


df.to_csv("../../../data/linear_model/linear_model_processed_"+ baseline + ".csv", index = False)