import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_scatterplot(df, model):

    models_available = df["model"].drop_duplicates().values
    assert model in models_available

    # only focus on this model
    df = df.loc[df["model"] == model, :]

    fig, ax = plt.subplots(figsize=(6, 6))
    ptsize = 20
    sns.scatterplot(data=df, x="beta_th1", y="beta_th2", s=ptsize, hue="category")

    # add symmetric lines
    lw = 1
    ax.axhline(0, color = "k", lw  = lw, ls = "--")
    ax.axvline(0, color = "k", lw  = lw, ls = "--")

    ticks = [-1,0,1]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_title(model)

    plt.show()

    if model == "Th1/2":
        model = "Th12"
    fig.savefig("../../figures/FigS4_A_" + model + ".pdf")

# plotting params

sns.set(style = "ticks", context = "talk")

# read data
baseline = "degs_Th1_Th2_cor0.3"

df = pd.read_csv("../../data/linear_model/linear_model_processed_"+ baseline + ".csv")

def norm_coeff(df):
    df["beta_th1"] = df["beta_th1"] / df["beta_th1"].max()
    df["beta_th2"] = df["beta_th2"] / df["beta_th2"].max()
    return df
  
# normalize coefficients?
df = norm_coeff(df)

plot_scatterplot(df, "Th1/2")
plot_scatterplot(df, "Th0")
