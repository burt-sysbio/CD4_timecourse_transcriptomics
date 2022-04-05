import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style = "ticks", context = "talk")


baseline = "degs_Th1_Th2_cor0.3"
mydir = "../../data/linear_model/linear_model_processed_" + baseline + ".csv"

df = pd.read_csv(mydir)


df = df.groupby(["model"])["category"].value_counts(normalize = True).mul(100).rename('percent').reset_index()


# plot individually hclust and linear model
sns.set(context = "poster", style = "ticks")

df_bar = df.loc[df["model"] == "Th1/2"]

g = sns.catplot(data = df_bar, x = "category", y = "percent", kind = "bar",
                palette = ["k"])

g.set(xlabel = "", ylabel = "genes (%)")
g.set_xticklabels(rotation=90)

plt.show()

g.savefig("../../figures/Fig4B.pdf")
