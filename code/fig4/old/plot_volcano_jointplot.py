import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("../../data/minimal_correlation.csv")

# watch out for y axis mincor has to be adjusted for Th1 and Th2...
g = sns.JointGrid()
sns.scatterplot(data = df, x = "rsq", y = "mincor_idx", hue = "sig",
                palette = ["grey", "red"],
                s=8, alpha =0.5, ax = g.ax_joint)
sns.kdeplot(data = df, x = "rsq", ax = g.ax_marg_x, color = "grey")
sns.kdeplot(data = df, y = "cor_idx", ax = g.ax_marg_y, color = "grey")

g.ax_joint.set_xlim([0.2,1.01])
g.ax_joint.set_ylim([-0.05,1.5])

g.ax_joint.axhline(0.3, ls = "--", c = "k")
g.ax_joint.axvline(0.6, ls = "--", c= "k")

plt.show()