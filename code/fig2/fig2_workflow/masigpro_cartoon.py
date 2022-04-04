#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 11:40:55 2021

@author: burt
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

rc_dict = {
    "axes.linewidth" : 0.8,
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    "xtick.major.size" : 3,
    "ytick.major.size" : 3,
    "font.size" : 8,
    "xtick.labelsize" : 6,
    "ytick.labelsize" : 6,
    "lines.markersize" : 4,
    "lines.linewidth" : 1.0,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
}
sns.set_context("paper", rc = rc_dict)
sns.set_style("ticks")

x = np.linspace(0,10,10)

y = np.array([0.6,0.8,0.7,0.65,0.8,0.75,0.35,0.25, 0.2,0.1])
y2 = np.array([0, 0.1, 0,0.05,0.07,0.08,0.1,0.3,0.7,0.8])

# newly added y4
y4 = np.array([0.3,0.4,0.25,0.3,0.6,0.8, 0.9, 0.3,0.1,0.15])

degree = 4
w1 = np.polyfit(x, y, degree)
m1 = np.poly1d(w1)
w2 = np.polyfit(x, y2, degree)
m2 = np.poly1d(w2)

mycol1 = "0.3"
y3 = np.random.normal(loc = 0.3,scale = 0.05, size = 10)

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize = (5,1.5))
ax1.scatter(x,y, color = mycol1)
#ax1.scatter(x,y2)
ax1.scatter(x,y3, color = "white", edgecolor = mycol1, zorder = 100)
ax1.plot(x,y, color = mycol1, label = "kinetic")
#ax1.plot(x,y2)
ax1.plot(x,y3, color = mycol1, label = "not kinetic", ls = "dashed")
ax1.legend()


spl1 = UnivariateSpline(x, y)

ax2.scatter(x,y, color = mycol1)
#ax2.scatter(x,y2)

x2 = np.linspace(0,10,50)

# increase res for splines fit (looks smoother, works also with normal x)
ax2.plot(x2,m1(x2), color = mycol1)
#ax2.plot(x2,m2(x2))

#fig.savefig("../../figures/fig2/masigpro_cartoon.pdf")

res = 50
for i in range(res):
    lw = 0.5
    sd = 0.03
    y = y + np.random.normal(0,sd,10)
    weights = np.polyfit(x, y, degree)
    model = np.poly1d(weights)
    ax3.plot(x2,model(x2), color = "k", lw = lw, zorder = 1000)
    
    y2 = y2 + np.random.normal(0,sd,10)
    weights = np.polyfit(x, y2, degree)
    model = np.poly1d(weights)
    ax3.plot(x2,model(x2), lw = lw, color = "0.5")

    y4 = y4 + np.random.normal(0,sd,10)
    weights = np.polyfit(x, y4, degree)
    model = np.poly1d(weights)
    ax3.plot(x2,model(x2), color = "0.9", lw = lw, zorder = 500)

for ax in (ax1,ax2,ax3):
    ax.set_yticks([0,0.5,1])
    ax.set_xticks([0,5,10])
    ax.set_xlim([0,10])
    ax.set_xlabel("time (h)")
ax1.set_ylabel("expression (a.u.)")
ax2.set_ylabel("")
ax3.set_ylabel("")


ax1.set_title("identify kinetic genes")
#ax2.set_title(r"$y=\beta_0 + \beta_1 t + \beta_2 t^2 +$ ...")
ax2.set_title("regression fit kinetic genes")
ax3.set_title("identify kinetic patterns")

plt.tight_layout()
plt.show()

fig.savefig("../../../figures/fig2/mspro_cartoon.pdf")
fig.savefig("../../../figures/fig2/mspro_cartoon.svg")
