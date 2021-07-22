import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

# all planets (circle) vs analyzed planets (triangle) and all planets colored by Vmag

df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_tep.csv')
targets = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/target_list.csv')

sns.set_style("dark") # darkgrid, whitegrid, dark, white, and ticks.
g = sns.relplot(x="Period(day)", y="R_b_over_R_A_squared", hue="Vmag", palette="bwr", data=df, edgecolor="black", markers={
  0: "o",
  1: "v",}, style = 'Analyzed', linewidth=0.8) #ch:r=-.5,l=.75
#sns.relplot(x="Period(day)", y="R_b_over_R_A_squared", hue="Vmag", palette="bwr", data=targets, edgecolor="black", marker='v', linewidth=0.8) #ch:r=-.5,l=.75
g.fig.set_figwidth(16.27)
g.fig.set_figheight(11.7)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Orbital period (days)', fontsize=18)
plt.ylabel(r'$(\frac{R_{p}}{R_{*}})^2$', fontsize=18)
plt.xticks(fontsize= 18)
plt.yticks(fontsize= 18)
plt.savefig('/Users/kate/Documents/research/paper/2_figures/0_depth_vs_period/depth_vs_period_1.png')
plt.show()