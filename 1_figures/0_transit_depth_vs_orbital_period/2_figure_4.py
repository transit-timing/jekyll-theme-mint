import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

# not color by magnitude anything  
# left subplot: all planets (blue), analyzed planets (red)
# right subplot: analyzed planets: blue (spoc), white (no tess data), red (qlp)

df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_tep.csv')
targets = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_target_list.csv')
targets = targets.fillna(value='No TESS data') 

fig, axes = plt.subplots(1, 2, figsize=(22, 10))
#plt.subplot(1, 2, 1)
#sns.set_style("dark") # darkgrid, whitegrid, dark, white, and ticks.
g1 = sns.scatterplot(ax=axes[0],x="Period(day)", y="R_b_over_R_A_squared", hue="Analyzed", palette="bwr", data=df, edgecolor="black", linewidth=0.9) #ch:r=-.5,l=.75

axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlabel(r'Orbital period (days)', fontsize=18)
axes[0].set_ylabel('Transit depth', fontsize=18)
#axes[0].set_xticks(fontsize= 18)
#axes[0].set_yticks(fontsize= 18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.subplot(1, 2, 2)
g2 = sns.scatterplot(ax=axes[1], x="Period(day)", y="R_b_over_R_A_squared", hue="TESS source", palette='bwr', data=targets, edgecolor="black", linewidth=0.9) #ch:r=-.5,l=.75

axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel('Orbital period (days)', fontsize=18)
axes[1].set_ylabel('Transit depth', fontsize=18)
#axes[1].set_xticks(fontsize= 18)
#axes[1].set_yticks(fontsize= 18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
axes[0].tick_params(axis='both', which='major', labelsize=18)
axes[1].tick_params(axis='both', which='major', labelsize=18)
#g.fig.set_figwidth(16.27)
#g.fig.set_figheight(11.7)
plt.savefig('/Users/kate/Documents/research/paper/2_figures/0_depth_vs_period/depth_vs_period_4.png')
plt.show()