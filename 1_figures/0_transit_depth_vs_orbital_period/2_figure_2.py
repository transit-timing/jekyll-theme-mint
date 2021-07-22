import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

#  
# 1 plot: analyzed planets: blue (spoc), white (no tess data), red (qlp)

df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_target_list.csv')
#print(df["TESS source"] == 'NaN')
df = df.fillna(value='No TESS data') 
 
#sns.set_style("dark") # darkgrid, whitegrid, dark, white, and ticks.
g = sns.relplot(x="Period(day)", y="R_b_over_R_A_squared", hue="TESS source", palette='coolwarm', data=df, edgecolor="black", linewidth=0.8) #ch:r=-.5,l=.75
#sns.relplot(x="Period(day)", y="R_b_over_R_A_squared", hue="Vmag", palette="bwr", data=targets, edgecolor="black", marker='v', linewidth=0.8) #ch:r=-.5,l=.75
g.fig.set_figwidth(16.27)
g.fig.set_figheight(9.7)
plt.xscale('log')
plt.yscale('log')
#g.set_xticks([0.001, 0.005, 0.01, 1])
plt.xlabel('Orbital period (days)', fontsize=24)
plt.ylabel(r'$(\frac{R_{p}}{R_{*}})^2$', fontsize=18)
plt.xticks(fontsize= 18)
plt.yticks(fontsize= 18)
plt.savefig('/Users/kate/Documents/research/paper/2_figures/0_depth_vs_period/depth_vs_period_2.png')
plt.show()

#
#{
#  'No TESS data': "gray",
#  'SPOC': "blue", 'QLP': "red"}
#
#