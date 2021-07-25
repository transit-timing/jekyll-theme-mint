import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import batman
import emcee
import corner
import os, sys, time
from scipy.optimize import minimize
import seaborn as sns

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE) 


df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/P_T0_errors.csv')

lit_period_uncs = df['P error (literature data)'] * 24 * 60 * 60  
all_data_period_uncs = df['P error (all data)'] * 24 * 60 * 60  
lit_t0_uncs = df['T0 error (literature data)'] * 24 * 60 * 60  
all_data_t0_uncs = df['T0 error (all data)'] * 24 * 60 * 60  
targets = df['Target']

 
mask = lit_period_uncs <  all_data_period_uncs
print('targets for which precision in period did not improve from including TESS data: ', targets[mask])
print(f'Improved period precision for {targets[~mask].shape[0]} targets')
mask = lit_t0_uncs <  all_data_t0_uncs
print('targets for which precision in T0 did not improve from including TESS data: ', targets[mask])
print(f'Improved T0 precision for {targets[~mask].shape[0]} targets')
arr = np.linspace(min(lit_period_uncs), max(lit_period_uncs), 100)

fig = plt.figure()
plt.fill_between(arr, 0.5 * arr, 1.5 * arr,  color = 'gray', alpha = 0.2)
plt.plot(lit_period_uncs, all_data_period_uncs, '.k')
plt.plot(arr, arr)
plt.xlabel(r'$\sigma_{P_{lit}}$ [sec]')
plt.ylabel(r'$\sigma_{P_{all}}$ [sec]')
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/3_period/0_period_uncs_vs_lit_period_uncs.png')
plt.show()
plt.close(fig)

fig = plt.figure()
arr = np.linspace(min(lit_t0_uncs), max(lit_t0_uncs), 100)
plt.fill_between(arr, 0.5 * arr, 1.5 * arr,  color = 'gray', alpha = 0.2)
plt.plot(lit_t0_uncs, all_data_t0_uncs, '.k')
plt.plot(arr, arr)
plt.xlabel(r'$\sigma_{T_{0_{lit}}}$ [sec]')
plt.ylabel(r'$\sigma_{T_{0_{all}}}$ [sec]')
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/3_period/0_t0_uncs_vs_lit_t0_uncs.png')
plt.show()

 
