import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy.stats import chi2

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
 
# read in CSV files that contain chi-square and ndof for literature + TESS data and TESS data only
tess_chi_square =pd.read_csv('/Users/kate/Documents/research/paper/3_tables/tess_chi_square_vs_ndof.csv')
all_chi_square = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/all_chi_square_vs_ndof.csv')


ndofs = all_chi_square['ndof']
chi_squares = all_chi_square['Chi-square']
ndofs_tess = tess_chi_square['ndof']
chi_squares_tess = tess_chi_square['TESS chi-square'] 


# mask targets for which only two data points are available (all data)
mask_2dof = ndofs > 0
ndofs = ndofs[mask_2dof]
chi_squares = chi_squares[mask_2dof]

# mask targets for which only two data points are available (TESS data only)
mask_2dof = ndofs_tess > 0
ndofs_tess = ndofs_tess[mask_2dof]
chi_squares_tess = chi_squares_tess[mask_2dof]

print('Number of points used to create left subplot: ', chi_squares.shape[0])
print('Number of points used to create right subplot: ',chi_squares_tess.shape[0] )

###################################################################################################
# create dataframes with chi-squares after masking out targets with ndofs < 1

df = pd.DataFrame({'ndofs': ndofs, 'chi_squares': chi_squares})
df_tess = pd.DataFrame({'ndofs_tess': ndofs_tess, 'chi_squares_tess': chi_squares_tess})
 
# for plotting identity line on the left subplot
n_min = np.min(ndofs)
n_max = np.max(ndofs)
n_arr = np.linspace(n_min - 5, n_max + 50, 100)

# for plotting identity line on the left subplot
n_min_tess = np.min(ndofs_tess)
n_max_tess = np.max(ndofs_tess)
n_arr_tess = np.linspace(n_min_tess - 5, n_max_tess + 50, 100)

# boundary of the region encapsulating 2.5% to 97.5% of the chi-square distribution 
ndof_range_tess = np.arange(n_min_tess, n_max_tess+50)
y1_tess = chi2.ppf(0.025, ndof_range_tess) 
y2_tess = chi2.ppf(0.975, ndof_range_tess) 

ndof_range = np.arange(n_min, n_max+50)
y1 = chi2.ppf(0.025, ndof_range) 
y2 = chi2.ppf(0.975, ndof_range) 

###################################################################################################
# make left subplot
fig, axes = plt.subplots(1, 2, figsize=(22, 10))

g1 = sns.scatterplot(ax=axes[0],x='ndofs', y='chi_squares', data = df, color = 'blue', alpha = 0.9, edgecolor="black", linewidth=0.9)  
ndofs_arr = np.linspace( ndofs.min(), ndofs.max()+50, 300)
ndofs_tess_arr = np.linspace( ndofs_tess.min(), ndofs_tess.max()+50, 300)

axes[0].plot(n_arr, n_arr, '-', color = 'gray')
axes[0].set_ylabel(r'$\chi^2$',  fontsize=18)
axes[0].set_xlabel(r'$N$',  fontsize=18)
axes[0].fill_between(ndof_range, y1, y2,  color = 'gray', alpha = 0.2)
axes[0].set_xscale('log')
axes[0].set_yscale('log')

# make right subplot
g2 = sns.scatterplot(ax=axes[1],x='ndofs_tess', y='chi_squares_tess', data = df_tess, color = 'blue', alpha = 0.9, edgecolor="black", linewidth=0.9)

axes[1].plot(n_arr_tess, n_arr_tess, '-', color = 'gray')
axes[1].set_ylabel(r'$\chi^2$',  fontsize=18)
axes[1].set_xlabel(r'$N$',  fontsize=18)
axes[1].fill_between(ndof_range_tess, y1_tess, y2_tess,  color = 'gray', alpha = 0.2)
axes[1].set_xscale('log')
axes[1].set_yscale('log')

plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/1_chi_square/chi_square_vs_ndof.png')
plt.show()
plt.close(fig)


 