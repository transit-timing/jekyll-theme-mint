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
from scipy.stats import chi2


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

chi_squares = []
ndofs = []

chi_squares_tess = []
ndofs_tess = []


def sigma_clip(orbit_number, t, std, nsigma_clip = 5):
    niter_sigmaclip = 10

    outlier = np.full(t.shape[0], False, dtype=bool)
    idxs = []

    for i in range(niter_sigmaclip):
        n_outliers_before = np.sum(outlier)
        k, t0 = np.polyfit(orbit_number, t, 1, w=1./std**2)
        t_model = k*orbit_number + t0
        t_res = t-t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        #print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier))
        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
    #print('idxs ', idxs)
    idxs = np.unique(idxs)
    #idxs = np.unique(idxs)
   #idxs = idxs.flatten()
    print(idxs)
    return idxs, k, t0
 

def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)


for folder in os.listdir('/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/'):
  try:
    planet_name = folder.replace('-mid-times', '')
 

    path = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{folder}/'
    data = pd.read_csv(path + f'{planet_name}.csv')

    period = data['Period'].iloc[0]
    t = data['Mid-point']
    uncertainty = data['Uncertainty']
 
    t_arr = t.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
 
    orbit_number = fold_time(t, t0, period)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    sigma = np.mean(uncertainty)

    idxs, P, t0 = sigma_clip(orbit_number, t, uncertainty)



    m = np.ones(t.shape, bool) 
    m[idxs[:]] = False

    # select all mid-points that were not clipped. 
    o_c = t[m] - (t0 + P*orbit_number[m])

    chi_square = np.sum((o_c/uncertainty[m])**2)
    ndof = orbit_number[m].shape[0] - 2

    chi_squares.append(chi_square)
    ndofs.append(ndof)
    #print('Planet: ', planet_name, ' ndof: ', ndof)

    # tess data only
    mask = data['Source'] == 'Our work'  
    orbit_number = orbit_number[mask]
    t = t[mask]
    uncertainty = uncertainty[mask]

    idxs, P, t0 = sigma_clip(orbit_number, t, uncertainty)
    m = np.ones(t.shape, bool) 
    m[idxs[:]] = False

    if np.sum(m) > 0:
        o_c = t[m] - (t0 + P*orbit_number[m])

        chi_square = np.sum((o_c/uncertainty[m])**2)
        ndof = orbit_number[m].shape[0] - 2


        chi_squares_tess.append(chi_square)
        ndofs_tess.append(ndof)

  except:
    pass

ndofs= np.array(ndofs)
chi_squares = np.array(chi_squares)

# mask targets for which only two data points are available 
mask_2dof = ndofs > 0
ndofs = ndofs[mask_2dof]
chi_squares = chi_squares[mask_2dof]

ndofs_tess = np.array(ndofs_tess)
chi_squares_tess = np.array(chi_squares_tess)

mask_2dof = ndofs_tess > 0
ndofs_tess = ndofs_tess[mask_2dof]
chi_squares_tess = chi_squares_tess[mask_2dof]


#m = chi_squares < 4000
#ndofs= ndofs[m]
#chi_squares = chi_squares[m]


df = pd.DataFrame({'ndofs': ndofs, 'chi_squares': chi_squares})
#df.to_csv('/Users/kate/Desktop/check.csv')
df_tess = pd.DataFrame({'ndofs_tess': ndofs_tess, 'chi_squares_tess': chi_squares_tess})
print(df) 
n_min = np.min(ndofs)
n_max = np.max(ndofs)
n_arr = np.linspace(n_min - 5, n_max + 50, 100)


n_min_tess = np.min(ndofs_tess)
n_max_tess = np.max(ndofs_tess)
n_arr_tess = np.linspace(n_min_tess - 5, n_max_tess + 50, 100)

chi_max = ndofs + 3*np.sqrt(2*ndofs)
chi_min = ndofs - 3*np.sqrt(2*ndofs)
 

ndof_range_tess = np.arange(n_min_tess, n_max_tess+50)
y1_tess = chi2.ppf(0.025, ndof_range_tess) 
y2_tess = chi2.ppf(0.975, ndof_range_tess) 

ndof_range = np.arange(n_min, n_max+50)
y1 = chi2.ppf(0.025, ndof_range) 
y2 = chi2.ppf(0.975, ndof_range) 

fig, axes = plt.subplots(1, 2, figsize=(22, 10))

g1 = sns.scatterplot(ax=axes[0],x='ndofs', y='chi_squares', data = df, color = 'blue', alpha = 0.9, edgecolor="black", linewidth=0.9) #ch:r=-.5,l=.75
ndofs_arr = np.linspace( ndofs.min(), ndofs.max()+50, 300)
ndofs_tess_arr = np.linspace( ndofs_tess.min(), ndofs_tess.max()+50, 300)
#terms = np.linspace( ndofs_arr - 3*np.sqrt(2*ndofs_arr),  ndofs + 3*np.sqrt(2*ndofs_arr, 100))

#y1 = ndofs_arr - 3*np.sqrt(2*ndofs_arr)
#y2 = ndofs_arr + 3*np.sqrt(2*ndofs_arr)
#y1_tess = ndofs_tess_arr - 3*np.sqrt(2*ndofs_tess_arr)
#y2_tess =  ndofs_tess_arr + 3*np.sqrt(2*ndofs_tess_arr)
 
axes[0].plot(n_arr, n_arr, '-', color = 'gray')
#axes[0].plot(ndofs, y1, '.', color = 'gray')
#axes[0].plot(ndofs, y2, '.', color = 'gray')
#for term in terms:
#    axes[0].plot(ndofs_arr, term, '.', color = 'gray', alpha = 0.05)
 
#slope_arr = np.linspace(-3*np.sqrt(ndofs), 3*sqrt(ndofs), 100)
#for i in range(100):


axes[0].set_ylabel(r'$\chi^2$',  fontsize=18)
axes[0].set_xlabel(r'$N$',  fontsize=18)
#axes[0].set_xlim(n_min, ndofs.max()+50)
axes[0].fill_between(ndof_range, y1, y2,  color = 'gray', alpha = 0.2)
axes[0].set_xscale('log')
axes[0].set_yscale('log')

g2= sns.scatterplot(ax=axes[1],x='ndofs_tess', y='chi_squares_tess', data = df_tess, color = 'blue', alpha = 0.9, edgecolor="black", linewidth=0.9)


axes[1].plot(n_arr_tess, n_arr_tess, '-', color = 'gray')
axes[1].plot(n_arr_tess, n_arr_tess, '-', color = 'gray')
#axes[1].plot(ndofs_tess, chi_max_tess, '-', color = 'gray')
#axes[1].plot(ndofs_tess, chi_min_tess, '-', color = 'gray')
#for term in terms_tess:
#    axes[1].plot(ndofs_tess, term, '-', color = 'gray', alpha = 0.05)
 

axes[1].set_ylabel(r'$\chi^2$',  fontsize=18)
axes[1].set_xlabel(r'$N$',  fontsize=18)
#axes[1].set_xlim(n_min_tess - 1, n_max_tess + 1)
axes[1].fill_between(ndof_range_tess, y1_tess, y2_tess,  color = 'gray', alpha = 0.2)
axes[1].set_xscale('log')
axes[1].set_yscale('log')

plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/1_chi_square/chi_square_vs_ndof.png')
plt.show()
plt.close(fig)


 