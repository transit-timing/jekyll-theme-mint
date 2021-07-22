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

def sigma_clip(orbit_number, t, y_err, period_guess, t0_guess, nsigma_clip = 5):
    niter_sigmaclip = 10

    outlier = np.full(t.shape[0], False, dtype=bool)
    idxs = []

    for i in range(niter_sigmaclip):
        n_outliers_before = np.sum(outlier)



        # learn this is a magical function - it makes exactly what we want for this design matrix
        X = np.vander(orbit_number, N=2, increasing=True)
        # OR:
        # X = np.vstack((np.ones_like(x), x)).T

        Cov = np.diag(y_err**2)
        Cinv = np.linalg.inv(Cov) # we need the inverse covariance matrix

        # get the parameter covariance matrix ...
        theta_Cov = np.linalg.inv(X.T @ Cinv @ X)

        # and the best parameters using the new Python matrix operator
        theta_best = theta_Cov @ (X.T @ Cinv @ t)
        #print ("a, b = {:.3f}, {:.3f}".format(*theta_best))

        # add MLE estimate
        X_ = np.vander(orbit_number, N=2, increasing=True)
        t_model = X_ @ theta_best

        t0, per = theta_best

        t0_unc = np.sqrt(np.diag(theta_Cov))[0]
        per_unc = np.sqrt(np.diag(theta_Cov))[1]

        t_res = t-t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        #print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier))
        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
 
    idxs = np.unique(idxs)
 

    return idxs, per, t0, per_unc, t0_unc
 


def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)

lit_period_uncs = []
tess_period_uncs = []
lit_t0_uncs = []
tess_t0_uncs = []

for folder in os.listdir('/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/'):
  try:
    planet_name = folder.replace('-mid-times', '')
    #print('Planet: ', planet_name)

    path = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{folder}/'
    data = pd.read_csv(path + f'{planet_name}.csv')
    
    period = data['Period'].iloc[0]
    t = data['Mid-point']
    uncertainty = data['Uncertainty']

    #####################################
    # literature data only
    #####################################
    mask = data['Source'] != 'Our work'
    t = t[mask]
    uncertainty = uncertainty[mask]
 
    t_arr = t.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
 
    orbit_number = fold_time(t, t0, period)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    sigma = np.mean(uncertainty)

    idxs, P, t0, per_unc, t0_unc = sigma_clip(orbit_number, t, uncertainty, period, t0)

    lit_period_uncs.append(per_unc)
    lit_t0_uncs.append(t0_unc)


    period = data['Period'].iloc[0]
    t = data['Mid-point']
    uncertainty = data['Uncertainty']
    #####################################
    # literature and TESS data only
    #####################################
    t_arr = t.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
 
    orbit_number = fold_time(t, t0, period)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    sigma = np.mean(uncertainty)

    idxs, P, t0, per_unc, t0_unc = sigma_clip(orbit_number, t, uncertainty, period, t0)

    tess_period_uncs.append(per_unc)
    tess_t0_uncs.append(t0_unc)
  except Exception as e: print(e)
 
lit_period_uncs = np.array(lit_period_uncs) * 24 * 60 * 60  
tess_period_uncs = np.array(tess_period_uncs) * 24 * 60 * 60  
lit_t0_uncs = np.array(lit_t0_uncs) * 24 * 60 * 60  
tess_t0_uncs = np.array(tess_t0_uncs) * 24 * 60 * 60  


arr = np.linspace(min(lit_period_uncs), max(lit_period_uncs), 100)

plt.plot(lit_period_uncs, tess_period_uncs, '.k')
plt.plot(arr, arr)
plt.xlabel(r'$\sigma_{P_{lit}}$ [sec]')
plt.ylabel(r'$\sigma_{P_{all data}}$ [sec]')
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/3_period/0_period_uncs_vs_lit_period_uncs.png')
plt.show()


arr = np.linspace(min(lit_t0_uncs), max(lit_t0_uncs), 100)
plt.plot(lit_t0_uncs, tess_t0_uncs, '.k')
plt.plot(arr, arr)
plt.xlabel(r'$\sigma_{T_{0_{lit}}}$ [sec]')
plt.ylabel(r'$\sigma_{T_{0_{all data}}}$ [sec]')
plt.xscale('log')
plt.yscale('log')
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/3_period/0_t0_uncs_vs_lit_t0_uncs.png')
plt.show()

'''
    #fig, axs = plt.subplots(1)
    fig = plt.figure()

    plt.errorbar(orbit_number[m], o_c[m]*24*60, yerr = uncertainty[m]*24*60, fmt='o', mew = 1, mfc = '#069AF3', alpha = 1,  markersize = 9, ecolor = '#069AF3')
    plt.errorbar(orbit_number[~m], o_c[~m]*24*60, yerr = uncertainty[~m]*24*60, fmt='o', mew = 1, mfc = 'white', alpha = 0.8,  markersize = 9, ecolor = '#808080')
 
    plt.xlabel('Epoch')
    plt.ylabel(r'$\Delta t$ [min]')
    plt.title(f'{planet_name.upper()}')
    plt.xlim(t_min - 50, t_max + 50)
    #legend = f't0 = %.4f Period = %.4f d' % (t0, period)
    #axs[0].legend([f'{legend}'], handletextpad=-2.0, handlelength=0)

 

    plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/o_c_figures/o_c_linear_regression/{planet_name}_o_c.png')
    #plt.show()
    plt.close(fig)
'''

