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

 
def sigma_clip(orbit_number, t, y_err, nsigma_clip = 5):
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

        t0, p = theta_best


        t_model = p*orbit_number + t0
        t_res = t-t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        t0_unc = np.sqrt(np.diag(theta_Cov))[0]
        per_unc = np.sqrt(np.diag(theta_Cov))[1]

        #print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier))
        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
    idxs = np.unique(idxs)
 
    return idxs, p, t0 #, t0_unc, per_unc

 
def sigma_clip_quadratic(orbit_number, t, y_err, nsigma_clip = 5):
    niter_sigmaclip = 10

    outlier = np.full(t.shape[0], False, dtype=bool)
    idxs = []

    for i in range(niter_sigmaclip):
        n_outliers_before = np.sum(outlier)



        # learn this is a magical function - it makes exactly what we want for this design matrix
        X = np.vander(orbit_number, N=3, increasing=True)
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
        X_ = np.vander(orbit_number, N=3, increasing=True)
        t_model = X_ @ theta_best

        t0, p, dpde_over_2 = theta_best


        t_model = t0 + p * orbit_number + dpde_over_2 * orbit_number**2
        t_res = t - t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        t0_unc = np.sqrt(np.diag(theta_Cov))[0]
        per_unc = np.sqrt(np.diag(theta_Cov))[1]
        dpde_over_2_unc = np.sqrt(np.diag(theta_Cov))[2]


        #print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier))
        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
    idxs = np.unique(idxs)
 
    return idxs, p, dpde_over_2, t0 #, t0_unc, per_unc, dpde_over_2_unc


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
    t_err = data['Uncertainty']
    source = data['Source']
 

    ###########################################
    # if uncertainty < 26 sec, set it to 30 sec
    small_uncertainty_mask = t_err < 0.0003
    t_err[small_uncertainty_mask] = 0.0003
    ###########################################

    t_arr = t.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
 
    orbit_number = fold_time(t, t0, period)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    orbit_number_arr = np.linspace(min(orbit_number), max(orbit_number), 100)


    idxs, P_init, t0_init = sigma_clip(orbit_number, t, t_err)
    m = np.ones(t.shape, bool) 
    m[idxs[:]] = False

 
    if np.sum(source == 'Our work') == 0:
        # after sigma clipping with the linear model
        orbit_number_after_sigma_clip = orbit_number[m]
        t_err_after_sigma_clip = t_err[m]
        t_after_sigma_clip = t[m]

        # if chi square > 1, multiply uncertainties by sqrt(chi square/ndof).
        # when calculating chi square, we exclude 5-sigma outliers

        model = t0_init + P_init * orbit_number_after_sigma_clip
        chi_square = np.sum((t_after_sigma_clip - model)**2/t_err_after_sigma_clip**2)
        ndof = t_after_sigma_clip.shape[0] - 2

        if chi_square/ndof > 1: factor = np.sqrt(chi_square/ndof)
        else: factor = 1

        # enlarge uncertainty
        t_err_after_sigma_clip = t_err_after_sigma_clip * factor

        # re-fit linear model after enlarging uncertainties
        idxs, P, t0 = sigma_clip(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)
        linear_model = t0 + orbit_number_arr*P
        o_c = t - (t0 + P*orbit_number)
 
        idxs, P, dpde_over_2, t0 = sigma_clip_quadratic(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)

        # plot o-c diagram     
        fig = plt.figure()

        plt.errorbar(orbit_number[m], o_c[m]*24*60, yerr = t_err[m]*24*60, fmt='o', mew = 1, mfc = '#069AF3', alpha = 1,  markersize = 9, ecolor = '#069AF3')
        plt.plot(orbit_number_arr, ((t0 + P*orbit_number_arr + dpde_over_2 * orbit_number_arr**2) - linear_model)*24*60, 'black') 
        plt.errorbar(orbit_number[~m], o_c[~m]*24*60, yerr = t_err[~m]*24*60, fmt='o', mew = 1, mfc = 'white', alpha = 0.8,  markersize = 9, ecolor = '#808080')
     
        plt.xlabel('Epoch')
        plt.ylabel(r'$\Delta t$ [min]')
        plt.title(f'{planet_name.upper()}')
        plt.xlim(t_min - 50, t_max + 50)
     

        plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/o_c_figures/o_c_quadratic/{planet_name}_o_c.png')
        #plt.show()
        plt.close(fig)

    else:
        ##############################################################################################
        # tess transit times
        src =  source == 'Our work'
     
        mask = src.to_numpy() & m
        t_tess = t.to_numpy()[mask]
        uncertainty_tess = t_err.to_numpy()[mask]
        orbit_number_tess = orbit_number[mask]
     
        model_t_tess = t0_init + orbit_number_tess * P_init

        chi_tess = np.sum((t_tess - model_t_tess)**2/uncertainty_tess**2)

        ##############################################################################################
        # literature transit times
        nonsrc =  source != 'Our work'
        mask_lit = nonsrc.to_numpy() & m
        t_lit = t.to_numpy()[mask_lit]
        uncertainty_lit = t_err.to_numpy()[mask_lit]
        orbit_number_lit = orbit_number[mask_lit]
        model_t_lit = t0_init + orbit_number_lit * P_init


        # chi square of literature values
        chi_lit = np.sum((t_lit - model_t_lit)**2/uncertainty_lit**2)  
        ##############################################################################################

        ndof = t[m].shape[0] - 2 # number of degrees of freedom = number of data points after sigma clipping - 2 free parameters (linear model)
   
        # calculate factor by which to enlarge literature uncertainties
        if ndof > 0:
          if t_tess.shape[0]>=1:
            if chi_lit/(ndof - chi_tess) > 1: factor = np.sqrt(chi_lit/(ndof - chi_tess))
            else: factor = 1
          else:
            if chi_lit/ndof > 1: factor = np.sqrt(chi_lit/ndof)
            else: factor = 1    
        else:
          factor = 1     

        # enlarge literature uncertainties
        t_err[nonsrc] = t_err[nonsrc] * factor

        # after sigma clipping with the linear model
        orbit_number_after_sigma_clip = orbit_number[m]
        t_err_after_sigma_clip = t_err[m]
        t_after_sigma_clip = t[m]

        # re-fit linear model after enlarging uncertainties
        idxs, P, t0 = sigma_clip(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)
        linear_model = t0 + orbit_number_arr*P
        o_c = t - (t0 + P*orbit_number)
 
        idxs, P_quadratic, dpde_over_2, t0_quadratic = sigma_clip_quadratic(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)

        data = pd.read_csv(path + f'{planet_name}.csv')
        t_error = data['Uncertainty']

        fig, axs = plt.subplots(2)
        axs[0].plot(orbit_number_arr, ((t0_quadratic + P_quadratic*orbit_number_arr + dpde_over_2 * orbit_number_arr**2)-linear_model)*24*60, 'black')
        axs[0].errorbar(orbit_number[m], o_c[m]*24*60, yerr = t_error[m]*24*60, fmt='o', mew = 1, mfc = '#069AF3', ecolor = '#069AF3', alpha = 0.6,  markersize = 9)
        axs[0].errorbar(orbit_number[~m], o_c[~m]*24*60, yerr = t_error[~m]*24*60, fmt='o', mew = 1, mfc ='white', alpha = 0.8,  markersize = 9, ecolor = '#808080')#'#808080', alpha = 0.8,  markersize = 9, ecolor = '#808080')
     

        mask = data['Source'] == 'ETD'
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='red', alpha = 0.6, markersize = 9)

        mask = data['Source'] == 'Our work'  
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 9)

        mask = data['Source'] == 'Our work (sector 1)' 
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 9)

        mask = data['Source'] == 'Our work (sector 2)' 
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 9)
        axs[0].plot(t_arr, zero_line, '--', color = 'gray')
 
        axs[0].set_ylabel(r'$\Delta t$ [min]',  fontsize=18)
        axs[0].set_title(f'{planet_name.upper()}')
        axs[0].set_xlim(t_min - 50, t_max + 50)
 

        mask = data['Source'] == 'Our work'  
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        t_min = np.min(orbit_number_etd)
        t_max = np.max(orbit_number_etd)
        t_arr = np.linspace(t_min - 5, t_max + 5, 100)
        zero_line = np.zeros(t_arr.shape[0])

        axs[1].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 10)

        mask = data['Source'] == 'Our work (sector 1)' 
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[1].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 10)

        mask = data['Source'] == 'Our work (sector 2)' 
        o_c_etd = o_c[mask]
        orbit_number_etd = orbit_number[mask]
        uncertainty_etd = t_error[mask]

        axs[1].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', ecolor='green', alpha = 0.6, markersize = 10)
        axs[1].plot(t_arr, zero_line, '--', color = 'gray')

        axs[1].set_xlabel('Epoch',  fontsize=18)
        axs[1].set_ylabel(r'$\Delta t$ [min]',  fontsize=18)
        axs[1].set_xlim(t_min - 1, t_max + 1)
     

        plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/o_c_figures/o_c_quadratic/{planet_name}_o_c.png')
        #plt.show()
        plt.close(fig)


  except Exception as e: print(e)
 
