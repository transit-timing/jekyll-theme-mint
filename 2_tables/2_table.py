import numpy as np
import pandas as pd
import os, sys, time

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
 
    return idxs, p, t0, t0_unc, per_unc

 
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
 
    return idxs, p, dpde_over_2, t0, t0_unc, per_unc, dpde_over_2_unc


def fold_time(t, t0, period):
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)
 


planet_names = []
t0s = []
periods = []
per_uncs = []
t0_uncs = []
chis = []
chis_tess = []
ndof_tess = []
ndofs = []
num_points = []
factors = []
dpdes = []
dpde_over_2_errs = []



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


    # choose the 0-th epoch to be closest the the middle of the time span of observations
    t_arr = t.to_numpy()
    a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
    t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
 
    orbit_number = fold_time(t, t0, period)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    orbit_number_arr = np.linspace(min(orbit_number), max(orbit_number), 100)


    idxs, P_init, t0_init, P_init_err, t0_init_err = sigma_clip(orbit_number, t, t_err)
    m = np.ones(t.shape, bool) 
    m[idxs[:]] = False

    # the case when there's no TESS data available
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
        idxs, P, t0, t0_err, P_err = sigma_clip(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)
        linear_model = t0 + orbit_number_arr*P
        o_c = t - (t0 + P*orbit_number)
 
        idxs, P_quadratic, dpde_over_2, t0_quadratic, t0_quadratic_err, P_quadratic_err, dpde_over_2_err = sigma_clip_quadratic(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)

        scale = 86400000 * 365 / P

        # append best-fitting parameters and their uncertainties
        planet_names.append(planet_name.upper())
        t0s.append(t0)
        periods.append(P)
        per_uncs.append(P_err)
        t0_uncs.append(t0_err)
        factors.append(factor)
        chis.append(chi_square)
        num_points.append(ndof + 2)
        ndofs.append(ndof)

        if ndof > 0:
            dpdes.append(2 * dpde_over_2 * scale)
            dpde_over_2_errs.append(dpde_over_2_err * scale)
        else:
            dpdes.append(np.nan)
            dpde_over_2_errs.append(np.nan)

    else:
        # the case when we have TESS data
        ##############################################################################################
        # tess transit times
        src =  source == 'Our work'
     
        mask = src.to_numpy() & m # mask to select TESS transit times after sigma-clipping was performed
        t_tess = t.to_numpy()[mask]
        uncertainty_tess = t_err.to_numpy()[mask]
        orbit_number_tess = orbit_number[mask]
     
        model_t_tess = t0_init + orbit_number_tess * P_init

        chi_tess = np.sum((t_tess - model_t_tess)**2/uncertainty_tess**2)

        ##############################################################################################
        # re-fit TESS data separately
        idxs_tess, P_tess, t0_tess, t0_err_tess, P_err_tess = sigma_clip(orbit_number.to_numpy()[src], t.to_numpy()[src], t_err.to_numpy()[src])
        fit_tess = t0_tess + orbit_number.to_numpy()[src] * P_tess
        chi_square_fit_tess = np.sum((t.to_numpy()[src] - fit_tess)**2/t_err.to_numpy()[src]**2)

        # append chi-square of TESS data and ndof for TESS data
        chis_tess.append(chi_square_fit_tess)
        ndof_tess.append(orbit_number.to_numpy()[src].shape[0] - 2)
        ##############################################################################################



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
        # calculate chi-square (literature and tess data after sigma clipping)
        chi_square = chi_lit + chi_tess

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
        idxs, P, t0, t0_err, P_err = sigma_clip(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)
        linear_model = t0 + orbit_number_arr*P
        o_c = t - (t0 + P*orbit_number)
 
        idxs, P_quadratic, dpde_over_2, t0_quadratic, t0_quadratic_err, P_quadratic_err, dpde_over_2_err = sigma_clip_quadratic(orbit_number_after_sigma_clip, t_after_sigma_clip, t_err_after_sigma_clip)
 
        scale = 86400000 * 365 / P

        # append best-fitting parameters and their uncertainties
        planet_names.append(planet_name.upper())
        t0s.append(t0)
        periods.append(P)
        per_uncs.append(P_err)
        t0_uncs.append(t0_err)
        factors.append(factor)
        chis.append(chi_square)
        num_points.append(ndof+2)
        ndofs.append(ndof)

        if ndof > 0:
            dpdes.append(2 * dpde_over_2 * scale)
            dpde_over_2_errs.append(dpde_over_2_err * scale)
        else:
            dpdes.append(np.nan)
            dpde_over_2_errs.append(np.nan)

  except Exception as e: print(e)
 


tess_chi_square = pd.DataFrame(data={'TESS chi-square': chis_tess, 'ndof': ndof_tess})
all_chi_square =  pd.DataFrame(data={'Chi-square': chis, 'ndof': ndofs})

tess_chi_square.to_csv('/Users/kate/Documents/research/paper/3_tables/tess_chi_square_vs_ndof.csv', index = False)
all_chi_square.to_csv('/Users/kate/Documents/research/paper/3_tables/all_chi_square_vs_ndof.csv', index = False)


df = pd.DataFrame(data={'Target': planet_names, 'T0 (BJD TDB)': t0s, r'T0 error': t0_uncs, 'Period (days)': periods, 'Period error': per_uncs, r'$\chi^2$': chis, 'Factor': factors, 'dP/dE': dpdes, 'dP/dE error': dpde_over_2_errs, 'Number of points': num_points})
df.to_csv('/Users/kate/Documents/research/paper/3_tables/2_table.csv', index = False)


