import numpy as np
import pandas as pd
import os, sys, time
import warnings
import matplotlib.pyplot as plt 
from helpers import *
warnings.filterwarnings("ignore")
  
ephemerides_linear = []
ephemerides_quadratic = []
quadratic_fit_preferred = []
significant_dpde = [] 
chi2s = []

directories = ['planets_tess/tt_output_nov15']

ephemerides_linear_tess_targets = [] # only store P, P_err for TESS targets
planets = ['WASP-161', 'XO-6']

for directory in directories:
  for planet_name in planets:
    try:
      path = f'/Users/kate/Documents/research/paper/6_database/{directory}/{planet_name}/'
      data = pd.read_csv(path + f'{planet_name}.csv')
      data = data.sort_values(["Mid-point"], ascending=True)
      data = data.reset_index(drop=True)
      t = data['Mid-point']
      t_err = data['Uncertainty']
      source = data['Source']

      mask = source == 'TESS'
      t = t[mask]
      t_err = t_err[mask]

 
      # if uncertainty < 26 sec, set it to 30 sec
      small_uncertainty_mask = t_err < 0.0003
      t_err[small_uncertainty_mask] = 0.0003

      # choose the 0-th epoch to be closest the the middle of the time span of observations
      if directory == 'planets_>1_no_tess':
        t_arr = t.to_numpy()
        a0 = 0.5 * (np.min(t_arr) + np.max(t_arr))
        t0 = t_arr.flat[np.abs(t_arr - a0).argmin()]
        ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/init_ephemeris/2_table.csv')
      else:
        ephemerides =  pd.read_csv('/Users/kate/Desktop/nov15/tables/ephemerides.csv')
        t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]

      period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]
      orbit_number = fold_time(t, t0, period)

      if t.shape[0] > 2:
        idxs, P_init, t0_init, P_init_err, t0_init_err = sigma_clip(orbit_number, t, t_err)

        m = np.ones(t.shape, bool) 
        m[idxs[:]] = False


        # after sigma clipping with the linear model
        # asc = after sigma clipping
        orbit_number_asc = orbit_number[m]
        t_err_asc = t_err[m]
        t_asc = t[m]

        # if chi square > 1, inflate uncertainties to force chi2 = ndof
        model = t0_init + P_init * orbit_number_asc
        chi_square = np.sum((t_asc - model)**2/t_err_asc**2)
        ndof = t_asc.shape[0] - 2

        if ndof > 0:
          if chi_square/ndof > 1: 
            root = find_root(t_asc, model, t_err_asc)
          else: 
            root = 0
        else:
          root = 0     

        # enlarge uncertainty
        t_err_enlarged_asc = (t_err_asc**2 + root**2)**0.5

        # re-fit linear model after enlarging uncertainties
        idxs, P, t0, t0_err, P_err = sigma_clip(orbit_number_asc, t_asc, t_err_enlarged_asc)
        linear_model = t0 + orbit_number_asc*P
        chi2_lit = np.sum((t_asc - linear_model)**2/t_err_asc**2)

        chi2s.append([planet_name, chi2_lit, t_asc.shape[0]-2, chi2_lit/(t_asc.shape[0]-2)])
        ephemerides_linear.append([planet_name, t_asc.shape[0], chi2_lit, root, P, P_err, t0, t0_err])

        if directory == 'planets_tess/tt_output_nov15':
          ephemerides_linear_tess_targets.append([planet_name, t_asc.shape[0], chi2_lit, root, P, P_err, t0, t0_err])

        if ndof > 0:
          # fit quadratic model to transit times after sigma-clipping and after enlarging literature uncertainties
          idxs, P_quadratic, dpde_over_2, t0_quadratic, t0_quadratic_err, P_quadratic_err, dpde_over_2_err = sigma_clip_quadratic(orbit_number_asc, t_asc, t_err_enlarged_asc)
          scale = 86400000 * 365 / P_quadratic
          # check if quadratic model is preffered
          fit_all_quadratic = t0_quadratic + P_quadratic * orbit_number_asc + dpde_over_2 * orbit_number_asc**2 
          chi2_quadratic = np.sum((t_asc - fit_all_quadratic)**2/t_err_asc**2)

          bic_quadratic = chi2_quadratic + 3 *np.log(t_asc.shape[0])
          bic_linear = chi2_lit + 2 *np.log(t_asc.shape[0])

          if bic_quadratic < bic_linear:
            quadratic_fit_preferred.append([planet_name, bic_quadratic, bic_linear])

          # check if dP/dE is significant at 3-sigma level
          positive_derivative_at_3_sigma = (dpde_over_2 > 0) & (dpde_over_2 > 3 * dpde_over_2_err) # positive dP/dt cases at 3-sigma significance
          negative_derivative_at_3_sigma = (dpde_over_2 < 0) & (dpde_over_2 < - 3 * dpde_over_2_err) # negative dP/dt cases at 3-sigma significance
          if positive_derivative_at_3_sigma:
            significant_dpde.append([planet_name,  dpde_over_2 * 2 * scale, dpde_over_2_err * 2 * scale, 1])
          elif negative_derivative_at_3_sigma:
             significant_dpde.append([planet_name,  dpde_over_2 * 2 * scale, dpde_over_2_err * 2 * scale, -1])

          ephemerides_quadratic.append([planet_name, P_quadratic, P_quadratic_err, dpde_over_2 * 2 * scale, dpde_over_2_err * 2 * scale, t0_quadratic, t0_quadratic_err])
      else:
        slope = (t[1]-t[0])/(orbit_number[1]-orbit_number[0])
        slope_err = (t_err[1]**2 + t_err[0]**2)**0.5/(orbit_number[1]-orbit_number[0])
        b = t[1] - slope * orbit_number[1]
        b_err = (t_err[1]**2 + (slope_err*orbit_number[1])**2)**0.5
        ephemerides_linear.append([planet_name, 2, 0, 0, slope, slope_err, b, b_err])
        ephemerides_quadratic.append([planet_name, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
        chi2s.append([planet_name, np.nan, t.shape[0]-2, np.nan])

        if directory == 'planets_tess/tt_output_nov15':
          ephemerides_linear_tess_targets.append([planet_name, 2, 0, 0, slope, slope_err, b, b_err])
    
    except Exception as e: 
      print(e, ': ', planet_name)
      ephemerides_quadratic.append([planet_name, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

print(len(ephemerides_linear))
print(len(ephemerides_quadratic))
 
 
P_T0_linear = pd.DataFrame(columns = ['System', 'N', 'Chi-square', 'sigma_0', 'Period (days)', 'P error', 'T0 (BJD TDB)', 'T0 error'], data = ephemerides_linear)
P_T0_linear.to_csv('/Users/kate/Documents/research/paper/5_code/2_tables/2_ephemeris/wasp_161_xo_6/wasp_161_xo_6.csv', index=False)

print(P_T0_linear)