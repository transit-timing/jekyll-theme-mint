import numpy as np
import pandas as pd
import os, sys, time
import warnings
import matplotlib.pyplot as plt 
from helpers import * 
warnings.filterwarnings("ignore")
 
tess_chi2 = []
ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/2_ephemerides/ephemerides.csv')

for planet_name in os.listdir('/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15'):
  try:
    path = f'/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/{planet_name}/'
    data = pd.read_csv(path + f'{planet_name}.csv')

    t = data['Mid-point']
    t_err = data['Uncertainty']
    source = data['Source']
    original_errors = np.copy(t_err) # we won't modify this array of original errors

    # if uncertainty < 26 sec, set it to 30 sec
    small_uncertainty_mask = t_err < 0.0003
    t_err[small_uncertainty_mask] = 0.0003
 

    # choose the 0-th epoch to be closest the the middle of the time span of observations
    t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
    period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]
    orbit_number = fold_time(t, t0, period)

    # tess transit times
    src =  source == 'TESS'
    orbit_number_tess = orbit_number[src]
    t_tess =  t[src]
    t_err_tess = t_err[src]
    t_err_tess_original = original_errors[src]
 

    # here timing errors are floored at 26 sec; we fit all TESS data
    if t_tess.shape[0]>2:
      idxs_tess, P_tess, t0_tess, t0_err_tess, P_err_tess = sigma_clip(orbit_number_tess, t_tess, t_err_tess)
      # after sigma clipping with the linear model
      m_tess= np.ones(t_tess.shape, bool) 
      m_tess[idxs_tess[:]] = False
      #asc stands for after_sigma_clip
      orbit_number_tess_asc = orbit_number_tess[m_tess]
      t_err_tess_asc = t_err_tess[m_tess]
      t_err_tess_original_asc = t_err_tess_original[m_tess]
      t_tess_asc = t_tess[m_tess]

      fit_tess_asc = t0_tess + orbit_number_tess_asc * P_tess
      chi2_tess = np.sum((t_tess_asc - fit_tess_asc)**2/t_err_tess_original_asc**2)

      # append chi-square of TESS data and ndof for TESS data
      tess_chi2.append([planet_name, chi2_tess, t_tess_asc.shape[0]-2, chi2_tess/(t_tess_asc.shape[0]-2)])
    else:
      tess_chi2.append([planet_name, np.nan, t_tess.shape[0]-2, np.nan])

  except Exception as e: 
    print(e, ': ', planet_name)
 
 
tess_chi_square = pd.DataFrame(columns=['System', 'chi-square TESS','N dof', 'Reduced chi-square'], data = tess_chi2)
sorted_df = tess_chi_square.sort_values(["Reduced chi-square"], ascending=False)
sorted_df.to_csv('/Users/kate/Documents/research/paper/3_tables/0_ttv/tess_chi_square_vs_ndof.csv', index = False)
 