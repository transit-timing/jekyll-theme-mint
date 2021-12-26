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
import emcee
from matplotlib.ticker import MaxNLocator
SMALL_SIZE = 24
MEDIUM_SIZE = 24
BIGGER_SIZE = 24

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels


path2database = f'/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/'


def round_to_uncertainty(value, uncertainty):
    # round the uncertainty to 1-2 significant digits
    significant_digit = 0
    for k in range(15):
      u = uncertainty * 10**k  
      if u > 1:
        significant_digit = k
        break
    unc, d = divmod(uncertainty*10**(significant_digit+1), 1)
    #print('significant_digit+1 = ', significant_digit+1)
    v, d = divmod(value*10**(significant_digit+1), 1)
    return str(v/10**(significant_digit+1)), str(int(unc)), significant_digit+1

def formatNumber(n, digits):
    formatter = formatter = '{:.' + '{}'.format(digits) + 'f}'
    x = round(n, digits)
    return formatter.format(x)

def sigma_clip(orbit_number, t, std, period_guess, t0_guess, dpdt_over_2_guess):
  params =  np.array([period_guess, t0_guess, dpdt_over_2_guess])
  pos = params + 1e-10 * np.random.randn(100, 3)
  nwalkers, ndim = pos.shape

  sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (orbit_number, t, std, t0_guess))
  sampler.run_mcmc(pos, 1000, progress=True)
  flat_samples = sampler.get_chain(discard=200, thin=15, flat=True)
  period_result = np.percentile(flat_samples[:, 0], [16, 50, 84])
  q = np.diff(period_result)
  #print('flat samples ', flat_samples)
  t0_result = np.percentile(flat_samples[:, 1], [16, 50, 84])
  q_t0 = np.diff(t0_result)

  dpde_result = np.percentile(flat_samples[:, 2], [16, 50, 84])
  q_dpde = np.diff(dpde_result)


  return period_result[1], t0_result[1], q[1], q_t0[1], dpde_result[1], q_dpde[1]


# Priors
def lnprior(theta, t0_init): 
  per, t0, dpdt_over_2 = theta
  if (t0_init - 0.25 < t0 < t0_init + 0.25):
    return 0
  return -np.inf


def lnlike(theta, orbit_number, t, sigma):
  per, t0, dpdt_over_2  = theta   
  model = t0 + orbit_number*per  + dpdt_over_2 * orbit_number**2
  inv_sigma2 = 1.0 / (sigma**2)
  return -0.5*(np.sum((t-model)**2*inv_sigma2))


# Define log of probability function.
def lnprob(theta, orbit_number, t, sigma, t0_init):
  lp = lnprior(theta, t0_init)
  if not np.isfinite(lp):
    return -np.inf
  return lp + lnlike(theta, orbit_number, t, sigma)

def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)

path2database = f'/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/'
#path2database = f'/Users/kate/Documents/research/paper/6_database/planets_>1_no_tess/'
ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/2_ephemerides/ephemerides.csv')

#for planet_name in os.listdir(path2database):
if True:
    planet_name = 'XO-3'
    path = path2database + planet_name + '/' #planets_>1_no_tess
    data = pd.read_csv(path + f'{planet_name}.csv')


    t = data['Mid-point']
    t_err = data['Uncertainty']
    t_err_original = np.copy(t_err)
    source = data['Source']

    t0_linear = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
    u_t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 error'].iloc[0]
    per_linear = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]
    u_per = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['P error'].iloc[0]
    t0_quadratic = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (quadraitc model)'].iloc[0]
    period_quadratic = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['P (quadraitc model)'].iloc[0]
    dpdt = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['dP/dE (quadraitc model) msec/yr'].iloc[0]
    scale = 86400000 * 365 / period_quadratic
    dpde_over_2 = dpdt/(2*scale)
    ###########################################
    # if uncertainty < 26 sec, set it to 30 sec
    small_uncertainty_mask = t_err < 0.0003

    t_err[small_uncertainty_mask] = 0.0003
    ###########################################

 
    orbit_number = fold_time(t, t0_linear, per_linear)
 
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)

    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    # re-fit linear model after enlarging uncertainties
    period, t0, P_err, t0_err, dpde_over_2, dpde_over_2_err = sigma_clip(orbit_number, t, t_err, period_quadratic, t0_quadratic, dpde_over_2)
    scale = 86400000 * 365 / period
    dpdt, dpdt_err = dpde_over_2 * 2 * scale, dpde_over_2_err * 2 * scale

    orbit_number_arr = np.linspace(min(orbit_number), max(orbit_number), 100)
    quadratic_model = dpde_over_2 * orbit_number_arr**2 +  period *orbit_number_arr + t0
    linear_model = per_linear*orbit_number_arr + t0_linear
    quadratic_residuals = quadratic_model - linear_model

    print(f'{planet_name}: P = {period} +- {P_err} and t0 = {t0} +- {t0_err} and dP/dE = {dpdt} +- {dpdt_err} [msec/yr]')


    o_c = t - (t0_linear + per_linear*orbit_number)

    mask = data['Source'] == 'TESS'  
    o_c_etd = o_c[mask]
    orbit_number_etd = orbit_number[mask]
    uncertainty_etd = t_err_original[mask]

  

    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)


    per, u_per, keep_index_per = round_to_uncertainty(per_linear, u_per)
    t0_val, u_t0, keep_index_t0 = round_to_uncertainty(t0_linear, u_t0)
    #print('P, per_unc: ', P, per_unc, 'per, u_per: ', per, u_per)
    per = formatNumber(per_linear, keep_index_per)
    t0 = formatNumber(t0_linear, keep_index_t0)
  


    # plot o-c diagram     
    fig, axs = plt.subplots(2)
    axs[0].plot(orbit_number_arr, quadratic_residuals*24*60, 'k')
    axs[0].errorbar(orbit_number[~mask], o_c[~mask]*24*60, yerr = t_err_original[~mask]*24*60, fmt='o', mew = 1, mfc = 'blue', alpha = 0.45,  markersize = 24, ecolor = 'black')
    
    axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='red', ecolor='black', alpha = 0.55, markersize = 18)
    axs[0].plot(t_arr, zero_line, '--', color = 'gray')


    axs[0].set_ylabel(r'$\Delta t$ [min]',  fontsize=24)
    axs[0].set_title(f'{planet_name}:' + rf' $P = {per}({u_per}), T_0 = {t0}({u_t0})$', fontsize=24)
    axs[0].set_xlim(t_min - 50, t_max + 50)
    axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))





    # bottom subplot (TESS data only)
    t_min = np.min(orbit_number_etd)
    t_max = np.max(orbit_number_etd)
    t_arr = np.linspace(t_min - 10, t_max + 10, 100)
    zero_line = np.zeros(t_arr.shape[0])

    axs[1].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='red', ecolor='black', alpha = 0.6, markersize = 18)
    axs[1].plot(t_arr, zero_line, '--', color = 'gray')

    axs[1].set_xlabel('Orbit Number',  fontsize=24)
    axs[1].set_ylabel(r'$\Delta t$ [min]',  fontsize=24)
    axs[1].set_xlim(t_min - 10, t_max + 10)
    axs[1].xaxis.set_major_locator(MaxNLocator(integer=True))


    plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/5_o_c_figures/o_c_quadratic/{planet_name}_o_c.png')
    #plt.show()
    plt.close(fig)



  #except Exception as e: print(e)
 
