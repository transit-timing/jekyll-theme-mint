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
from decimal import Decimal
from scipy import optimize
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings("ignore")

MEDIUM_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE) 

# the name of the folder where O-C diagrams will be saved.


# if looking at targets with TESS data
output_folder = f'o_c'
path2database = f'/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/'
#path2database = f'/Users/kate/Documents/research/paper/6_database/planets_>1_no_tess/'

# if looking at targets without tess data but with >1 data point
#path2database = '/Users/kate/Documents/research/paper/6_database/planets_>1_no_tess/'
#output_folder = f'o_c_26102021_>1_no_tess'

 
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

nsigma_clip = 5
niter_sigmaclip = 10

def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)

def formatNumber(n, digits):
    formatter = formatter = '{:.' + '{}'.format(digits) + 'f}'
    x = round(n, digits)
    return formatter.format(x)
    
########################################################################################################

ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/2_ephemerides/ephemerides.csv')

planet_name = 'WASP-004'

path = path2database + planet_name + '/' #planets_>1_no_tess
data = pd.read_csv(path + f'{planet_name}.csv')


t = data['Mid-point']
t_err = data['Uncertainty']
source = data['Source']

t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
t0_unc = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 error'].iloc[0]
period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]
per_unc = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['P error'].iloc[0]
t0_quadratic = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (quadraitc model)'].iloc[0]
period_quadratic = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['P (quadraitc model)'].iloc[0]

orbit_number = fold_time(t, t0, period)

t_min = np.min(orbit_number)
t_max = np.max(orbit_number)

t_arr = np.linspace(t_min - 50, t_max + 50, 100)
zero_line = np.zeros(t_arr.shape[0])

orbit_number_arr = np.linspace(min(orbit_number), max(orbit_number), 100)

outlier = np.full(t.shape[0], False, dtype=bool)
idxs = []

for i in range(niter_sigmaclip):
    n_outliers_before = np.sum(outlier)

    t_model = period*orbit_number + t0
    t_res = t-t_model
    scatter = np.std(t_res[~outlier])
    outlier = (np.abs(t_res/scatter) > nsigma_clip)
    idx = np.where(outlier)[0]
    idxs.append(idx)

    if (np.sum(outlier) == n_outliers_before): break
idxs = np.concatenate(idxs).ravel()
idxs = np.unique(idxs)
m = np.ones(t.shape, bool) 
m[idxs[:]] = False

o_c = t - (t0 + period*orbit_number)

if np.sum(source == 'TESS') == 0:   
    y_min = (o_c[m] - t_err[m]).min() - 5/(60*24)
    y_max = (o_c[m]+t_err[m]).max() + 5/(60*24)
    plt.ylim(y_min*60*24, y_max*60*24)      
    n = o_c[(o_c  > y_max) | (o_c < y_min)].shape[0]
    if  n > 1:
      plt.annotate(f'{n} points are out of axis limits', xy=(0.05, 0.9), xycoords='axes fraction')
    elif n==1:
      plt.annotate(f'{n} point is out of axis limits', xy=(0.05, 0.9), xycoords='axes fraction')


    per, u_per, keep_index_per = round_to_uncertainty(period, per_unc)
    t0_val, u_t0, keep_index_t0 = round_to_uncertainty(t0, t0_unc)
    per = formatNumber(period, keep_index_per)
    t0 = formatNumber(t0, keep_index_t0)
  

    fig = plt.figure()

    #plt.errorbar(orbit_number[m], o_c[m]*24*60, yerr = t_err[m]*24*60, fmt='o', mew = 1, mfc = '#069AF3', alpha = 1,  markersize = 9, ecolor = '#069AF3')
    plt.errorbar(orbit_number[m], o_c[m]*24*60, yerr = t_err[m]*24*60, fmt='o', mew = 1, mfc = 'blue', alpha = 0.45,  markersize = 24, ecolor = 'black')
    plt.errorbar(orbit_number[~m], o_c[~m]*24*60, yerr = t_err[~m]*24*60, fmt='o', mew = 1, mfc = 'white', alpha = 0.8,  markersize = 24, ecolor = '#808080')
    
    # plot zero line
    t_min = np.min(orbit_number)
    t_max = np.max(orbit_number)
    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])
    plt.plot(t_arr, zero_line, '--', color = 'gray')
    plt.title(f'{planet_name}:' + rf' $P = {per}({u_per}), T_0 = {t0}({u_t0})$')
    plt.xlabel('Orbit Number')
    plt.ylabel(r'$\Delta t$ [min]')
    plt.xlim(t_min - 50, t_max + 50)
    plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/5_o_c_figures/{output_folder}/{planet_name}_o_c.png')
    #plt.show()
    plt.close(fig)

else:
    ##############################################################################################
    # tess transit times
    src =  source == 'TESS'
    mask = data['Source'] != 'TESS'
    o_c_etd = o_c[mask & m]
    orbit_number_etd = orbit_number[mask & m]
    uncertainty_etd = t_err[mask & m]

    fig, axs = plt.subplots(2)
    # upper subplot (all of the data: literature + TESS)
    #axs[0].errorbar(orbit_number[m], o_c[m]*24*60, yerr = t_err[m]*24*60, fmt='o', mew = 1, mfc = '#069AF3', ecolor = '#069AF3', alpha = 0.6,  markersize = 9)
    axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc = 'blue', ecolor = 'black', alpha = 0.45,  markersize = 18)        
    axs[0].errorbar(orbit_number[~m], o_c[~m]*24*60, yerr = t_err[~m]*24*60, fmt='o', mew = 1, mfc ='white', alpha = 0.8,  markersize = 18, ecolor = '#808080')#'#808080', alpha = 0.8,  markersize = 9, ecolor = '#808080')
 

    mask = data['Source'] == 'ETD'
    o_c_etd = o_c[mask]
    orbit_number_etd = orbit_number[mask]
    uncertainty_etd = t_err[mask]

    axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='green', alpha = 0.6, markersize = 18)

    mask = data['Source'] == 'TESS'  
    o_c_etd = o_c[mask]
    orbit_number_etd = orbit_number[mask]
    uncertainty_etd = t_err[mask]


    t_min = np.min(orbit_number[m])
    t_max = np.max(orbit_number[m])
    t_arr = np.linspace(t_min - 50, t_max + 50, 100)
    zero_line = np.zeros(t_arr.shape[0])

    per, u_per, keep_index_per = round_to_uncertainty(period, per_unc)
    t0_val, u_t0, keep_index_t0 = round_to_uncertainty(t0, t0_unc)

    per = formatNumber(period, keep_index_per)
    t0 = formatNumber(t0, keep_index_t0)
    
    #print('P, per_unc: ', P, per_unc, 'per, u_per: ', per, u_per)

    axs[0].errorbar(orbit_number_etd, o_c_etd*24*60, yerr = uncertainty_etd*24*60, fmt='o', mew = 1, mfc='red', ecolor='black', alpha = 0.55, markersize = 18)
    axs[0].plot(t_arr, zero_line, '--', color = 'gray')
    
    axs[0].set_ylabel(r'$\Delta t$ [min]',  fontsize=24)
    axs[0].set_title(f'{planet_name}:' + rf' $P = {per}({u_per}), T_0 = {t0}({u_t0})$',  fontsize=24)
     
    axs[0].set_xlim(t_min - 50, t_max + 50)
    axs[0].xaxis.set_major_locator(MaxNLocator(integer=True))

    # determine y-limits based on unclipped datapoints
    y_min = (o_c[m] - t_err[m]).min() - 0.5/(60*24)
    y_max = (o_c[m]+t_err[m]).max() + 0.5/(60*24)
    axs[0].set_ylim(y_min*60*24, y_max*60*24)
    n = o_c[(o_c  > y_max) | (o_c < y_min)].shape[0]
    if  n > 1:
      axs[0].annotate(f'{n} points are out of axis limits', xy=(0.05, 0.9), xycoords='axes fraction')
    elif n==1:
      axs[0].annotate(f'{n} point is out of axis limits', xy=(0.05, 0.9), xycoords='axes fraction')


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
 

    plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/5_o_c_figures/{output_folder}/{planet_name}_o_c.png')
    #plt.show()
    plt.close(fig)


  #except Exception as e: print(e)
 
