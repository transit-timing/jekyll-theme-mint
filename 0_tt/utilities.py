import sys
from datetime import datetime
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib.transforms import blended_transform_factory
from mandelagol import occultquad
from scipy import optimize
from scipy.optimize import minimize_scalar
#from lmfit import minimize, Parameters, fit_report
import lmfit as lmfit 



# This script contains helpher functions like fold time, mcmc, transit model etc


plt.rcParams['figure.figsize'] = (16.0, 8.0)
plt.rcParams['lines.markersize'] = np.sqrt(30)
plt.rcParams['font.size'] = 18
plt.rcParams['font.serif'] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
 
 
# Estimate the windowed scatter in a lightcurve
def estimate_scatter_with_mask(mask, flux):
    f = np.sum(flux[:, mask], axis=-1)
    #smooth data
    smooth = savgol_filter(f, 501, polyorder=5)
    return 1e6 * np.sqrt(np.median((f / smooth - 1) ** 2))

def fold_time(t, t0, period):
    
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero
 
    print(np.divmod( t-(t0-0.5*period), period ))
     

    orbit_number, t_fold = np.divmod( t-(t0-0.5*period), period )
    t_fold = t_fold - 0.5*period
    return(t_fold, orbit_number.astype(np.int))

def avg_y_in_bins_of_x(x,y,npts):

    # npts = how many points to average
    # no need to sort by x in advance
    
    nbins = np.int(len(x)/npts)
    xbin = np.zeros(nbins)
    ybin = np.zeros(nbins)
   
    i_sorted = np.argsort(x)
    x_sorted = x[i_sorted]
    y_sorted = y[i_sorted]
    
    i1 = 0
    i2 = npts
    for i in range(nbins):
        xbin[i] = np.mean(x_sorted[i1:i2])
        ybin[i] = np.mean(y_sorted[i1:i2])
        i1 = i2
        i2 = i1+npts
    
    return(xbin, ybin)

def detrend(t, f, deg_max):

    ndata = len(t)
    coeff = np.polyfit(t,f,1)
    f_calc = np.polyval(coeff,t)
    f_res = f-f_calc
    sigma = np.std(f_res)
    bic = np.sum((f_res/sigma)**2) + 1.0*np.log(ndata)
    
    if (deg_max > 1):
        for n in range(2,deg_max+1):
            coeff_trial = np.polyfit(t,f,n)
            f_calc = np.polyval(coeff_trial,t)
            f_res = f-f_calc
            bic_trial = np.sum((f_res/sigma)**2) + n*np.log(ndata)
            if (bic_trial < bic): coeff = coeff_trial
        
    return coeff
        
def transit_lightcurve(params, t):

    t0 = params['t0']
    period = params['period']
    u1 = params['limbdark1']
    u2 = params['limbdark2']
    k = params['radratio']
    b = params['impactparam']
    a = params['a_over_r']
 

    phi = 2.*np.pi*(t-t0)/period
    x = a*np.sin(phi)
    y = b*np.cos(phi)
    s = np.sqrt(x*x + y*y)
    f_calc, f_no_ld = occultquad(s,u1,u2,k)

 
    
    return f_calc

def transit_lightcurve_30min(params, t):
    ts = np.linspace(t.min() - 0.015, t.max() + 0.015, num = t.shape[0] * 100)
    f_calc = transit_lightcurve(params, ts)
    avg_f = np.zeros(t.shape[0])

    for i in range(t.shape[0]):
        bound1 = t[i] + 0.01
        bound2 = t[i] - 0.01
        m1 = ts < bound1
        m2 = ts > bound2
        mask = m1 & m2
        avg_f[i] = np.mean(f_calc[mask])
    return avg_f


def transit_lightcurve_residuals(params, t, f, f_unc, cadence):
    if int(cadence) == 120:       
        f_calc = transit_lightcurve(params, t)
        res = (f-f_calc)/f_unc
    else:
        avg_f = transit_lightcurve_30min(params, t)
        res = (f-avg_f)/f_unc        
    return res


def transit_lightcurve_chisq(params, t, f, f_unc, cadence):
    if int(cadence) == 120:  
        f_calc = transit_lightcurve(params, t)
        res = (f-f_calc)/f_unc
    else:
        avg_f = transit_lightcurve_30min(params, t)
        res = (f-avg_f)/f_unc            
    return np.sum(res*res)


 
def estimate_timing_uncertainty(params, n, sigma):
    
    k = params['radratio']
    b = params['impactparam']
    if (b>0.9): b=0.9
    a = params['a_over_r']
    per = params['period']
    durtot = per/a/np.pi*np.sqrt(1.0-b*b)
    Q = k*k/sigma*np.sqrt(n)
    durpar = 2.*k*per/a/np.pi/np.sqrt(1.0-b*b)
    theta = 2.*k/(1.0 - b*b)
    dt_unc_est = durtot/Q * np.sqrt(theta/2.)
    return dt_unc_est


def mcmc(params, t0_guess, t, f, f_detrended, f_unc, coeffs, oot, tmid, n_links, cadence, duration):

    print('  Starting MCMC with n_links = ', n_links)
    dt_unc_est = estimate_timing_uncertainty(params, len(t), f_unc[0])
    print('Theoretically  estimated error = ', dt_unc_est)
    
    t0 = np.zeros(n_links)
     
    coefficients = np.zeros((n_links, coeffs.shape[0]))
    for j in range(coeffs.shape[0]):
        coefficients[0, j] = coeffs[j]
    t0[0] = t0_guess
    print('starting MCMC with t0 = ', t0_guess)
    chisq = np.zeros(n_links)
    chisq[0] = transit_lightcurve_chisq(params, t, f_detrended, f_unc, cadence)
    n_accept=0
    for i in range(1,n_links):
        t0_trial = t0[i-1] + np.random.normal(0, 1/40 * duration)
        #print('t0_trial: ', t0_trial)
        coeffs_trial = coefficients[i-1, :] + np.random.normal(0,0.00025, size = coeffs.shape[0])

        f_baseline = np.polyval(coeffs_trial, t-tmid)
        fp = f/f_baseline
        scatter = np.std(fp[oot])
        f_unc = scatter
        params.add('t0', value = t0_trial, vary = True)
 
        chisq_trial = transit_lightcurve_chisq(params, t, fp, f_unc, cadence)
        del_chisq = chisq_trial-chisq[i-1]
        if (del_chisq < 0):
            n_accept = n_accept + 1
            chisq[i] = chisq_trial
            t0[i] = t0_trial
            for j in range(coeffs.shape[0]):
                coefficients[i, j] = coeffs_trial[j]
        else:
            prob = np.exp(-0.5*del_chisq)
            random_number = np.random.uniform(0,1)
            if (random_number <= prob):
                n_accept = n_accept + 1
                t0[i] = t0_trial
                chisq[i] = chisq_trial
                for j in range(coeffs.shape[0]):
                    coefficients[i, j] = coeffs_trial[j]
            else:
                chisq[i] = chisq[i-1]
                t0[i] = t0[i-1]
                for j in range(coeffs.shape[0]):
                    coefficients[i, j] = coefficients[i-1, j]


    print('   Done, acceptance rate 0 = ', n_accept/n_links)

    
    return t0, coefficients 


'''


def transit_lightcurve_residuals_30min(params, t, f, f_unc):

    avg_f = transit_lightcurve_30min(params, t)
    res = (f-avg_f)/f_unc
    return res


def transit_lightcurve_chisq_30min(params, t, f, f_unc):
    ts = np.linspace(t.min() - 0.015, t.max() + 0.015, num = t.shape[0] * 30)
    f_calc = transit_lightcurve(params, ts)
    avg_f = np.zeros(f.shape[0])

    for i in range(t.shape[0]):
        bound1 = t[i] + 0.01
        bound2 = t[i] - 0.01
        m1 = ts < bound1
        m2 = ts > bound2
        mask = m1 & m2
        #print('Number of points in a bin: ', f_calc[mask].shape[0])
        avg_f[i] = np.mean(f_calc[mask])

    res = (f-avg_f)/f_unc
    return np.sum(res*res)


def transit_lightcurve_chisq_30min(params, t, f, f_unc):
    avg_f = transit_lightcurve_30min(params, t)
    res = (f-avg_f)/f_unc
    return np.sum(res*res)




def mcmc_30min(params, t, f, f_unc, n_links):

    print('  Starting MCMC with n_links = ', n_links)
    dt_unc_est = estimate_timing_uncertainty(params, len(t), f_unc[0])
    print('  estimated error = ', dt_unc_est)
    
    t0 = np.zeros(n_links)
    t0[0] = params['t0']
    chisq = np.zeros(n_links)
    chisq[0] = transit_lightcurve_chisq_30min(params, t, f, f_unc)
    n_accept=0
    for i in range(1,n_links):
        t0_trial = t0[i-1] + 2.0*dt_unc_est*np.random.normal(0,1)
        params.add('t0', value=t0_trial, vary=True)
        chisq_trial = transit_lightcurve_chisq_30min(params, t, f, f_unc)
        del_chisq = chisq_trial-chisq[i-1]
        if (del_chisq < 0):
            n_accept = n_accept + 1
            chisq[i] = chisq_trial
            t0[i] = t0_trial
        else:
            prob = np.exp(-0.5*del_chisq)
            random_number = np.random.uniform(0,1)
            if (random_number <= prob):
                n_accept = n_accept + 1
                chisq[i] = chisq_trial
                t0[i] = t0_trial
            else:
                chisq[i] = chisq[i-1]
                t0[i] = t0[i-1]
 

    print('   Done, acceptance rate = ', n_accept/n_links)
    
    t0_center = np.mean(t0)
    t0_sigma = np.std(t0)
    
    return t0

'''