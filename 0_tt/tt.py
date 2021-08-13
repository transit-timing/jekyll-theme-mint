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
import numdifftools as nd
from astroquery.mast import Observations
from astroquery.mast import Catalogs
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
import pickle5 as pickle
from scipy.signal import savgol_filter
import eleanor 
from astropy import units as u
from astropy.coordinates import SkyCoord
import lightkurve as lk
 
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
        
    return(coeff)
        
def transit_lightcurve(params, t):

    t0 = params['t0']
    period = params['period']
    u = params['limbdark']
    k = params['radratio']
    b = params['impactparam']
    a = params['a_over_r']
    c0 = params['c0']
    c1 = params['c1']

    phi = 2.*np.pi*(t-t0)/period
    x = a*np.sin(phi)
    y = b*np.cos(phi)
    s = np.sqrt(x*x + y*y)
    f_calc, f_no_ld = occultquad(s,u,0.0,k)

    f_calc = f_calc/(c0 + c1*(t-t0))
    
    return f_calc

def transit_lightcurve_residuals(params, t, f, f_unc):

    f_calc = transit_lightcurve(params, t)
    res = (f-f_calc)/f_unc
    return res

def transit_lightcurve_chisq(params, t, f, f_unc):

    f_calc = transit_lightcurve(params, t)
    res = (f-f_calc)/f_unc
    return np.sum(res*res)

def transit_lightcurve_residuals_30min(params, t, f, f_unc):

    ts = np.linspace(t.min() - 0.015, t.max() + 0.015, num = t.shape[0] * 10)
    f_calc = transit_lightcurve(params, ts)
    avg_f = np.zeros(f.shape[0])

    for i in range(t.shape[0]):
        bound1 = t[i] + 0.01
        bound2 = t[i] - 0.01
        m1 = ts < bound1
        m2 = ts > bound2
        mask = m1 | m2
        avg_f[i] = np.mean(f_calc[mask])

    res = (f-avg_f)/f_unc
    return res

def transit_lightcurve_chisq_30min(params, t, f, f_unc):
    f_calc = transit_lightcurve(params, t)
    res = (f-f_calc)/f_unc
    return np.sum(res*res)

    
'''
    ts = np.linspace(t.min() - 0.015, t.max() + 0.015, num = t.shape[0] * 10)
    f_calc = transit_lightcurve(params, ts)
    avg_f = np.zeros(f.shape[0])

    for i in range(t.shape[0]):
        bound1 = t[i] + 0.01
        bound2 = t[i] - 0.01
        m1 = ts < bound1
        m2 = ts > bound2
        mask = m1 and m2
        avg_f[i] = np.mean(f_calc[mask])

    res = (f-avg_f)/f_unc
    return np.sum(res*res)
'''
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

def mcmc(params, t, f, f_unc, n_links):

    print('  Starting MCMC with n_links = ', n_links)
    dt_unc_est = estimate_timing_uncertainty(params, len(t), f_unc[0])
    print('  estimated error = ', dt_unc_est)
    
    t0 = np.zeros(n_links)
    t0[0] = params['t0']
    chisq = np.zeros(n_links)
    chisq[0] = transit_lightcurve_chisq(params, t, f, f_unc)
    n_accept=0
    for i in range(1,n_links):
        t0_trial = t0[i-1] + 2.0*dt_unc_est*np.random.normal(0,1)
        params.add('t0', value=t0_trial, vary=True)
        chisq_trial = transit_lightcurve_chisq(params, t, f, f_unc)
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
#        print(t0[i])

    print('   Done, acceptance rate = ', n_accept/n_links)
    
    t0_center = np.mean(t0)
    t0_sigma = np.std(t0)
    
    return(t0_center, t0_sigma)

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
#        print(t0[i])

    print('   Done, acceptance rate = ', n_accept/n_links)
    
    t0_center = np.mean(t0)
    t0_sigma = np.std(t0)
    
    return(t0_center, t0_sigma)

def analyze_system(params_system, flux, time):


    # Unpack the parameters
    
    name, period_ref, t0_ref, duration, depth = params_system
    t0_ref = t0_ref - 2457000.0 # convert to TESS time
    
    # Set initial guesses and other parameters
    
    radratio = np.sqrt(depth)
    impactparam = 0.5
    limbdark = 0.5
    a_over_r = period_ref/(duration*np.pi)*np.sqrt(1.-impactparam**2)
    
    n_durations = 3
    oot_factor = 1.3
    min_npts_tra = 1#0.67*duration*24.*60./2. # 2/3 of the transit should be covered
    min_npts_oot = 1#0.5*duration*24.*60./2.
    nsigma_clip = 4
    method='nelder'
    n_links = 100 # if you are doing MCMC
    niter_sigmaclip = 5
    ndeg_max = 3 # maximum degree of polynomial for detrending

    # Make a directory for all the output
    

    dirname = '/Users/kate/tt/3_systems_lk/' + name + '_dir/'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        
    # Open the logfile
    
    original_stdout = sys.stdout
    sys.stdout = open(dirname + name + '_log.txt','w')
 
    print(name)
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print('Starting at ', date_time)

    # Open the key results file

    results_file = open(dirname + name + '_results.txt', 'w')

 
    ok = np.isfinite(flux)
    time=time[ok]
    flux=flux[ok]
 
    time_folded, orbit_number = fold_time(time, t0_ref, period_ref)
 
    orbits = np.unique(orbit_number) # a list of the orbits with valid timestamps

    oibeo = -1.0 + 0.0*time
    #  oibeo tells us:
    #   -1 is not analyzed, too far from transit
    #    0 is in-transit
    #    1 is pre-transit
    #    2 is post-transit
    #   >0 is OOT

    for n in orbits:

        this_orbit = (orbit_number == n) & (np.abs(time_folded) <= 0.5*n_durations*duration)
        pre = this_orbit & (time_folded <= -0.5*oot_factor*duration)
        post= this_orbit & (time_folded >=  0.5*oot_factor*duration)
        oot = pre | post
        tra = this_orbit & (np.abs(time_folded) <= 0.5*duration)
        n_pre = np.sum(pre)
        n_post = np.sum(post)
        n_tra = np.sum(tra)
        if (n_pre >= min_npts_oot) & (n_tra >= min_npts_tra) & (n_post >= min_npts_oot):
            oibeo[this_orbit] = 0
            oibeo[pre] = 1
            oibeo[post] = 2
        else:
            print('Rejecting data from orbit ', n)
            print('  n_tra, min_npts_tra = ', n_tra, min_npts_tra, ' and n_pre, n_post, min_npts_oot = ', n_pre, n_post, min_npts_oot)
            
    keep = (oibeo >= 0)
    if (np.sum(keep) == 0):
        print('No transits occurred during the timespan of TESS observations.\n')
        plt.figure(figsize=(10,10))
        plt.plot([0,0],[1,1],'w.')
        plt.text(0.5,0.5, name + ': no transits this sector')
        plt.savefig(dirname + name + '_NoData.png',bbox_inches='tight')
        plt.close()
        
        sys.stdout = original_stdout
        return
    
    t = time[keep]
    f = flux[keep]
    oibeo = oibeo[keep]
    orbit_number = orbit_number[keep]
    orbits = np.unique(orbit_number)
            
    # Plot the time series with transits identified

    print('\nPlotting the time series.')
    
    fig, ax = plt.subplots()
    plt.plot(time,flux,'c.',ms=5,alpha=0.2)
    plt.plot(t,f,'r.',ms=5,alpha=0.2)
    time_label = r'Time [BJD$_{\rm TDB}$]'
    plt.xlabel(time_label)
    plt.ylabel(r'Flux [e$^{-}$/s]')
    plt.title(name + ': PDCSAP Time Series')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_a_TimeSeries.png',bbox_inches='tight')
    plt.close()                    

    # Plot the individual transits and detrending function

    f_unc = 0.*t # we will set this equal to the std of oot after detrending
    
    n_cols = 3
    n_rows = np.int(np.round(len(orbits)/n_cols))+1
    fig = plt.figure(figsize=(15,3*n_rows))
    j=0
    for n in orbits:

        this_orbit = (orbit_number == n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        oot = (oibeo[this_orbit] > 0.)

        j=j+1
        plt.subplot(n_rows,n_cols,j)
        plt.plot(tp, fp, 'k.')

        tmid = np.mean(tp[oot])
        coeff = detrend(tp[oot]-tmid, fp[oot], ndeg_max)
        f_baseline = np.polyval(coeff,tp-tmid)
        plt.plot(tp[oot], f_baseline[oot], 'r.')

        fp = fp/f_baseline
        f[this_orbit] = fp
        scatter = np.std(fp[oot])
        f_unc[this_orbit] = scatter
        print('   Detrended orbit ', n, ' with polynomial of order ', len(coeff)-1, ', scatter = ', scatter)

        plt.title(str(int(n)))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    plt.savefig(dirname + name + '_b_IndividualTransits.png',bbox_inches='tight')
    plt.close()
    
    # Plot the individual transits after detrending

    fig = plt.figure(figsize=(15,3*n_rows))
    j=0
    for n in orbits:

        this_orbit = (orbit_number == n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        
        j=j+1
        plt.subplot(n_rows,n_cols,j)
        plt.plot(tp, fp, 'k.')
        plt.title(str(int(n)))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    plt.savefig(dirname + name + '_c_IndividualTransitsDetrended.png',bbox_inches='tight')
    plt.close()
     
    # Find a preliminary fit to folded light curve

    print('\nPlotting folded light curve.')
    
    t_folded, orbit_number = fold_time(t, t0_ref, period_ref)
    
    params = lmfit.Parameters()
    params.add('period', value=period_ref, vary=False)
    params.add('t0', value=0.0, vary=True, min=0, max = 1000000)
    params.add('radratio', value=radratio, vary=True, min=0.0)
    params.add('a_over_r', value=a_over_r, vary=True)
    params.add('impactparam', value=impactparam, vary=True, min=0.0, max=1.0)
    params.add('limbdark', value=limbdark, vary=True, min=0.0, max=1.0)
    params.add('c0', value=1.0, vary=False)
    params.add('c1', value=0.0, vary=False)
        
    print('\nFitting the folded light curve.\n')

    outlier = np.full(len(f), False, dtype=bool)

    mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                           fcn_args=(t_folded, f, f_unc))
    
    out = mini.minimize(method=method)
    popt = out.params
    f_calc = transit_lightcurve(popt, t_folded)
        
    print(lmfit.fit_report(out))

    # Plot the preliminary fit
    
    plt.figure(figsize=(15,10))

    ax = plt.subplot(2,1,1)
    plt.plot(t_folded, f, 'k.',alpha=0.3)
    #plt.plot(t_folded, f, '.', color='magenta', ms=15, alpha=0.5)
    plt.ylabel('Rectified Flux')
    plt.title('After rectification, before sigma clipping')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    f_calc_init = transit_lightcurve(popt, t_folded)
    q = np.argsort(t_folded)
    plt.plot(t_folded[q], f_calc_init[q], 'b-', linewidth=2)
    plt.plot(t_folded[q], f_calc[q], 'r-', linewidth=2)
    
    ax = plt.subplot(2,1,2)
    plt.plot(t_folded, f-f_calc, 'k.', alpha=0.3)
    #plt.plot(t_folded, f-f_calc, '.', color='magenta',ms=10, alpha=0.5)
    plt.xlabel('Folded time [days]')
    plt.ylabel('Residual flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_d_FoldedLightCurve.png',bbox_inches='tight')
    plt.close()
        
    # Find a fit after sigma-clipping

 

    print('\nFitting the folded light curve after sigma-clipping.\n')
    out = lmfit.minimize(transit_lightcurve_residuals, popt,
                   args=(t_folded, f, f_unc),
                   method=method)
    popt = out.params

    ######### If the fit fails, halt

    if (out.success == False):
        print('LM declares the fit a failure.')
        print('LM declares the fit a failure.', file=results_file)
        sys.stdout = original_stdout
        return
    
    ######### If the fit is bad, issue warning
    #print('ndof: ', out.nfree, ' num points: ', t_folded.shape[0])
    
    chisqr_red = out.chisqr/out.nfree
    chisqr_max = out.nfree + 3.*np.sqrt(2.*out.nfree)
    
    if (out.chisqr > chisqr_max):
        print('WARNING: Fit to folded light curve may be bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', out.chisqr, chisqr_max)
        print('WARNING: Fit to folded light curve may be bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', out.chisqr, chisqr_max, file=results_file)
    
    t0_revised = t0_ref + popt['t0']
    delta_t0 = t0_revised - t0_ref
    print(lmfit.fit_report(out))

    ########  Plot after sigma clipping
    
    plt.figure(figsize=(15,10))

    f_calc = transit_lightcurve(popt, t_folded)

    ax = plt.subplot(2,1,1)
    plt.plot(t_folded, f, 'k.', alpha=0.3)
    plt.ylabel('Rectified flux')
    plt.title('After rectification and sigma clipping')
    if (out.chisqr > chisqr_max): plt.title('After rectification and sigma clipping - BAD FIT?', color='red')
    q = np.argsort(t_folded)
    plt.plot(t_folded[q], f_calc[q], 'r-', linewidth=2)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    ax = plt.subplot(2,1,2)
    plt.plot(t_folded, f-f_calc, 'k.', alpha=0.3)
    plt.xlabel('Folded time [days]')
    plt.ylabel('Residual flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    yrmin, yrmax = ax.get_ylim()
    plt.savefig(dirname + name + '_e_FoldedLightCurve.png',bbox_inches='tight')
    plt.close()

    ######### Now for the main event, timing each transit

    print('\nTiming the transits.')

    orbits = np.unique(orbit_number)
    no = len(orbits)
    t_obs = np.zeros((3,no))
    t_unc = np.zeros((3,no))

    n_cols = 3
    n_rows = np.int(np.round(no/n_cols))+1
    
    fig = plt.figure(figsize=(15,3*n_rows))
    plt.axis('off')
    outer = gridspec.GridSpec(n_rows, n_cols, figure=fig)

    i = 0
    for n in orbits:

        print('\nWorking on transit serial number ', i, ', orbit number ', n) 
        this_orbit = (orbit_number==n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        fp_unc = f_unc[this_orbit]
                
        params = popt

        dt_unc_est = estimate_timing_uncertainty(params, len(tp), fp_unc[0])
        print('   estimated timing uncertainty [days,min]   = ', dt_unc_est, 24.*60.*dt_unc_est)
        
        t_n_guess = t0_revised + n*period_ref + dt_unc_est*np.random.normal(0.,1.)

        params = lmfit.Parameters()
        params.add('t0', value=t_n_guess, vary=True, min=t_n_guess-duration, max=t_n_guess+duration)
        params.add('period', value=popt['period'], vary=False)
        params.add('radratio', value=popt['radratio'].value, min=0.0, vary=False)
        params.add('a_over_r', value=popt['a_over_r'].value, vary=False)
        params.add('impactparam', value=popt['impactparam'].value, vary=False)
        params.add('limbdark', value=popt['limbdark'].value, vary=False)
        params.add('c0', value=1.0, vary=True)
        params.add('c1', value=0.0, vary=True)
        
        mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                       fcn_args=(tp, fp, fp_unc))

        out = mini.minimize(method=method)
        print(lmfit.fit_report(out))
        popt = out.params

        # Obtain results from LMFIT

        t_obs[0,i] = popt['t0'].value
        t_unc[0,i] = popt['t0'].stderr

        # Obtain results from LMFIT with CI routine

#        ci = lmfit.conf_interval(mini, out, sigmas=[1.0])
#        lmfit.printfuncs.report_ci(ci)
#
#        ci_t0 = ci['t0']
#        sigma, t0_lo = ci_t0[0]
#        sigma, t0_mid = ci_t0[1]
#        sigma, t0_hi = ci_t0[2]
#        t_obs[1,i] = 0.5*(t0_hi+t0_lo)
#        t_unc[1,i] = 0.5*(t0_hi-t0_lo)

#        chisqr = out.chisqr
#        chisqr_max = out.nfree + 3.*np.sqrt(2.*out.nfree)
#        if (out.success == False):
#            print('LM declares the fit a failure.')
#            print('LM declares the fit a failure.', file=results_file)
#            t_obs[i] = 0.0
#            t_unc[i] = -999
#        elif (chisqr > chisqr_max):
#            print('Fit is declared bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', chisqr, chisqr_max)
#            print('Fit is declared bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', chisqr, chisqr_max, file=results_file)
#            fac = np.sqrt(chisqr/chisqr_max)
#            print('Scaling up timing uncertainty by ', fac)
#            print('Scaling up timing uncertainty by ', fac, file=results_file)
#            t_unc[i] = t_unc[i] * fac
#            print('t0 = ' , t_obs[i] , ' +/- ', t_unc[i], file=results_file)
#            print('Ndata, Nparam, Ndof, chisqr = ' , out.ndata , out.nvarys , out.nfree, out.chisqr, '\n', file=results_file)
#        else:
#            print('t0 = ' , t_obs[i] , ' +/- ', t_unc[i], file=results_file)
#            print('Ndata, Nparam, Ndof, chisqr = ' , out.ndata , out.nvarys , out.nfree, out.chisqr, '\n', file=results_file)
#

        t0, t0_unc = mcmc(popt, tp, fp, fp_unc, n_links)
        print('   MCMC results for t0 = ', t0, t0_unc)
        t_obs[1,i] = t0
        t_unc[1,i] = t0_unc

        t0, t0_unc = mcmc(popt, tp, fp, fp_unc, 3*n_links)
        print('   MCMC-long results for t0 = ', t0, t0_unc)
        t_obs[2,i] = t0
        t_unc[2,i] = t0_unc

        # Plot the result using the LMFIT model

        tp_fit = np.linspace(tp.min(), tp.max(), num=500)
        
        fp_calc = transit_lightcurve(popt, tp_fit)

        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i])
        
        ax = plt.Subplot(fig, inner[0])
        ax.plot(tp-popt['t0'],fp,'.')
        ax.plot(tp_fit-popt['t0'],fp_calc,'r')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Orbit ' + str(np.int(n)))
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        fig.add_subplot(ax)

        fp_calc = transit_lightcurve(popt, tp)
    
        ax = plt.Subplot(fig, inner[1])
        ax.plot(tp-popt['t0'],fp-fp_calc,'.')
        ax.axhline(0,color='red')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(yrmin, yrmax)

        fig.add_subplot(ax) 
        i=i+1
        
    plt.savefig(dirname + name + '_f_IndividualTransitsWithFit.png',bbox_inches='tight')
    plt.close()
        
    #### Transit timing deviations

    print('\nCalculating transit timing deviations.')

    color = ['red','green','blue']
    for i in range(3):

        ok = np.isfinite(t_unc[i,:])
        nok = np.sum(ok)
        tn = t_obs[i,ok]
        tn_unc = t_unc[i,ok]
        n = orbits[ok]

        print('\nTimes with uncertainty method ', i, file=results_file)
        for j in range(len(tn)):
            print(tn[j] + 2457000.0, tn_unc[j], file=results_file)
        
        if (nok < 2.):
            print('\nNot enough timings to calculate ephemeris.')
            np.savez(dirname + name + '_results.npz', popt=popt, t=t_obs, t_unc=t_unc)
            sys.stdout = original_stdout
            return
        
        print('  Uncertainty method ', i, ': number of valid transit times = ', nok)
        ephemeris = np.polyfit(n, tn, 1, w=1./tn_unc)
        tn_calc = np.polyval(ephemeris, n)
        chisqr = np.sum( ((tn-tn_calc)/tn_unc)**2 )
        print('Chisqr, Ndof, Ndata, Np = ', chisqr, nok-2, nok, 2)
        print('\nResult of linear fit to transit times using uncertainty method ', i, file=results_file)
        print('  Period = ' , ephemeris[0], file=results_file)
        print('  t0     = ' , ephemeris[1] + 2457000, file=results_file)
        print('  Chisqr, Ndof, Ndata, Nparam = ', chisqr, nok-2, nok, 2, file=results_file)
        
        plt.errorbar(n+0.1*i,24.*60.*(tn-tn_calc),yerr=24.*60.*tn_unc,fmt='.',color=color[i],ms=10)
        
    plt.xlabel('Orbit Number')
    plt.ylabel('Timing deviation [min]')
    plt.savefig(dirname + name + '_g_TimingResiduals.png',bbox_inches='tight')
    plt.tight_layout()
    plt.close()

    sys.stdout = original_stdout

    return


def analyze_system_30min(params_system, flux, time):

    # Unpack the parameters
    
    name, period_ref, t0_ref, duration, depth = params_system
    t0_ref = t0_ref - 2457000.0 # convert to TESS time
    
    # Set initial guesses and other parameters
    
    radratio = np.sqrt(depth)
    impactparam = 0.5
    limbdark = 0.5
    a_over_r = period_ref/(duration*np.pi)*np.sqrt(1.-impactparam**2)
    
    n_durations = 3
    oot_factor = 1.3
    min_npts_tra = 1#0.67*duration*24.*60./2. # 2/3 of the transit should be covered
    min_npts_oot = 1#0.5*duration*24.*60./2.
    nsigma_clip = 4
    method='nelder'
    n_links = 100 # if you are doing MCMC
    niter_sigmaclip = 5
    ndeg_max = 3 # maximum degree of polynomial for detrending

    # Make a directory for all the output
    

    dirname = '/Users/Kate/tt/3_systems_lk/' + name + '_dir/'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        
    # Open the logfile
    
    original_stdout = sys.stdout
    sys.stdout = open(dirname + name + '_log.txt','w')
 
    print(name)
    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    print('Starting at ', date_time)

    # Open the key results file

    results_file = open(dirname + name + '_results.txt', 'w')


    ok = np.isfinite(flux)
    time=time[ok]
    flux=flux[ok]
    time_folded, orbit_number = fold_time(time, t0_ref, period_ref)
    orbits = np.unique(orbit_number) # a list of the orbits with valid timestamps

    oibeo = -1.0 + 0.0*time
    #  oibeo tells us:
    #   -1 is not analyzed, too far from transit
    #    0 is in-transit
    #    1 is pre-transit
    #    2 is post-transit
    #   >0 is OOT

    for n in orbits:

        this_orbit = (orbit_number == n) & (np.abs(time_folded) <= 0.5*n_durations*duration)
        pre = this_orbit & (time_folded <= -0.5*oot_factor*duration)
        post= this_orbit & (time_folded >=  0.5*oot_factor*duration)
        oot = pre | post
        tra = this_orbit & (np.abs(time_folded) <= 0.5*duration)
        n_pre = np.sum(pre)
        n_post = np.sum(post)
        n_tra = np.sum(tra)
        if (n_pre >= min_npts_oot) & (n_tra >= min_npts_tra) & (n_post >= min_npts_oot):
            oibeo[this_orbit] = 0
            oibeo[pre] = 1
            oibeo[post] = 2
        else:
            print('Rejecting data from orbit ', n)
            print('  n_tra, min_npts_tra = ', n_tra, min_npts_tra, ' and n_pre, n_post, min_npts_oot = ', n_pre, n_post, min_npts_oot)
            
    keep = (oibeo >= 0)
    if (np.sum(keep) == 0):
        print('No transits occurred during the timespan of TESS observations.\n')
        plt.figure(figsize=(10,10))
        plt.plot([0,0],[1,1],'w.')
        plt.text(0.5,0.5, name + ': no transits this sector')
        plt.savefig(dirname + name + '_NoData.png',bbox_inches='tight')
        plt.close()
        
        sys.stdout = original_stdout
        return
    
    t = time[keep]
    f = flux[keep]
    oibeo = oibeo[keep]
    orbit_number = orbit_number[keep]
    orbits = np.unique(orbit_number)
            
    # Plot the time series with transits identified

    print('\nPlotting the time series.')
    
    fig, ax = plt.subplots()
    plt.plot(time,flux,'c.',ms=5,alpha=0.7)
    plt.plot(t,f,'r.',ms=5,alpha=0.7)
    time_label = r'Time [BJD$_{\rm TDB}$]'
    plt.xlabel(time_label)
    plt.ylabel(r'Flux [e$^{-}$/s]')
    plt.title(name + ': Time Series')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_a_TimeSeries.png',bbox_inches='tight')
    plt.close()                    

    # Plot the individual transits and detrending function

    f_unc = 0.*t # we will set this equal to the std of oot after detrending
    
    n_cols = 3
    n_rows = np.int(np.round(len(orbits)/n_cols))+1
    fig = plt.figure(figsize=(15,3*n_rows))
    j=0
    for n in orbits:

        this_orbit = (orbit_number == n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        oot = (oibeo[this_orbit] > 0.)

        j=j+1
        plt.subplot(n_rows,n_cols,j)
        plt.plot(tp, fp, 'k.')

        tmid = np.mean(tp[oot])
        coeff = detrend(tp[oot]-tmid, fp[oot], ndeg_max)
        f_baseline = np.polyval(coeff,tp-tmid)
        plt.plot(tp[oot], f_baseline[oot], 'r.')

        fp = fp/f_baseline
        f[this_orbit] = fp
        scatter = np.std(fp[oot])
        f_unc[this_orbit] = scatter
        print('   Detrended orbit ', n, ' with polynomial of order ', len(coeff)-1, ', scatter = ', scatter)

        plt.title(str(int(n)))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    plt.savefig(dirname + name + '_b_IndividualTransits.png',bbox_inches='tight')
    plt.close()
    
    # Plot the individual transits after detrending

    fig = plt.figure(figsize=(15,3*n_rows))
    j=0
    for n in orbits:

        this_orbit = (orbit_number == n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        
        j=j+1
        plt.subplot(n_rows,n_cols,j)
        plt.plot(tp, fp, 'k.')
        plt.title(str(int(n)))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    plt.savefig(dirname + name + '_c_IndividualTransitsDetrended.png',bbox_inches='tight')
    plt.close()
     
    # Find a preliminary fit to folded light curve

    print('\nPlotting folded light curve.')
    
    t_folded, orbit_number = fold_time(t, t0_ref, period_ref)
    
    params = lmfit.Parameters()
    params.add('period', value=period_ref, vary=False)
    params.add('t0', value=-0, vary=True, min=-1, max = 1)
    params.add('radratio', value=radratio, vary=True, min=0.0)
    params.add('a_over_r', value=a_over_r, vary=True)
    params.add('impactparam', value=impactparam, vary=True, min=0.0, max=1.0)
    params.add('limbdark', value=limbdark, vary=True, min=0.0, max=1.0)
    params.add('c0', value=1.0, vary=False)
    params.add('c1', value=0.0, vary=False)
        
    print('\nFitting the folded light curve.\n')

     

    mini = lmfit.Minimizer(transit_lightcurve_residuals_30min, params,
                           fcn_args=(t_folded, f, f_unc))
    
    out = mini.minimize(method=method)
    popt = out.params
    f_calc = transit_lightcurve(popt, t_folded)
    print('my report')    
    print(lmfit.fit_report(out))

    # Plot the preliminary fit
    
    plt.figure(figsize=(15,10))

    ax = plt.subplot(2,1,1)
    plt.plot(t_folded, f, 'k.',alpha=0.3)
    #plt.plot(t_folded, f, '.', color='magenta', ms=15, alpha=0.5)
    plt.ylabel('Rectified Flux')
    plt.title('After rectification, before sigma clipping')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    f_calc_init = transit_lightcurve(popt, t_folded)
    q = np.argsort(t_folded)
    #plt.plot(t_folded[q], f_calc[q], 'r.', linewidth=2)
    plt.plot(t_folded[q], f_calc_init[q], 'r-', linewidth=2)
    
    ax = plt.subplot(2,1,2)
    plt.plot(t_folded, f-f_calc, 'k.', alpha=0.3)
    #plt.plot(t_folded, f-f_calc, '.', color='magenta',ms=10, alpha=0.5)
    plt.xlabel('Folded time [days]')
    plt.ylabel('Residual flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.savefig(dirname + name + '_d_FoldedLightCurve.png',bbox_inches='tight')
    plt.close()
   
        

    print('\nFitting the folded light curve after sigma-clipping.\n')
    out = lmfit.minimize(transit_lightcurve_residuals_30min, popt,
                   args=(t_folded, f, f_unc),
                   method=method)
    popt = out.params

    ######### If the fit fails, halt

    if (out.success == False):
        print('LM declares the fit a failure.')
        print('LM declares the fit a failure.', file=results_file)
        sys.stdout = original_stdout
     
    
    ######### If the fit is bad, issue warning
    
    chisqr_red = out.chisqr/out.nfree
    chisqr_max = out.nfree + 3.*np.sqrt(2.*out.nfree)
    
    if (out.chisqr > chisqr_max):
        print('WARNING: Fit to folded light curve may be bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', out.chisqr, chisqr_max)
        print('WARNING: Fit to folded light curve may be bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', out.chisqr, chisqr_max, file=results_file)
    
    t0_revised = t0_ref + popt['t0']
    delta_t0 = t0_revised - t0_ref
    print(lmfit.fit_report(out))

    ########  Plot after sigma clipping
    
    plt.figure(figsize=(15,10))

    f_calc = transit_lightcurve(popt, t_folded)

    ax = plt.subplot(2,1,1)
    plt.plot(t_folded, f, 'k.', alpha=0.3)
    plt.ylabel('Rectified flux')
    plt.title('After rectification and sigma clipping')
    if (out.chisqr > chisqr_max): plt.title('After rectification and sigma clipping - BAD FIT?', color='red')
    q = np.argsort(t_folded)
    plt.plot(t_folded[q], f_calc[q], 'r-', linewidth=2)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    ax = plt.subplot(2,1,2)
    plt.plot(t_folded, f-f_calc, 'k.', alpha=0.3)
    plt.xlabel('Folded time [days]')
    plt.ylabel('Residual flux')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    yrmin, yrmax = ax.get_ylim()
    plt.savefig(dirname + name + '_e_FoldedLightCurve.png',bbox_inches='tight')
    plt.close()

    ######### Now for the main event, timing each transit

    print('\nTiming the transits.')

    orbits = np.unique(orbit_number)
    no = len(orbits)
    t_obs = np.zeros((3,no))
    t_unc = np.zeros((3,no))

    n_cols = 3
    n_rows = np.int(np.round(no/n_cols))+1
    
    fig = plt.figure(figsize=(15,3*n_rows))
    plt.axis('off')
    outer = gridspec.GridSpec(n_rows, n_cols, figure=fig)

    i = 0
    for n in orbits:

        print('\nWorking on transit serial number ', i, ', orbit number ', n) 
        this_orbit = (orbit_number==n)
        tp = t[this_orbit]
        fp = f[this_orbit]
        fp_unc = f_unc[this_orbit]
                
        params = popt

        dt_unc_est = estimate_timing_uncertainty(params, len(tp), fp_unc[0])
        print('   estimated timing uncertainty [days,min]   = ', dt_unc_est, 24.*60.*dt_unc_est)
        
        t_n_guess = t0_revised + n*period_ref + dt_unc_est*np.random.normal(0.,1.)

        params = lmfit.Parameters()
        params.add('t0', value=t_n_guess, vary=True, min=t_n_guess-duration, max=t_n_guess+duration)
        params.add('period', value=popt['period'], vary=False)
        params.add('radratio', value=popt['radratio'].value, min=0.0, vary=False)
        params.add('a_over_r', value=popt['a_over_r'].value, vary=False)
        params.add('impactparam', value=popt['impactparam'].value, vary=False)
        params.add('limbdark', value=popt['limbdark'].value, vary=False)
        params.add('c0', value=1.0, vary=True)
        params.add('c1', value=0.0, vary=True)
        
        mini = lmfit.Minimizer(transit_lightcurve_residuals, params,
                       fcn_args=(tp, fp, fp_unc))

        out = mini.minimize(method=method)
        print(lmfit.fit_report(out))
        popt = out.params

        # Obtain results from LMFIT

        t_obs[0,i] = popt['t0'].value
        t_unc[0,i] = popt['t0'].stderr

        # Obtain results from LMFIT with CI routine

#        ci = lmfit.conf_interval(mini, out, sigmas=[1.0])
#        lmfit.printfuncs.report_ci(ci)
#
#        ci_t0 = ci['t0']
#        sigma, t0_lo = ci_t0[0]
#        sigma, t0_mid = ci_t0[1]
#        sigma, t0_hi = ci_t0[2]
#        t_obs[1,i] = 0.5*(t0_hi+t0_lo)
#        t_unc[1,i] = 0.5*(t0_hi-t0_lo)

#        chisqr = out.chisqr
#        chisqr_max = out.nfree + 3.*np.sqrt(2.*out.nfree)
#        if (out.success == False):
#            print('LM declares the fit a failure.')
#            print('LM declares the fit a failure.', file=results_file)
#            t_obs[i] = 0.0
#            t_unc[i] = -999
#        elif (chisqr > chisqr_max):
#            print('Fit is declared bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', chisqr, chisqr_max)
#            print('Fit is declared bad: Chisqr, Ndof+3*sqrt(2Ndof) = ', chisqr, chisqr_max, file=results_file)
#            fac = np.sqrt(chisqr/chisqr_max)
#            print('Scaling up timing uncertainty by ', fac)
#            print('Scaling up timing uncertainty by ', fac, file=results_file)
#            t_unc[i] = t_unc[i] * fac
#            print('t0 = ' , t_obs[i] , ' +/- ', t_unc[i], file=results_file)
#            print('Ndata, Nparam, Ndof, chisqr = ' , out.ndata , out.nvarys , out.nfree, out.chisqr, '\n', file=results_file)
#        else:
#            print('t0 = ' , t_obs[i] , ' +/- ', t_unc[i], file=results_file)
#            print('Ndata, Nparam, Ndof, chisqr = ' , out.ndata , out.nvarys , out.nfree, out.chisqr, '\n', file=results_file)
#

        t0, t0_unc = mcmc_30min(popt, tp, fp, fp_unc, n_links)
        print('   MCMC results for t0 = ', t0, t0_unc)
        t_obs[1,i] = t0
        t_unc[1,i] = t0_unc

        t0, t0_unc = mcmc_30min(popt, tp, fp, fp_unc, 3*n_links)
        print('   MCMC-long results for t0 = ', t0, t0_unc)
        t_obs[2,i] = t0
        t_unc[2,i] = t0_unc

        # Plot the result using the LMFIT model
        
        tp_fit = np.linspace(tp.min(), tp.max(), num=500)
        
        fp_calc = transit_lightcurve(popt, tp_fit)

        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i])
        
        ax = plt.Subplot(fig, inner[0])
        ax.plot(tp-popt['t0'],fp,'.')
        ax.plot(tp_fit-popt['t0'],fp_calc,'r')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Orbit ' + str(np.int(n)))
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        fig.add_subplot(ax)

        fp_calc = transit_lightcurve(popt, tp)
    
        ax = plt.Subplot(fig, inner[1])
        ax.plot(tp-popt['t0'],fp-fp_calc,'.')
        ax.axhline(0,color='red')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(yrmin, yrmax)

        fig.add_subplot(ax) 
        i=i+1
        
    plt.savefig(dirname + name + '_f_IndividualTransitsWithFit.png',bbox_inches='tight')
    plt.close()
        
    #### Transit timing deviations

    print('\nCalculating transit timing deviations.')

    color = ['red','green','blue']
    for i in range(3):

        ok = np.isfinite(t_unc[i,:])
        nok = np.sum(ok)
        tn = t_obs[i,ok]
        tn_unc = t_unc[i,ok]
        n = orbits[ok]

        print('\nTimes with uncertainty method ', i, file=results_file)
        for j in range(len(tn)):
            print(tn[j] + 2457000, tn_unc[j], file=results_file)
 
        
        print('  Uncertainty method ', i, ': number of valid transit times = ', nok)
        ephemeris = np.polyfit(n, tn, 1, w=1./tn_unc)
        tn_calc = np.polyval(ephemeris, n)
        chisqr = np.sum( ((tn-tn_calc)/tn_unc)**2 )
        print('Chisqr, Ndof, Ndata, Np = ', chisqr, nok-2, nok, 2)
        print('\nResult of linear fit to transit times using uncertainty method ', i, file=results_file)
        print('  Period = ' , ephemeris[0], file=results_file)
        print('  t0     = ' , ephemeris[1] + 2457000, file=results_file)
        print('  Chisqr, Ndof, Ndata, Nparam = ', chisqr, nok-2, nok, 2, file=results_file)
        
        plt.errorbar(n+0.1*i,24.*60.*(tn-tn_calc),yerr=24.*60.*tn_unc,fmt='.',color=color[i],ms=10)
        
    plt.xlabel('Orbit Number')
    plt.ylabel('Timing deviation [min]')
    plt.savefig(dirname + name + '_g_TimingResiduals.png',bbox_inches='tight')
    plt.tight_layout()
    plt.close()
    

    sys.stdout = original_stdout

    return

def main():

    # Open the DataFrame with the sample information
    sample = pd.read_csv(os.path.dirname(os.getcwd()) + '/3_tables/1_target_list.csv')
    
    #sample = pd.read_pickle("sample.pkl")
    path = os.path.dirname(os.getcwd()) + '/5_data/'
    favorite = 'XO-7'
    idx = 4 # row to be analyzed in the searchresult table



    sample = sample[(sample['System'] == favorite)]
     
    print('\nWorking on: ', favorite)
    print('period: ', sample['period'].iloc[0])

 
 
 

    result = lk.search_lightcurve(f'{favorite}', mission = 'TESS')
    #print(result)
    cadence = result.exptime[idx].value
    print('cadence: ', cadence)
    lc = result[idx].download()

    t = lc.time.to_value('jd', 'float')
    time = t - 2457000.0 
    flux = lc.flux.value

    #for i in range(flux.shape[0]):
    #    print('flux ', flux[i])
    nan_array = np.isnan(flux)
    not_nan_array = ~nan_array
    flux = flux[not_nan_array]
    time = time[not_nan_array]

    nsigma_clip = 15
    scatter = np.std(flux)
    print('scatter ', scatter)
    print('median ',  np.median(flux))

    mask = np.abs(flux - np.median(flux)) > nsigma_clip * scatter 

    fig, axs = plt.subplots(2)
    axs[0].plot(time, flux, '.c')
    axs[1].plot(time[~mask], flux[~mask], '.c')
    plt.show()
 
                                            
    params = (favorite + f'_{idx}',
              sample['period'].iloc[0],
              sample['t0'].iloc[0],
              sample['duration'].iloc[0],
              sample['depth'].iloc[0]/100.)

    flux = flux[~mask]
    time = time[~mask]
    #mask1 = time < 1610
   #mask2 = time < 1813
    #mask3 = time < 1789
    #mask = (mask1 | mask2) #& mask3
    #flux = flux[mask2]
    #time = time[mask2]


    if cadence == 120.0:
        print('analyzing 2-minute TESS data')
        analyze_system(params, flux, time)

    else:
        print('analyzing 30-minute TESS data')
        analyze_system_30min(params, flux, time)
      

       
if __name__ == "__main__":
    main()
