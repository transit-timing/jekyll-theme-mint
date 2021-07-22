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


# This script creates a histogram of the difference between MCMC TESS timing uncertainties and theoretical predictions

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels


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

chi_squares = []
ndofs = []

chi_squares_tess = []
ndofs_tess = []
 
uncs = []
estimated_uncs = []
systems = []

for folder in os.listdir('/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/'):
  try:
    planet_name = folder.replace('-mid-times', '')

    path = f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{folder}/'
    data = pd.read_csv(path + f'{planet_name}.csv')

    period = data['Period'].iloc[0]
    t = data['Mid-point']
    uncertainty = data['Uncertainty']
    mask = data['Source'] == 'Our work' 
    tess_uncertainty = uncertainty[mask]

    for file in os.listdir(f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{folder}/src/'):
        if file.find('results') != -1:
            results =  np.loadtxt(f'/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/{folder}/src/' + f'{file}', delimiter='\n', dtype=str)
            system = file.split("results.txt",1)[0]
            folder_name = system + 'dir'
            text = np.loadtxt(f'/Users/kate/tt/systems_merged/{folder_name}/' + system + 'log.txt',  delimiter='\n', dtype=str)
 
            stop = results == 'Result of linear fit to transit times using uncertainty method  0'
            start = results == 'Times with uncertainty method  0'

            idx2 = np.where(stop == True)[0][0]
            idx1 = np.where(start == True)[0][0]
            substring = results[idx1+1:idx2]

            num_mesurements = 0
            for measurement in substring:
                unc = measurement.split()[1]
                uncs.append(float(unc)*24*60)
                num_mesurements += 1
            num_estimates = 0
            for line in text:
                if line.find('estimated timing uncertainty [days,min] ') != -1:
                    unc_d_m = line.split("estimated timing uncertainty [days,min]   =  ",1)[1]

                    est_unc = unc_d_m.split()[1]
                    estimated_uncs.append(float(est_unc))
                    systems.append(system)

 
                    if float(unc)*24*60 - float(est_unc) > 10:
                
                        print(system)
                         

                    num_estimates += 1
            if num_mesurements != num_estimates:
                print('num_mesurements != num_estimates ', system)

  except:
    pass
print('estimated: ' ,len(estimated_uncs))
print('true: ', len(uncs))
estimated_uncs = np.array(estimated_uncs) 
uncs = np.array(uncs)
#np.savetxt('/Users/kate/Desktop/est.csv', estimated_uncs)
#np.savetxt('/Users/kate/Desktop/unc.csv', uncs)

diff = uncs - estimated_uncs
print(np.sum(diff>10))


# select mid-points with uncertainty < 10
mask = uncs < 10
uncs = uncs[mask]
estimated_uncs = estimated_uncs[mask]
#systems = systems[mask]
diff = diff[mask]

fig, axes = plt.subplots(1, 1, figsize=(22, 10))
#df = pd.DataFrame(data = {'system': systems, 'difference': diff, 'Theoretical uncertainty': estimated_uncs, 'MCMC uncertainty': uncs})
df = pd.DataFrame(data = {'difference': diff, 'Theoretical uncertainty': estimated_uncs, 'MCMC uncertainty': uncs})

straight_line = np.linspace(min(estimated_uncs), max(estimated_uncs), 100)

axes.plot(straight_line, straight_line, '-', color = 'gray')
#axes[0].plot(ndofs, y1, '.', color = 'gray')
#axes[0].plot(ndofs, y2, '.', color = 'gray')
#for term in terms:
#    axes[0].plot(ndofs_arr, term, '.', color = 'gray', alpha = 0.05)
 
#slope_arr = np.linspace(-3*np.sqrt(ndofs), 3*sqrt(ndofs), 100)
#for i in range(100):


axes.set_ylabel(r"$\delta_{MCMC}$ [min]")
axes.set_xlabel(r'$\delta_{theory}$[min]')
#axes[0].set_xlim(n_min, ndofs.max()+50)

axes.set_xscale('log')
axes.set_yscale('log')

#g = sns.histplot(data=df, x='estimated', bins=100) #r"$\delta_{MCMC} - \delta_{theory}$"
df.to_csv('/Users/kate/Desktop/df.csv')
g1 = sns.scatterplot(ax=axes, x="Theoretical uncertainty", y="MCMC uncertainty", data=df, edgecolor="black", linewidth=0.9) 
 
 
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/2_theoretical_uncs/2_mcmc_uncs_vs_theory_uncs.png')
plt.show()


#ch:r=-.5,l=.75

'''
plt.hist(diff, bins = 100)
plt.xlabel(r"$\delta_{MCMC} - \delta_{theory}$ [min]")
plt.ylabel('Count')
plt.xlim(diff.min()-2,10)
plt.savefig(f'/Users/kate/Documents/research/paper/2_figures/2_theoretical_uncs/1_hist_uncs.png')
plt.show()
''' 