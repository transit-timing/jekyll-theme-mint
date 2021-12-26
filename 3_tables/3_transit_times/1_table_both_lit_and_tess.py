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

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 20

plt.rcParams["figure.figsize"] = (24,8)
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels

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
    return str(v/10**(significant_digit+1)), int(unc), significant_digit+1

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

ts = []
uncertainties = []
sources = []
planet_names = []
orbit_numbers = []
time_systems = []
notes = []
num_trs = []
#ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/v0/3_tables/2_table.csv')

ephemerides =  pd.read_csv('/Users/kate/Documents/research/paper/3_tables/2_ephemerides/ephemerides.csv')


for planet_name in os.listdir('/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/'):
  path = f'/Users/kate/Documents/research/paper/6_database/planets_tess/tt_output_nov15/{planet_name}/'
  data = pd.read_csv(path + f'{planet_name}.csv')
  data = data.sort_values(["Mid-point"], ascending=True)
  data = data.reset_index(drop=True)


  t = data['Mid-point']
  uncertainty = data['Uncertainty']
  source = data['Source']
  time_system = data['Time System']
  if 'Note' in data.columns:
    note = data['Note']
  elif 'Notes' in data.columns:
    note = data['Notes']
  else:
    note = np.full(t.shape[0], '')
  
  num_tr = data['#']

  t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
  period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]

  orbit_number = fold_time(t, t0, period)

  for i in range(t.shape[0]):
    midpoint = t.iloc[i]
    unc = uncertainty.iloc[i]

    t0_val, u_t0,  keep_index_t0 = round_to_uncertainty(midpoint, unc)
    t0_mid = formatNumber(midpoint, keep_index_t0)

    ts.append(float(t0_mid))
    uncertainties.append(float(u_t0/10**keep_index_t0))

    src = source.iloc[i]
    time_sys = time_system.iloc[i]
    time_systems.append(str(time_sys))
 
    sources.append(str(src))
    planet_names.append(planet_name)
    orbit_numbers.append(int(orbit_number[i]))
    notes.append(note[i])
    num_trs.append(num_tr.iloc[i])
      
  #except Exception as e: print(e, ': ', planet_name)
for planet_name in os.listdir('/Users/kate/Documents/research/paper/6_database/planets_>1_no_tess/'):
  path = f'/Users/kate/Documents/research/paper/6_database/planets_>1_no_tess/{planet_name}/'
  data = pd.read_csv(path + f'{planet_name}.csv')


  t = data['Mid-point']
  uncertainty = data['Uncertainty']
  source = data['Source']
  time_system = data['Time System']
  if 'Note' in data.columns:
    note = data['Note']
  elif 'Notes' in data.columns:
    note = data['Notes']
  else:
    note = np.full(t.shape[0], '')
  
  num_tr = data['#']

  t0 = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['T0 (BJD TDB)'].iloc[0]
  period = ephemerides[(ephemerides['System'].str.upper() == planet_name.upper())]['Period (days)'].iloc[0]

  orbit_number = fold_time(t, t0, period)

  for i in range(t.shape[0]):
    midpoint = t.iloc[i]
    unc = uncertainty.iloc[i]
 
    t0_val, u_t0,  keep_index_t0 = round_to_uncertainty(midpoint, unc)
    t0_mid = formatNumber(midpoint, keep_index_t0)

    ts.append(float(t0_mid))
    uncertainties.append(float(u_t0/10**keep_index_t0))

    src = source.iloc[i]
    time_sys = time_system.iloc[i]
    time_systems.append(str(time_sys))

    sources.append(str(src))
    planet_names.append(planet_name)
    orbit_numbers.append(int(orbit_number[i]))
    notes.append(note[i])
    num_trs.append(num_tr.iloc[i])

print('======================')
print(len(ts))
print(len(time_systems))
 
print('Total number of mid-transit measurements: ', len(planet_names)) 
print('Total number of papers: ', len(set(sources)))
print('Total number of targets: ', len(set(planet_names)))

df = pd.DataFrame(data={'System': planet_names, 'Orbit number': orbit_numbers, r'T_mid': ts, 'Uncertainty (days)': uncertainties, 'Time System': time_systems, '#': num_trs, 'Note': notes, 'Source': sources})
sorted_df = df.sort_values(["System"], ascending=True)
print(sorted_df)
sorted_df.to_csv('/Users/kate/Documents/research/paper/3_tables/3_transit_times/1_table.csv', index = False)

