import numpy as np
import pandas as pd
import os, sys, time

 
observables = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/observables.csv')
 

systems = observables['System             '].str.lower() 
systems = systems.str.strip()
observables_p = observables[' Period(day)'].to_numpy()
observables_t0 = observables[' T0 (HJD or BJD)'].to_numpy()
observables_p_err = observables[' Perioderr  '].to_numpy() 
observables_t0_err = observables[' T0err    '].to_numpy() 
observables_ref = observables[' Ephemeris_reference'].str.strip().to_numpy() 

print(observables_ref)
planet_names = []
periods = []
t0s= []
period_errs = []
t0_errs = []
references = []
  
for folder in os.listdir('/Users/kate/Documents/research/paper/1_database/planets_1_data_point/'):
    planet_name = folder.replace('-mid-times', '')
    print(planet_name)
    current_name = systems == planet_name
    p = observables_p[current_name][0]
    t0 = observables_t0[current_name][0]
    p_err = observables_p_err[current_name][0]
    t0_err = observables_t0_err[current_name][0]
    ref = observables_ref[current_name][0]
 
     
    periods.append(p)
    t0s.append(t0)
    period_errs.append(p_err)
    t0_errs.append(t0_err) 
    references.append(ref)
    planet_names.append(planet_name.upper())
 
print(len(periods)) 
print(len(t0s)) 
ephemerides = pd.DataFrame(data={'Target': planet_names, 'Period': periods, 'Period error': period_errs, 'T0 (HJD or BJD)': t0s, 'T0 error': t0_errs,'Reference': references})

ephemerides.to_csv('/Users/kate/Documents/research/paper/3_tables/2_table_1_data_point.csv')
 