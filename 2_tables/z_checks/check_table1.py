import numpy as np
import pandas as pd
import os, sys, time

 
df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_target_list_tic_tmag_snr_mr_ratio.csv')
names = df['System'].str.lower() 

planet_names = []

count = 0
for folder in os.listdir('/Users/kate/Documents/research/paper/1_database/planets_1_data_point/'):
    planet_name = folder.replace('-mid-times', '')
    #names = names.to_numpy()
     
    if np.sum(np.isin(planet_name, names)) == 0:
    #if planet_name.upper() in names == False:
        print(planet_name)
        planet_names.append(planet_name)
        count += 1
print(count)
np.savetxt('check.txt', planet_names, fmt='%s')