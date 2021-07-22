import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import batman
import emcee
import corner
import os, sys, time

# sample_30_snr.csv was created using Sample Selection for Transit Timing.ipynb code. 
# This script checks which targets on the SNR list are not present on our final target list and visa versa. 


df1 = pd.read_csv('/Users/kate/Desktop/1_target_list.csv')
df2 = pd.read_csv('/Users/kate/Desktop/sample_25_snr.csv')


our_systems = df1['System'].to_list()
tepcat = df2['name'].str.strip().to_list()

#print(our_systems) 


print('target not on the snr list: ', list(set(our_systems) - set(tepcat)))


print('target not on our list: ', list(set(tepcat) - set(our_systems)))

print('target not on the snr list length: ', len(list(set(our_systems) - set(tepcat))))
print('target not on our list: ', len(list(set(tepcat) - set(our_systems))))

np.savetxt('/Users/kate/Desktop/target_not_on_snr_list.csv', list(set(our_systems) - set(tepcat)), fmt = '%s')
np.savetxt('/Users/kate/Desktop/target_not_on_our_list.csv', list(set(tepcat) - set(our_systems)), fmt = '%s')