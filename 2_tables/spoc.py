import os
import numpy as np 
import pandas as pd 

# This script was used to create spoc.csv file which contains the list of targets analyzed in our work 
# for which we also analyzed TESS data

df = pd.read_csv('/Users/kate/Documents/research/paper/3_tables/1_target_list.csv')
ls = []

for file in os.listdir('/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean'):
	#print(file)
	path = '/Users/kate/Documents/research/paper/1_database/reviewed_planets_clean/' + file + '/src'
	for fname in os.listdir(path):
		if fname.endswith('.txt') and fname.find('results') != -1:
			print(file)
			ls.append(file)
			break
print(len(ls))
ls = np.array(ls)
np.savetxt('/Users/kate/Desktop/spoc.csv', ls, fmt='%s')