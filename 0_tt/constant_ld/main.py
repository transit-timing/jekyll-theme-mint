import sys
import os
import pandas as pd
import numpy as np
from astropy.io import fits
from tt import *


plt.rcParams['figure.figsize'] = (16.0, 8.0)
plt.rcParams['lines.markersize'] = np.sqrt(30)
plt.rcParams['font.size'] = 18
plt.rcParams['font.serif'] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
 
 

def main():


    # Open the DataFrame with the sample information
    #sample = pd.read_csv(os.path.dirname(os.getcwd()) + '/3_tables/1_target_list.csv')
    table = pd.read_csv('/scratch/gpfs/eivshina/tt2della/tables/1_target_list.csv')
    ephemerides =  pd.read_csv('/scratch/gpfs/eivshina/tt2della/nov15/tables/ephemerides.csv')
    df = pd.read_csv('/scratch/gpfs/eivshina/tt2della/nov15/tables/0_tess_merged.csv')
    
    omitted_orbits = [] # keep track whether MCMC converged, if there were orbits with 
    # too few points or if some orbits had high MAD

    i = int(sys.argv[-1]) 

    favorite = df['Target'].iloc[i]
    sector = df['Sector'].iloc[i]
    cadence = df['Cadence'].iloc[i] 
    sample = table[(table['System'].str.upper() == favorite.upper())]

    t0_updated = ephemerides[(ephemerides['Target'].str.upper() == favorite.upper())]['T0 (BJD TDB)'].iloc[0]
    period_updated = ephemerides[(ephemerides['Target'].str.upper() == favorite.upper())]['Period (days)'].iloc[0]
     
    #print('\nWorking on: ', favorite)
    #print('period: ', sample['period'].iloc[0])


    hdul = fits.open(f'/scratch/gpfs/eivshina/tt2della/nov15/data/{favorite}/{favorite}_{sector}.fits')
    data = hdul[1].data
    t = data['TIME']
    flux = data['FLUX']
                                            
    params = (favorite + f'_Sector_{sector}',
              period_updated,
              t0_updated,
              sample['duration'].iloc[0],
              sample['depth'].iloc[0]/100,
              cadence)

 
    mcmc_check, orbits_with_few_pts, orbits_with_high_mad, model_choice, ld1, ld2  = analyze_system(params, flux, t)
    omitted_orbits.append([favorite, sector, cadence, mcmc_check, orbits_with_few_pts, orbits_with_high_mad, model_choice, ld1, ld2])
    print(favorite, mcmc_check, orbits_with_few_pts, orbits_with_high_mad)

    #np.savetxt(f'/scratch/gpfs/eivshina/tt2della/tt_output_status_oct28/{favorite}.txt', omitted_orbits, fmt='%s')
    df = pd.DataFrame(columns = ['Target', 'Sector', 'Cadence', 'MCMC', 'Min pts', 'High Mad', 'Linear vs quadratic ld', 'LD1', 'LD2'],data = omitted_orbits)
    df.to_csv(f'/scratch/gpfs/eivshina/tt2della/nov15/status_subset/{favorite}_Sector_{sector}.csv', index=False)
           
if __name__ == "__main__":
    main()
