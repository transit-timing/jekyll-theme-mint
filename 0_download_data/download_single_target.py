import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import os 
import shutil

target = 'NGTS-11'
row = 4
#shutil.rmtree(f'/Users/kate/tt/data/{target}/')

dirname = '/Users/kate/tt/data/'+ target
if not os.path.exists(dirname):
  os.mkdir(dirname)

result = lk.search_lightcurve(f'{target}', mission = 'TESS')
sector = result[row].mission[0] 
 
s = int(sector.replace('TESS Sector ', ''))
print(f'{target}, Sector: {s}')
lc = result[row].download()
lc.to_fits(f'/Users/kate/tt/data/{target}/{target}_{s}.fits')

t = lc.time.to_value('jd', 'long')
f = lc.flux

plt.plot(t-2457000.0 , f, '.c')
plt.show()

hdul = fits.open(f'/Users/kate/tt/data/{target}/{target}_{s}.fits')
data = hdul[1].data
t = data['TIME']
f = data['FLUX']
plt.plot(t,f, '.')
plt.show()