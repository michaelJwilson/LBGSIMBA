import matplotlib;  matplotlib.use('PDF')

import pylab as pl
import numpy as np


files = ['pklin_z202.txt', 'pklin_z300.txt', 'pklin_z396.txt', 'pklin_z502.txt']

for file in files:
  zee            = file.split('.txt')[0].split('_')[-1]
    
  k, Plin, Pnlin = np.loadtxt('dat/{}'.format(file), unpack=True)

  pl.loglog(k, Plin, label=zee)

pl.savefig('plots/pklin.pdf')
