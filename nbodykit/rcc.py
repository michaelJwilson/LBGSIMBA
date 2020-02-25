import matplotlib; matplotlib.use('PDF')

import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

zs      = [2.024621, 3.003070, 3.963392, 5.0244]
bs      = [5.500000, 10.50000, 18.00000, 30.000]

for i, x in enumerate(zs):
  print('Solving for redshift: {:.2f}'.format(x))

  # Galaxies.
  fpath         = '/home/mjwilson/LBGSIMBA/nbodykit/dat/gpk_{:.5f}.txt'.format(x)
  gk, gP, shot  = np.loadtxt(fpath, unpack=True)

  # Add shotnoise back in. 
  gP           += shot
  
  print('Galaxies solved.')
  
  # Cross
  fpath         = '/home/mjwilson/LBGSIMBA/nbodykit/dat/cpk_{:.5f}.txt'.format(x)
  xk, xP, shot  = np.loadtxt(fpath, unpack=True)

  print('Cross solved.')
  
  # Dark Matter. 
  fpath         = '/home/mjwilson/LBGSIMBA/nbodykit/dat/dmpk_{:.5f}.txt'.format(x)
  dk, dP, shot  = np.loadtxt(fpath, unpack=True)

  dP           += shot
  
  print('Dark matter solved.')
  
  assert  np.all(gk == xk)
  assert  np.all(dk == xk)

  rcc           = xP / np.sqrt(gP * dP)

  pl.loglog(gk, gP, 'k-')
  pl.loglog(xk, xP, 'b-')
  pl.loglog(dk, dP, 'r-')
  
  # pl.semilogx(dk, rcc, label=r'$z={:.2f}$'.format(x))

  break
  
plt.legend(loc=1, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$r_{\rm{cc}}(k)$")

plt.xlim(0.01,  1.0)
#plt.ylim(0.00,  1.0)

plt.tight_layout()

pl.savefig('plots/rcc.pdf'.format(x))


