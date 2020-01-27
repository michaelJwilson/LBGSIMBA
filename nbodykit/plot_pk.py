import numpy                              as np

import matplotlib; matplotlib.use('PDF')
import matplotlib.pyplot                  as plt

import pylab as pl

from   get_data import snaps


for x in snaps.keys():
  fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/pk_{:.5f}.txt'.format(x)
  k, P  = np.loadtxt(fpath, unpack=True)
    
  plt.semilogy(k, P, label=r'{:.2f}'.format(x), marker='^')

  # linear theory.
  k, P  = np.loadtxt('/home/mjwilson/LBGSIMBA/linpk/linpk_{:.5f}.txt'.format(x), unpack=True)

  plt.semilogy(k, P, label=r'{:.2f}'.format(x))
  
# format the axes                                                                                                                                                                                                                                                                                                       
plt.legend(loc=0, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.1,   0.6)
plt.ylim(1.e0, 2.e4)

plt.tight_layout()

pl.savefig('plots/pk.pdf'.format(x))


