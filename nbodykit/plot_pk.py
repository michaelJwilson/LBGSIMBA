import matplotlib; matplotlib.use('PDF')

import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify


latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)


colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, x in enumerate(snaps.keys()):
  fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/pk_{:.5f}.txt'.format(x)
  k, P  = np.loadtxt(fpath, unpack=True)
    
  plt.semilogy(k, P, marker='^', color=colors[i], lw=0, markersize=3)

  # linear theory.
  iz      = int(100 * x + 0.001)
  k, P, _ = np.loadtxt('/home/mjwilson/LBGSIMBA/linpk/dat/pklin_z{:03d}.txt'.format(iz), unpack=True)
  
  plt.semilogy(k, P, label=r'{:.2f}'.format(x), color=colors[i])
  
# format the axes                                                                                                                                                                                                                                                                                                       
plt.legend(loc=0, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.1,   0.6)
plt.ylim(1.e0, 2.e4)

plt.tight_layout()

pl.savefig('plots/pk.pdf'.format(x))


