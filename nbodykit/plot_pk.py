import matplotlib; matplotlib.use('PDF')

import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify


latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)

colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

zs      = [2.024621, 3.00307, 3.963392, 5.0244]

for i, x in enumerate(zs):
  fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/gpk_{:.5f}.txt'.format(x)
  k, P, shot  = np.loadtxt(fpath, unpack=True)

  
  pl.axhline(shot[0], xmin=0, xmax=1, color=colors[i], alpha=0.6, linestyle='-')
  plt.semilogy(k, P, marker='^', color=colors[i], lw=0, markersize=3)
  
  fpath       = '/home/mjwilson/LBGSIMBA/nbodykit/dat/dmpk_{:.5f}.txt'.format(x)
  k, P, shot  = np.loadtxt(fpath, unpack=True)

  pl.axhline(shot[0], xmin=0, xmax=1, color=colors[i], alpha=0.3)
  plt.semilogy(k, P, marker='s', color=colors[i], lw=0, markersize=3)
  
  # linear theory.
  iz      = int(100 * x + 0.001)
  k, P, _ = np.loadtxt('/home/mjwilson/LBGSIMBA/linpk/dat/pklin_z{:03d}.txt'.format(iz), unpack=True)
  
  plt.loglog(k, P, label=r'{:.2f}'.format(x), color=colors[i])
  
# format the axes                                                                                                                                                                                                                                                                                                       
plt.legend(loc=0, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.01,  1.0)
plt.ylim(1.e1, 2.e4)

plt.tight_layout()

pl.savefig('plots/pk.pdf'.format(x))


