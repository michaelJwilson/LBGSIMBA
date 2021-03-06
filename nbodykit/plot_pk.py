import matplotlib; matplotlib.use('PDF')

import numpy                              as     np
import matplotlib.pyplot                  as     plt
import pylab                              as     pl

from   snaps                              import snaps
from   sutils                             import latexify


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']

zs       = [2.024621, 3.003070, 3.963392, 5.0244]
bs       = [2.200000, 3.200000, 4.200000, 5.5000]

insample = 1

for i, x in enumerate(zs):
  # Theory.                                                                                                                                                                                                                             
  iz      = int(100 * x + 0.001)
  k, P, H = np.loadtxt('/home/mjwilson/LBGSIMBA/linpk/dat/pklin_z{:03d}.txt'.format(iz), unpack=True)

  # Linear                                                                                                                                                                                                                              
  plt.loglog(k, P, label=r'', color=colors[i], alpha=0.3, lw=0.5)

  # Halofit                                                                                                                                                                                                                             
  plt.loglog(k, H, label=r'{:.2f}'.format(x), color=colors[i], lw=0.5)

  # Biased Halofit.                                                                                                                                                                                                                     
  plt.loglog(k, bs[i] * bs[i] * H, color=colors[i], lw=0.5)


  # Galaxies.
  fpath       = '/home/mjwilson/LBGSIMBA/nbodykit/dat/gpk_{:.5f}_insample_{}.txt'.format(x, insample)
  k, P, shot  = np.loadtxt(fpath, unpack=True)
  
  pl.axhline(shot[0], xmin=0, xmax=1, color=colors[i], alpha=0.6, linestyle='-')
  plt.semilogy(k, P, marker='^', color=colors[i], lw=0, markersize=3)

  # Cross
  fpath       = '/home/mjwilson/LBGSIMBA/nbodykit/dat/cpk_{:.5f}.txt'.format(x)
  k, P, shot  = np.loadtxt(fpath, unpack=True)

  plt.semilogy(k, np.abs(P), marker='x', color=colors[i], lw=0, markersize=4)
  
  # Dark Matter. 
  fpath       = '/home/mjwilson/LBGSIMBA/nbodykit/dat/dmpk_{:.5f}.txt'.format(x)
  k, P, shot  = np.loadtxt(fpath, unpack=True)
  
  pl.axhline(shot[0], xmin=0, xmax=1, color=colors[i], alpha=0.3)
  plt.semilogy(k, P, marker='s', color=colors[i], lw=0, markersize=2)
  
plt.legend(loc=3, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.07,  1.0)
plt.ylim(1.e0, 2.e4)

plt.tight_layout()

pl.savefig('plots/pk_insample_{}.pdf'.format(insample))


