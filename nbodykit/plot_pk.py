import numpy                              as np

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot                  as plt

from   get_data import snaps


for x in snaps.keys():
  fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/pk_{:.5f}.txt'.format(x)
  k, P  = np.loadtxt(fpath, unpack=True)
    
  plt.semilogy(k, P, label=r'{:.2f}'.format(x), marker='-')

# format the axes                                                                                                                                                                                                                                                                                                       
plt.legend(loc=0, frameon=False)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.1,   0.6)
plt.ylim(1.e1, 1.e4)

plt.tight_layout()

pl.savefig('plots/pk.pdf'.format(x))


