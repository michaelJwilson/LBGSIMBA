import matplotlib; matplotlib.use('PDF')

import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify


latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)


colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/abacus_pk_720.txt'
k, P  = np.loadtxt(fpath, unpack=True)

plt.loglog(k, P, marker='^', lw=0, markersize=3)

fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/abacus_pk_100.txt'
k, P  = np.loadtxt(fpath, unpack=True)
    
plt.loglog(k, P, marker='^', lw=0, markersize=3)

plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0$ [$h^{-3} \ \mathrm{Mpc}^3$]")

plt.xlim(0.01,  1.0)
plt.ylim(1.e1, 2.e4)

plt.tight_layout()

pl.savefig('plots/abacuspk_100.pdf')


