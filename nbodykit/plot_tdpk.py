import matplotlib; matplotlib.use('PDF')

import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify
from   scipy.interpolate                  import Rbf


latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)


colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, x in enumerate(['tdpk_2.02462.txt']):
  fpath = '/home/mjwilson/LBGSIMBA/nbodykit/dat/{}'.format(x)
  kp, kt, P, N = np.loadtxt(fpath, unpack=True)

  sampled      = N > 1 

  kp           = kp[sampled][:30]
  kt           = kt[sampled][:30]
  P            =  P[sampled][:30]

  # https://stackoverflow.com/questions/33704428/plot-contours-of-a-given-set-of-points

  # initialize radial basis function
  rb = Rbf(kt, kp, P)

  # interpolate onto a 100x100 regular grid
  X, Y = np.mgrid[:2*np.pi:100j, :2*np.pi:100j]
  Z    = rb(X.ravel(), Y.ravel()).reshape(X.shape)
  
  # plotting
  fig, ax = plt.subplots(1, 1)
  ax.set_aspect('equal')
  ax.hold(True)

  m = ax.contourf(X, Y, Z, 20, cmap=plt.cm.Greens)

  ax.scatter(x, y, c=z, s=60, cmap=m.cmap, vmin=m.vmin, vmax=m.vmax)

  cb = fig.colorbar(m)

  cb.set_label('$f(x, y)$', fontsize='xx-large')
  ax.set_xlabel('$x$', fontsize='xx-large')

  ax.set_ylabel('$y$', fontsize='xx-large')

  ax.margins(0.05)

  fig.tight_layout()

pl.savefig('plots/tdpk.pdf'.format(x))


