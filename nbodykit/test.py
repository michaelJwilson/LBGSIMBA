import matplotlib

import pylab             as pl
import numpy             as np
import matplotlib.cm     as cm
import matplotlib.pyplot as plt

from   scipy             import interpolate


tracer         = 'g'  # ['g', 'dm']
x              = '2.02462'

k              = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_k.npy'.format(tracer, x))
mu             = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_mu.npy'.format(tracer, x))
Pk             = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_Pk.npy'.format(tracer, x))

mask           = np.isnan(k) | np.isnan(mu) | np.isnan(Pk) | 

k[mask]        = 1.e-9
mu[mask]       = 1.e-9
Pk[mask]       = 1.e-9

#f             = interpolate.interp2d(np.log10(k), mu, np.log10(Pk), kind='linear')

x              = np.arange(-0.3, 0.0, 1.5)
y              = np.arange(+0.0, 1.0, 1.5)

X, Y           = np.meshgrid(x, y)

Z              = f(X,Y)

fig, ax        = plt.subplots()
CS             = ax.contour(X, Y, X)

ax.clabel(CS, inline=1, fontsize=10)

pl.savefig('plots/contour.png')
