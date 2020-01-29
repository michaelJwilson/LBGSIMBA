import matplotlib; matplotlib.use('PDF')

import numpy.ma                           as     ma
import numpy                              as     np
import matplotlib.pyplot                  as     plt
import pylab                              as     pl

from   get_data                           import snaps
from   utils                              import latexify
from   scipy                              import interpolate
from   mpl_toolkits.axes_grid1            import make_axes_locatable



latexify(columns=1, equal=True, fontsize=6, ggplot=True, usetex=True)


colors         = plt.rcParams['axes.prop_cycle'].by_key()['color']

tracer         = 'g'  # ['g', 'dm']

for i, redshift in enumerate(snaps.keys()):
  pl.clf()

  k            = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{:.5f}_k.npy'.format(tracer,  redshift))
  mu           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{:.5f}_mu.npy'.format(tracer, redshift))
  Pk           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{:.5f}_Pk.npy'.format(tracer, redshift))

  k            = np.log10(k)
  Pk           = np.log10(Pk)
  
  k            = ma.MaskedArray( k, np.isnan(k)  | np.isinf(k)) 
  Pk           = ma.MaskedArray(Pk, np.isnan(Pk) | np.isinf(Pk))
  mu           = ma.MaskedArray(mu, np.isnan(mu) | np.isinf(mu))

  mask         = k.mask | mu.mask | Pk.mask

  k.mask       = mask
  mu.mask      = mask
  Pk.mask      = mask
    
  # plotting
  ax           = pl.gca()
  fig          = pl.gcf()

  ax.grid(False)
  
  x            = np.arange(-0.3, 0.0, 0.01)
  y            = np.arange( 0.0, 1.0, 0.01)

  X, Y         = np.meshgrid(x, y)

  pct          = 0.1
  
  vmin         = np.percentile(Pk.data[~mask].flatten(),  1)
  vmax         = np.percentile(Pk.data[~mask].flatten(), 99)

  print(vmin, vmax)
  
  Z            = interpolate.griddata((k.data[~mask], mu.data[~mask]), Pk.data[~mask], (X, Y), method='cubic')
  
  #  m         = ax.scatter(k.data[~mask], mu.data[~mask], c=Pk.data[~mask], vmin=vmin, vmax=vmax, s=4)
  m            = ax.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax)
  
  divider      = make_axes_locatable(ax)
  cax          = divider.append_axes("right", size="3%", pad=0.01)
  
  cb           = fig.colorbar(m, cax=cax)

  # labelpad=-.1
  cb.set_label(r'$\log_{\rm{10}}|P(k)|$', fontsize=8)
  
  ax.set_xlabel('$\log_{10}|k|$',         fontsize=8)
  ax.set_ylabel('$\mu$',                  fontsize=8)

  ax.set_xlim([-0.3, 0.0])
  
  # ax.margins(0.05)
  
  # ax.hold(True)
  
  # fig.tight_layout()

  pl.savefig('plots/ztdpk_{}_{}.pdf'.format(tracer, redshift))


