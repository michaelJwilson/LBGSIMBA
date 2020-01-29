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

for i, x in enumerate(['2.02462']):
  k            = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_k.npy'.format(tracer, x))
  mu           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_mu.npy'.format(tracer, x))
  Pk           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/ztdpk_{}_{}_Pk.npy'.format(tracer, x))

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
  
  x            = np.arange(-0.3, 0.0, 0.01)
  y            = np.arange( 0.0, 1.0, 0.01)

  X, Y         = np.meshgrid(x, y)

  pct          = 1
  
  vmin         = np.percentile(Pk.data[~mask].flatten(), pct)
  vmax         = np.percentile(Pk.data[~mask].flatten(), 100 - pct)

  print(vmin, vmax)
  
  Z            = interpolate.griddata((k.data[~mask], mu.data[~mask]), Pk.data[~mask], (X, Y), method='cubic')
  
  #  m         = ax.scatter(k.data[~mask], mu.data[~mask], c=Pk.data[~mask], vmin=vmin, vmax=vmax, s=4)
  m            = ax.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax)
  
  divider      = make_axes_locatable(ax)
  cax          = divider.append_axes("right", size="3%", pad=0.01)
  
  cb           = fig.colorbar(m, cax=cax)

  # cb.set_label(r'$\log_{\rm{10}}|P(k)|$', fontsize=6)

  ax.grid(False)
  
  ax.set_xlabel('$\log_{10}|k|$',         fontsize=8)
  ax.set_ylabel('$\mu$',                  fontsize=8)

  ax.set_xlim([-0.3, 0.0])
  '''
  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')
  
  ax.margins(0.05)
  '''
  # ax.hold(True)
  
  # fig.tight_layout()

pl.savefig('plots/ztdpk_{}.pdf'.format(tracer, x))


