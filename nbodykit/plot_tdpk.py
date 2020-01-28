import matplotlib; matplotlib.use('PDF')

import matplotlib.tri                     as tri
import numpy                              as np
import matplotlib.pyplot                  as plt
import pylab                              as pl

from   get_data                           import snaps
from   utils                              import latexify
from   scipy.interpolate                  import Rbf
from   mpl_toolkits.axes_grid1            import make_axes_locatable


latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)


colors         = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, x in enumerate(['2.02462']):
  k            = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/tdpk_{}_k.npy'.format(x))
  mu           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/tdpk_{}_mu.npy'.format(x))
  Pk           = np.load('/home/mjwilson/LBGSIMBA/nbodykit/dat/tdpk_{}_Pk.npy'.format(x))

  mask         = np.isnan(k) | np.isnan(mu) | np.isnan(Pk)

  k[mask]      = 1.e-9
  mu[mask]     = 1.e-9
  Pk[mask]     = 1.e-9
    
  # plotting
  fig, ax      = plt.subplots(1, 1)
  
  m            = ax.scatter(np.log10(k), mu, c=np.log10(Pk), vmin=0.5, vmax=3., cmap=plt.cm.Greens, s=4) 

  divider      = make_axes_locatable(ax)
  cax          = divider.append_axes("right", size="5%", pad=0.05)
  
  cb           = fig.colorbar(m, cax=cax)

  cb.set_label(r'$\log_{\rm{10}}|P(k)|$', fontsize='large')
  
  ax.set_xlabel('$\log_{10}|k|$',         fontsize='large')
  ax.set_ylabel('$\mu$',                  fontsize='large')

  ax.set_xlim([-0.5, 0.0])
  
  #ax.margins(0.05)

  ax.hold(True)
  
  fig.tight_layout()

pl.savefig('plots/tdpk.pdf'.format(x))


