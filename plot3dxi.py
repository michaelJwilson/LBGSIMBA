import matplotlib; matplotlib.use('PDF')

import numpy                              as      np
import matplotlib.pyplot                  as      plt
import pylab                              as      pl

from   scipy                              import  interpolate
from   utils                              import  latexify
from   mpl_toolkits.axes_grid1            import  make_axes_locatable
from   scipy.ndimage                      import  gaussian_filter


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

ax         =  pl.gca()
fig        =  pl.gcf()

ax.grid(False)

bins       =  np.arange(0., 26., 1.)

rs         = (bins[:-1] + bins[1:]) / 2.
ps         =  np.arange(25.)

X, Y       =  np.meshgrid(rs, ps)

tracer     = 'g'

for space in ['z', '']:
  for redshift in [2.025, 3.003, 3.963, 5.024]: 
    _, pis, xi = np.loadtxt('dat/corrfunc3dxi_{}_{}_{:.3f}.txt'.format(tracer, space, redshift), unpack=True)

    X, Y       = np.meshgrid(rs, pis)

    Z          = interpolate.griddata((_, pis), xi, (X, Y), method='cubic')

    m          = ax.pcolormesh(X, Y, Z, vmin=-5., vmax=50., cmap='coolwarm')

    ##  CS     = ax.contour(X, Y, gaussian_filter(Z, sigma=1), levels=[6, 8], vmin=-2., vmax=10., colors='k')
    
    ax.set_xlabel(r'$r_p \ [h^{-1} \ \rm{Mpc}]$')
    ax.set_ylabel(r'$\pi \ [h^{-1} \ \rm{Mpc}]$')

    plt.xticks(np.arange(0, 26, 5))
    plt.yticks(np.arange(0, 26, 5))
    
    divider    = make_axes_locatable(ax)
    cax        = divider.append_axes("right", size="3%", pad=0.00)
   
    cb         = fig.colorbar(m, cax=cax, label=r'$\xi(r)$')

    plt.tight_layout()
    
    pl.savefig('plots/xi3d_{}_{}.pdf'.format(space, redshift))
