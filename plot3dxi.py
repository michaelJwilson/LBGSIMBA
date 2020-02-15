import matplotlib; matplotlib.use('PDF')

import numpy                              as      np
import matplotlib.pyplot                  as      plt
import pylab                              as      pl

from   scipy                              import  interpolate
from   utils                              import  latexify
from   mpl_toolkits.axes_grid1            import  make_axes_locatable


latexify(columns=1, equal=True, fontsize=6, ggplot=True, usetex=True)

ax         = pl.gca()
fig        = pl.gcf()

ax.grid(False)

bins       = np.logspace(0., 1.47, 25)
rs         = (bins[:-1] + bins[1:]) / 2.
ps         = np.arange(25.)

X, Y       = np.meshgrid(rs, ps)

for redshift in [2.025]: 
    _, pis, xi = np.loadtxt('dat/corrfunc3dxi_{}.txt'.format(redshift), unpack=True)

    X, Y       = np.meshgrid(rs, pis)
    Z          = interpolate.griddata((_, pis), 1. + xi, (X, Y), method='cubic')

    m          = ax.pcolormesh(X, Y, Z, vmin=0., vmax=10.)

    ax.set_xlabel(r'$r_p \ [h^{-1} \ \rm{Mpc}]$')
    ax.set_ylabel(r'$\pi \ [h^{-1} \ \rm{Mpc}]$')
    
    divider    = make_axes_locatable(ax)
    cax        = divider.append_axes("right", size="3%", pad=0.01)
   
    cb         = fig.colorbar(m, cax=cax, label=r'$1 + \xi(r)$')
    
    pl.savefig('plots/xi3d_{}d.pdf'.format(redshift))
