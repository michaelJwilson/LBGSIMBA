import matplotlib; matplotlib.use('PDF')

import numpy                              as      np
import matplotlib.pyplot                  as      plt
import pylab                              as      pl

from   scipy                              import  interpolate
from   utils                              import  latexify
from   mpl_toolkits.axes_grid1            import  make_axes_locatable


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

ax         =  pl.gca()
fig        =  pl.gcf()

ax.grid(False)

bins       =  np.arange(0., 26., 1.)

rs         = (bins[:-1] + bins[1:]) / 2.
ps         =  np.arange(25.)

X, Y       =  np.meshgrid(rs, ps)

for redshift in [2.025]: 
    _, pis, xi = np.loadtxt('dat/corrfunc3dxi_{}.txt'.format(redshift), unpack=True)

    X, Y       = np.meshgrid(rs, pis)

    Z          = interpolate.griddata((_, pis), xi, (X, Y), method='cubic')

    m          = ax.pcolormesh(X, Y, Z, vmin=-2., vmax=5., cmap='coolwarm')

    ax.set_xlabel(r'$r_p \ [h^{-1} \ \rm{Mpc}]$')
    ax.set_ylabel(r'$\pi \ [h^{-1} \ \rm{Mpc}]$')

    plt.xticks(np.arange(0, 26, 5))
    plt.yticks(np.arange(0, 26, 5))
    
    divider    = make_axes_locatable(ax)
    cax        = divider.append_axes("right", size="3%", pad=0.00)
   
    cb         = fig.colorbar(m, cax=cax, label=r'$\xi(r)$')

    plt.tight_layout()
    
    pl.savefig('plots/xi3d_{}d.pdf'.format(redshift))
