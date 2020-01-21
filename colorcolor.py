import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    get_data          import  get_data
from    utils             import  latexify
from    sphotometry       import  read_mags
from    fast_scatter      import  fast_scatter

##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]
##  Available snapshots: ['062',   '078',    '051',    '042']

redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, two,   two_nd    =  read_mags('078', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, three, three_nd  =  read_mags('062', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, four,  four_nd   =  read_mags('051', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, five,  five_nd   =  read_mags('042', infile=None, magcols=None, SUFF='app')

umg2         =  two['u'].values   - two['g'].values
gmr2         =  two['g'].values   - two['r'].values

umg3         =  three['u'].values - three['g'].values
gmr3         =  three['g'].values - three['r'].values

gmr4         =  four['g'].values  - four['r'].values 
rmi4         =  four['r'].values  - four['i'].values

imz5         =  five['i'].values  - five['z'].values
rmi5         =  five['r'].values  - five['i'].values

##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

fast_scatter(axes[0], gmr2, umg2, np.ones_like(gmr2), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[1], gmr3, umg3, np.ones_like(gmr3), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[2], rmi4, gmr4, np.ones_like(rmi4), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[3], imz5, rmi5, np.ones_like(imz5), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)

axes[0].set_xlabel(r'$g-r$')
axes[0].set_ylabel(r'$u-g$')

axes[1].set_xlabel(r'$g-r$')
axes[1].set_ylabel(r'$u-g$')

axes[2].set_xlabel(r'$r-i$')
axes[2].set_ylabel(r'$g-r$')

axes[3].set_xlabel(r'$i-z$')
axes[3].set_ylabel(r'$r-i$')

for ax in axes:
  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')

  ax.set_xlim(-.5,   .5)
  ax.set_ylim(-.5, 1.75)

  ax.legend(frameon=False, loc=1)
  
plt.tight_layout()

pl.savefig('plots/colorcolor.png')
