import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    get_data          import  get_data, get_pyloser
from    utils             import  latexify
from    sphotometry       import  read_mags
from    fast_scatter      import  fast_scatter
from    color_box         import  color_box
from    get_data          import  get_phys
from    insample          import  insample


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]
##  Available snapshots: ['062',   '078',    '051',    '042']

redshifts    = [2.024621, 3.00307, 3.963392, 5.0244]

nrows        =  -1
prop         =  'hmass'  ## 'smass'

wave, two    =  get_pyloser(boxsize, 2.024621, nrows=nrows)
wave, three  =  get_pyloser(boxsize, 3.003070, nrows=nrows)
wave, four   =  get_pyloser(boxsize, 3.963392, nrows=nrows)
wave, five   =  get_pyloser(boxsize, 5.024400, nrows=nrows)

phys2        =  get_phys(boxsize, 2.024621, nrows=nrows)['hmass']
phys3        =  get_phys(boxsize, 3.003070, nrows=nrows)['hmass']
phys4        =  get_phys(boxsize, 3.963392, nrows=nrows)['hmass']
phys5        =  get_phys(boxsize, 5.024400, nrows=nrows)['hmass']

##
umg2         =  two['LSST_u'].values   - two['LSST_g'].values
gmr2         =  two['LSST_g'].values   - two['LSST_r'].values

umg3         =  three['LSST_u'].values - three['LSST_g'].values
gmr3         =  three['LSST_g'].values - three['LSST_r'].values

gmr4         =  four['LSST_g'].values  - four['LSST_r'].values 
rmi4         =  four['LSST_r'].values  - four['LSST_i'].values

imz5         =  five['LSST_i'].values  - five['LSST_z'].values
rmi5         =  five['LSST_r'].values  - five['LSST_i'].values

##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

cmap = 'viridis'

axes[0].scatter(gmr2, umg2, marker='.', c=np.log10(phys2), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[1].scatter(gmr3, umg3, marker='.', c=np.log10(phys3), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[2].scatter(rmi4, gmr4, marker='.', c=np.log10(phys4), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[3].scatter(imz5, rmi5, marker='.', c=np.log10(phys5), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)

#fast_scatter(axes[0], gmr2, umg2, np.log10(phys2), 9.5, 14., 20, markersize=0.1, cmap=cmap, printit=False, alpha=1.0)
#fast_scatter(axes[1], gmr3, umg3, np.log10(phys3), 9.5, 14., 20, markersize=0.1, cmap=cmap, printit=False, alpha=1.0)
#fast_scatter(axes[2], rmi4, gmr4, np.log10(phys4), 9.5, 14., 20, markersize=0.1, cmap=cmap, printit=False, alpha=1.0)
#fast_scatter(axes[3], imz5, rmi5, np.log10(phys5), 9.5, 14., 20, markersize=0.1, cmap=cmap, printit=False, alpha=1.0)

color_box(axes[1], 'u')
color_box(axes[2], 'g')
color_box(axes[3], 'r')

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

  ax.set_xlim(-1.5, 2.00)
  ax.set_ylim(-2.0, 3.55)

  ax.legend(frameon=False, loc=1)

plt.tight_layout()

pl.savefig('plots/colorcolor.png')
