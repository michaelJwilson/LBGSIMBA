import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    get_data          import  get_data, get_pyloser, get_phys
from    utils             import  latexify
from    sphotometry       import  read_mags
from    fast_scatter      import  fast_scatter
from    color_box         import  color_box, color_box_steidel
from    insample          import  insample


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.
maglim       =  False
nodust       =  False

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]
##  Available snapshots: ['062',   '078',    '051',    '042']

redshifts    = [2.024621, 3.00307, 3.963392, 5.0244]

nrows        =  -1
prop         =  'hmass'  ## 'smass'

wave, two,   ids2   =  get_pyloser(boxsize, 2.024621, nrows=nrows, nodust=nodust, steidel=True)
wave, three, ids3   =  get_pyloser(boxsize, 3.003070, nrows=nrows, nodust=nodust)
wave, four,  ids4   =  get_pyloser(boxsize, 3.963392, nrows=nrows, nodust=nodust)
wave, five,  ids5   =  get_pyloser(boxsize, 5.024400, nrows=nrows, nodust=nodust)

phys2               =  get_phys(boxsize, 2.024621, nrows=nrows)['hmass']
phys3               =  get_phys(boxsize, 3.003070, nrows=nrows)['hmass']
phys4               =  get_phys(boxsize, 3.963392, nrows=nrows)['hmass']
phys5               =  get_phys(boxsize, 5.024400, nrows=nrows)['hmass']

##  umg2            =  two['LSST_u'].values   - two['LSST_g'].values
##  gmr2            =  two['LSST_g'].values   - two['LSST_r'].values

umg2                =  two['steidel_un'].values - two['steidel_g'].values
gmr2                =  two['steidel_g'].values  - two['steidel_rs'].values

umg3                =  three['LSST_u'].values - three['LSST_g'].values
gmr3                =  three['LSST_g'].values - three['LSST_r'].values

gmr4                =  four['LSST_g'].values  - four['LSST_r'].values 
rmi4                =  four['LSST_r'].values  - four['LSST_i'].values

imz5                =  five['LSST_i'].values  - five['LSST_z'].values
rmi5                =  five['LSST_r'].values  - five['LSST_i'].values

##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

if maglim == True:
  maglim2    = (two['steidel_rs'] > 23.5) & (two['steidel_rs'] < 25.5)
  maglim3    =  three['LSST_i'].values < 24.6
  maglim4    =   four['LSST_i'].values < 25.8
  maglim5    =   five['LSST_z'].values < 25.8

else:
   maglim2,  maglim3,  maglim4,  maglim5 = np.ones_like(two['LSST_u'].values, dtype=bool),\
                                           np.ones_like(three['LSST_u'].values, dtype=bool),\
                                           np.ones_like(four['LSST_u'].values, dtype=bool),\
                                           np.ones_like(five['LSST_u'].values, dtype=bool)
  
##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

cmap = 'viridis'

axes[0].scatter(gmr2[maglim2], umg2[maglim2], marker='.', c=np.log10(phys2[maglim2]), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[1].scatter(gmr3[maglim3], umg3[maglim3], marker='.', c=np.log10(phys3[maglim3]), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[2].scatter(rmi4[maglim4], gmr4[maglim4], marker='.', c=np.log10(phys4[maglim4]), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)
axes[3].scatter(imz5[maglim5], rmi5[maglim5], marker='.', c=np.log10(phys5[maglim5]), lw=0, s=3, cmap=cmap, vmin=9.5, vmax=14.)

color_box_steidel(axes[0], 'bx')
color_box(axes[1],         'u')
color_box(axes[2],         'g')
color_box(axes[3],         'r')

# axes[0].set_xlabel(r'$g-r$')
# axes[0].set_ylabel(r'$u-g$')

axes[0].set_xlabel(r'$G-R$')                                                                                                                                                                                                    
axes[0].set_ylabel(r'$U_n-G$')   

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

pl.savefig('plots/colorcolor_maglim_{:d}_nodust_{:d}.png'.format(maglim, nodust))
