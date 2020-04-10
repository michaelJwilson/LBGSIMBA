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
from    cosmo             import  cosmo


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.

redshifts    = [3.00307, 3.963392, 5.0244]

latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

fig, axes    = plt.subplots(nrows=2, ncols=3, sharey=False, figsize=(10, 7))
plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

for j, nodust in enumerate([False, True]):
 for i, (zz, band) in enumerate(zip(redshifts, ['i', 'i', 'z'])):
  print(zz, band)
  
  band              = 'LSST_{}'.format(band)

  chi               =  cosmo.comoving_distance(zz).value   ##  [Mpc]. 
  lumdist           = (1. + zz) * chi            
  distmod           =  5. * np.log10(1.e5 * lumdist)

  fpath             = '/home/mjwilson/LBGSIMBA/pylosers/m100n1024/s50/pyloser_m100n1024_051.hdf5'
  fpath             =  None
  
  wave, frame, ids  =  get_pyloser(boxsize, zz, nrows=-1, nodust=nodust, steidel=False, magtype='abs', fpath=fpath, allfilters=True)
  MUV               =  frame['i1500'].values

  wave, frame, ids  =  get_pyloser(boxsize, zz, nrows=-1, nodust=nodust, steidel=False, magtype='app', fpath=fpath, allfilters=True)
  muv               =  frame[band].values

  print(np.count_nonzero(muv <= 25.8))
  
  exp               =  muv - distmod + 2.5 * np.log10(1. + zz)
  
  axes[j][i].plot(muv, exp - MUV, marker='.', c='k', lw=0, markersize=0.5, label=r'$z$={:.2f}, {}'.format(zz, 'no dust' if nodust else 'dust'))

  axes[j][i].set_xlabel(r'{}'.format(band[-1]) + r'$_{AB}$')
  axes[j][i].set_ylabel(r'$\langle M_{UV} \rangle - M_{UV}$')

  axes[j][i].set_xlim(21., 32.0)
  # axes[j][i].set_ylim(1.0,  2.0)

  axes[j][i].legend(loc=3, frameon=False)

  ax = axes[j][i]
  
  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')
  
plt.tight_layout()

pl.savefig('plots/kcorr.pdf')
