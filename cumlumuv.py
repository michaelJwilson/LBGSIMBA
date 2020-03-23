import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  pandas             as      pd
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    scipy.spatial      import  KDTree
from    itertools          import  product
from    hod                import  get_data
from    get_data           import  snaps, get_pyloser
from    sphotometry        import  read_mags
from    insample           import  read_insample


dMUV         =  0.1
lims         =  np.arange(-23., -11.5, dMUV)

print('\n\nWelcome to Simba UV luminosity.')

colors       = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, axes    = plt.subplots(2, 1, figsize=(2.5, 5.))

##  Mpc/h.
for j, boxsize in enumerate([100.]):
 for i, (redshift, name, c) in enumerate(zip([2.024621, 3.00307, 3.963392, 5.0244], ['two', 'three', 'four', 'five'], colors)):
  for alpha, nodust in zip([1.0, 0.5], [True, False]):   
    wave, frame, ids = get_pyloser(boxsize, redshift, nrows=-1, magtype='abs', nodust=nodust)
  
    ##  Read sample selection.
    ##  lsst_sample  =  read_insample(redshift)
    ##  isin         =  lsst_sample['INSAMPLE'].values

    cum      =  np.array([np.count_nonzero(frame['i1500'] > x) for x in lims]) / (boxsize ** 3.)

    if nodust:
      label=r'$z={:.2f}$'.format(redshift)

    else:
      label=''
      
    axes[j].semilogy(lims, cum, c=c, alpha=alpha, label=label)
    axes[j].set_title('Boxsize: {}'.format(boxsize))
    
pl.legend(loc=2, frameon=False, handlelength=1)
  
pl.xlim(-22.,  -16.)
pl.ylim(1.e-4, 0.04)

pl.xlabel(r'$M_{UV}$')
pl.ylabel(r"$\bar n(> M_{UV})$")

pl.legend(loc=3, frameon=False)

plt.tight_layout()

pl.savefig('plots/cumlumuv.pdf')
