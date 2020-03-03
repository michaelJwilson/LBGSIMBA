import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  pandas            as      pd
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    hod               import  get_data
from    get_data          import  snaps, get_pyloser
from    sphotometry       import  read_mags
from    insample          import  read_insample


boxsize      =  100.              ##  Mpc/h.
vol          =  boxsize ** 3.

dMUV         =  0.3
bins         =  np.arange(-23., -11.5, dMUV)

print('\n\nWelcome to Simba UV luminosity.')

colors       = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, (redshift, name) in enumerate(zip([2.024621, 3.00307, 3.963392, 5.0244], ['two', 'three', 'four', 'five'])):
  wave, frame, ids = get_pyloser(boxsize, redshift, nrows=-1, magtype='abs')
  
  ##  Read sample selection.
  lsst_sample   =  read_insample(redshift)
  isin          =  lsst_sample['INSAMPLE'].values
  
  for x, alpha, label in zip([frame['i1500'], frame['i1500'][isin]], [0.5, 1.0], ['', '$z$ = %.2lf' % redshift]):
    blumuv      =  np.digitize(x[x >= bins[0]], bins=bins, right=False)

    ubins, cnts =  np.unique(blumuv, return_counts=True)

    cnts         =   cnts[ubins < len(bins)]
    ubins        =  ubins[ubins < len(bins)]
    
    missing      =  np.setdiff1d(np.arange(len(bins)), ubins)

    ubins        =  np.concatenate([ubins, missing])
    cnts         =  np.concatenate([cnts,  np.zeros_like(missing)])

    sortind      =  np.argsort(ubins)    

    ubins        =  ubins[sortind]
    cnts         =  cnts[sortind]

    cnts         =  cnts / dMUV
  
    assert  len(ubins) == len(bins)
  
    pl.plot(bins + dMUV/2., np.log10(cnts / vol), lw=1, alpha=alpha, label=label, c=colors[i])

    
pl.legend(loc=2, frameon=False, handlelength=1)
  
pl.xlim(-22., -16.)
pl.ylim(-5., -1.5)

pl.xlabel(r'$M_{UV}$')
pl.ylabel(r'$\log_{10}|d\bar n/dM_{UV}|$')

plt.tight_layout()

pl.savefig('plots/lumuv.pdf')
