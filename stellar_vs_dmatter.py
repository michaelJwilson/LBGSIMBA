import matplotlib;  matplotlib.use('PDF')

import numpy    as     np
import pylab    as     pl
import matplotlib.pyplot as plt

from   snaps    import snaps
from   get_data import get_pyloser, get_phys
from   insample import insample, insample_steidel


nrows              =   -1
boxsize            =  100.

zs                 =  np.sort(list(snaps.keys()))

print(zs)

fig, axes          =  plt.subplots(2, 3, figsize=(10, 10))

for i, (nodust, name) in enumerate(zip([True, False], ['nodust', 'dust'])):
  wave, two,   ids2    =  get_pyloser(boxsize, zs[0], nrows=nrows, snaps=snaps, nodust=nodust, steidel= True)
  wave, three, ids3    =  get_pyloser(boxsize, zs[1], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)
  wave, four,  ids4    =  get_pyloser(boxsize, zs[2], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)

  phys2                =     get_phys(boxsize, zs[0], nrows=nrows)
  phys3                =     get_phys(boxsize, zs[1], nrows=nrows)
  phys4                =     get_phys(boxsize, zs[2], nrows=nrows)
  phys5                =     get_phys(boxsize, zs[3], nrows=nrows)

  ##
  bxdrops              =  insample_steidel('BX', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=True, default=True)
  bmdrops              =  insample_steidel('BM', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=True, default=True)

  ##                                                                                                                                                                                                                            
  udrops               =  insample('u',  three['LSST_u'].values, three['LSST_g'].values, three['LSST_r'].values, three['LSST_i'].values, three['LSST_z'].values, three['LSST_y'].values, maglim=True, default=True)
  gdrops               =  insample('g',   four['LSST_u'].values,  four['LSST_g'].values,  four['LSST_r'].values,  four['LSST_i'].values,  four['LSST_z'].values,  four['LSST_y'].values, maglim=True, default=True)

  axes[i][0].loglog(phys2['hmass'],          phys2['stellarmass'],          lw=0., marker='.', markersize=2, c='k')
  axes[i][0].loglog(phys2['hmass'][bxdrops], phys2['stellarmass'][bxdrops], lw=0., marker='.', markersize=2, c='gold')
  axes[i][0].loglog(phys2['hmass'][bmdrops], phys2['stellarmass'][bmdrops], lw=0., marker='.', markersize=2, c='c')
  
  axes[i][1].loglog(phys3['hmass'],         phys3['stellarmass'],         lw=0., marker='.', markersize=2, c='k')
  axes[i][1].loglog(phys3['hmass'][udrops], phys3['stellarmass'][udrops], lw=0., marker='.', markersize=2, c='b')

  axes[i][2].loglog(phys4['hmass'],         phys4['stellarmass'],         lw=0., marker='.', markersize=2, c='k')
  axes[i][2].loglog(phys4['hmass'][gdrops], phys4['stellarmass'][gdrops], lw=0., marker='.', markersize=2, c='g')
  
  print(np.count_nonzero(bxdrops), np.count_nonzero(bmdrops), len(bmdrops), np.count_nonzero(udrops), len(udrops), np.count_nonzero(gdrops), len(gdrops))
  
  for _ in range(3):
    axes[i][_].axhline(5.8e8, c='k', lw=1.0)
    axes[i][_].legend(frameon=False)
    
    axes[i][_].set_xlabel(r'$M_h \ \  [M_{\odot}]$')

  axes[i][1].set_title('No dust? {}'.format(nodust))
  axes[i][0].set_ylabel(r'$M_* \ \  [M_{\odot}]$')

plt.tight_layout()

pl.savefig('plots/stellar_vs_halo.png')
