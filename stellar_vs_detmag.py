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

fig, axes          =  plt.subplots(2, 3, figsize=(20, 5))

for i, (nodust, name) in enumerate(zip([True, False], ['nodust', 'dust'])):
  wave, two,   ids2    =  get_pyloser(boxsize, zs[0], nrows=nrows, snaps=snaps, nodust=nodust, steidel= True)
  wave, three, ids3    =  get_pyloser(boxsize, zs[1], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)
  wave, four,  ids4    =  get_pyloser(boxsize, zs[2], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)

  phys2                =     get_phys(boxsize, zs[0], nrows=nrows)
  phys3                =     get_phys(boxsize, zs[1], nrows=nrows)
  phys4                =     get_phys(boxsize, zs[2], nrows=nrows)
  
  ##
  bxdrops              =  insample_steidel('BX', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=True, default=True)
  bmdrops              =  insample_steidel('BM', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=True, default=True)

  ##                                                                                                                                                                                                                            
  udrops               =  insample('u',  three['LSST_u'].values, three['LSST_g'].values, three['LSST_r'].values, three['LSST_i'].values, three['LSST_z'].values, three['LSST_y'].values, maglim=True, default=True)
  gdrops               =  insample('g',   four['LSST_u'].values,  four['LSST_g'].values,  four['LSST_r'].values,  four['LSST_i'].values,  four['LSST_z'].values,  four['LSST_y'].values, maglim=True, default=True)

  axes[i][0].semilogy(two['steidel_rs'].values,          phys2['stellarmass'],          lw=0., marker='.', markersize=2, c='k')
  axes[i][0].semilogy(two['steidel_rs'].values[bxdrops], phys2['stellarmass'][bxdrops], lw=0., marker='.', markersize=2, c='gold', label='BX')
  axes[i][0].semilogy(two['steidel_rs'].values[bmdrops], phys2['stellarmass'][bmdrops], lw=0., marker='.', markersize=2, c='c', label='BM')

  axes[i][1].semilogy(three['LSST_i'].values,          phys3['stellarmass'],         lw=0., marker='.', markersize=2, c='k')
  axes[i][1].semilogy(three['LSST_i'].values[udrops],  phys3['stellarmass'][udrops], lw=0., marker='.', markersize=2, c='b', label=r'$u$')

  axes[i][2].semilogy(four['LSST_i'].values,          phys4['stellarmass'],         lw=0., marker='.', markersize=2, c='k')
  axes[i][2].semilogy(four['LSST_i'].values[gdrops],  phys4['stellarmass'][gdrops], lw=0., marker='.', markersize=2, c='g', label=r'$g$')
  
  print(np.count_nonzero(bxdrops), np.count_nonzero(bmdrops), len(bmdrops), np.count_nonzero(udrops), len(udrops), np.count_nonzero(gdrops), len(gdrops))
  
  for _ in range(3):
    axes[i][_].axhline(5.8e8, c='k', lw=1.0)
    axes[i][_].legend(frameon=False)

  axes[i][0].set_ylabel(r'$M_* \ \  [M_{\odot}]$')    
  
  axes[i][0].set_xlabel(r'$R_s$')
  axes[i][1].set_xlabel(r'$i$')
  axes[i][2].set_xlabel(r'$i$')

plt.tight_layout()

pl.savefig('plots/stellar_vs_detmag.png')
