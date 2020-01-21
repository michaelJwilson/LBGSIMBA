import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    hod               import  get_data
from    get_data          import  snaps
from    sphotometry       import  read_mags


boxsize      =  100.              ##  Mpc/h.
vol          =  boxsize ** 3.

dMUV         =  0.5
bins         =  np.arange(-23., -11.5, dMUV)

print('\n\nWelcome to Simba UV luminosity.')

for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
  ##  'UV':  1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'
  _, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, mag, mag_nd = read_mags(snaps[redshift], infile=None, magcols=None, SUFF='abs')
  
  blumuv       =  np.digitize(mag['UV'][mag['UV'] >= bins[0]], bins=bins, right=False)

  ubins, cnts  =  np.unique(blumuv, return_counts=True)

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
  
  pl.plot(bins + dMUV/2., np.log10(cnts / vol), lw=1, alpha=0.8, label='$z$ = %.2lf' % redshift)
  
pl.legend(loc=4, frameon=False, handlelength=1)
  
pl.xlim(-22., -16.)
pl.ylim(-5., -1.5)

pl.xlabel(r'$M_{UV}$')
pl.ylabel(r'$\log_{10}|d\bar n dM_{UV} / (h^{-1} \rm{Mpc})^{-3}|$')

plt.tight_layout()

pl.savefig('plots/lumuv.pdf')
