import matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    hod               import  get_data


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244                                                                                                                                                                                                                                                                                                                                            
boxsize      = 100.
vol          = boxsize ** 3.

getredshift  = 3.00307

##                                                                                                                                                                                                                                                                                                                                                                                                     
print('\n\nWelcome to Simba UV luminosity.')

f, p         =  get_data(boxsize, getredshift)

##                                                                                                                                                                                                                                                                                                                                                                                                      
gid          =  f['galaxy_data']['GroupID'][:]
iscentral    =  f['galaxy_data']['central'][:]
haloindex    =  f['galaxy_data']['parent_halo_index'][:]

##  print(p['COLOR_INFO'][:])                                                                                                                                                                                                                                                                                                                                                                            
##                                                                                                                                                                                                                                                                                                                                                                                                     
cid          =  p['CAESAR_ID'][:]
LUMUV        =  p['absmag_19'][:]  ## 1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'                                                                                                                                                                                                                                                                                                   
lsstu        =  p['appmag_28'][:]
lsstg        =  p['appmag_29'][:]
lsstr        =  p['appmag_30'][:]
lssti        =  p['appmag_31'][:]
lsstz        =  p['appmag_33'][:]  ##  Note the switch in z and y ordering.                                                                                                                                                                                                                                                                                                                             
lssty        =  p['appmag_32'][:]

##  Assert on the CAESAR ordering:  i.e. equating GroupID to CAESAR_ID for galaxy catalogue.                                                                                                                                                                                                                                                                                                             
assert  np.all(gid == cid)

dMUV         =  0.5
bins         =  np.arange(-23., -11.5, dMUV)
blumuv       =  np.digitize(LUMUV, bins=bins)

ubins, cnts  =  np.unique(blumuv, return_counts = True)

assert  len(ubins) == (len(bins) - 1)

pl.plot(bins[:-1] + dMUV/2., np.log10(np.cumsum(cnts / vol)), c='darkcyan', lw=1, alpha=0.8)

pl.xlim(-22., -16.)
pl.ylim(-5., -1.5)

pl.xlabel(r'$M_{UV}$')
pl.ylabel(r'$\log_{10}|\bar n (M < M_{UV}) / (h^{-1} \rm{Mpc})^{-3}|$')

plt.tight_layout()

pl.savefig('plots/lumuv.pdf')
