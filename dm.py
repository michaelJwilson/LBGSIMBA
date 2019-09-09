import numpy           as     np
import astropy.io.fits as     fits

from   readgadget      import *
from   astropy.table   import Table


snap   = '/home/rad/data/m100n1024/s50/snap_m100n1024_062.hdf5'
hubble = readheader(snap, 'h')

dmpos  = readsnap(snap, 'pos', 'dm', units=1, nth=32)    # in comoving kpc/h.
dmpos /= 1.e3                                            # in comoving Mpc/h.   

dmpos  = Table(dmpos, names=('x', 'y', 'z'))

dmpos['x'].unit = 'Mpc/h'
dmpos['y'].unit = 'Mpc/h'
dmpos['z'].unit = 'Mpc/h'

print(dmpos)

dmpos.write('simba_dmpos_3.00307.fits', format='fits', overwrite=True)

##  print(hubble)
##  print(dmpos)

##  np.savetxt('simba_dmpos_3.00307.txt', dmpos)
