import numpy                        as      np
import time
import astropy.io.fits              as      fits

from   get_data                     import  snaps

from   nbodykit.lab                 import  *
from   nbodykit.io.fits             import  FITSFile
from   nbodykit.io.csv              import  CSVFile
from   nbodykit.source.catalog.file import  CSVCatalog, FITSCatalog
from   nbodykit.transform           import  StackColumns
from   nbodykit                     import  setup_logging

##  source activate nbodykit-env
##  python3 nkit.py

compute          = True
hubble           = 0.68

t0               = time.time()

cols             = ['x', 'y', 'z']

for x in snaps.keys():
    if compute:
      setup_logging()
    
      ##  Comoving Mpc/h.         
      ##  cat         = CSVCatalog(fpath, names=cols, delim_whitespace=True)

      fpath           = '/home/mjwilson/LBGSIMBA/dat/simba_dmpos_{:.5f}.fits'.format(x) 
      
      cat             = FITSCatalog(fpath)
      cat['Position'] = StackColumns(cat['x'], cat['y'], cat['z'])

      ##  Box size in Mpc/h.
      mesh            = cat.to_mesh(resampler='CIC', Nmesh=16, compensated=False, BoxSize=1.e2 * hubble)
      
      r               = FFTPower(mesh, mode='2d', dk=0.05, kmin=0.1, Nmu=120, los=[0,0,1], poles=[0,2,4])

      poles           = r.poles

      t1              = time.time()

      print(t1 - t0)

      print('Number of galaxies: {}'.format(len(cat['x'])))
      print('Shotnoise: {}'.format(poles.attrs['shotnoise']))

      for ell in [0]:    
        k      = poles['k']
        P      = poles['power_%d' % ell].real

        if ell == 0:
          P    = P - poles.attrs['shotnoise']

        np.savetxt('dat/pk_{:.5f}.txt'.format(x), np.c_[k, P])
