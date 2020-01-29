import pickle
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
##
##  export PYTHONPATH=/home/mjwilson/LBGSIMBA/:$PYTHONPATH
##
##  python3 nkit.py

compute          = True
hubble           = 0.68

t0               = time.time()

cols             = ['x', 'y', 'z']

tracer           = 'g' # ['dm', 'g']

for x in snaps.keys():
    if compute:
      setup_logging()
    
      ##  Comoving Mpc/h.         
      ##  cat         = CSVCatalog(fpath, names=cols, delim_whitespace=True)

      fpath           = '/home/mjwilson/LBGSIMBA/dat/simba_{}zpos_{:.5f}.fits'.format(tracer, x) 
      
      cat             = FITSCatalog(fpath)
      cat['Position'] = StackColumns(cat['x'], cat['y'], cat['z'])

      ##  Box size in Mpc/h.
      mesh            = cat.to_mesh(resampler='CIC', Nmesh=1024, compensated=False, BoxSize=1.e2)
      
      r               = FFTPower(mesh, mode='2d', dk=0.05, kmin=0.1, Nmu=60, los=[0,0,1], poles=[0,2,4], BoxSize=100., Nmesh=1024)

      poles           = r.poles

      t1              = time.time()

      print(t1 - t0)

      print('Number of galaxies: {}'.format(len(cat['x'])))
      print('Shotnoise: {}'.format(poles.attrs['shotnoise']))

      for _ in r.power.attrs:
          print("%s = %s" %(_, str(r.power.attrs[_])))

      # The 2D power spectrum results
      print("power = ", r.power)
      print("variables = ", r.power.variables)

      for name in r.power.variables:
          var = r.power[name]

          print("'%s' has shape %s and dtype %s" %(name, var.shape, var.dtype))

      np.save('dat/ztdpk_{}_{:.5f}_k.npy'.format(tracer, x),  r.power['k'])
      np.save('dat/ztdpk_{}_{:.5f}_mu.npy'.format(tracer, x), r.power['mu'])
      np.save('dat/ztdpk_{}_{:.5f}_Pk.npy'.format(tracer, x), r.power['power'].real - poles.attrs['shotnoise'])
      
      # The multipoles. 
      for ell in [0]:    
        k      = poles['k']
        P      = poles['power_%d' % ell].real

        if ell == 0:
          P    = P - poles.attrs['shotnoise']

        np.savetxt('dat/pk_{:.5f}.txt'.format(x), np.c_[k, P])            

      break
        
