import pickle
import pandas                        as      pd
import numpy                         as      np
import time
import astropy.io.fits               as      fits

from   get_data                      import  snaps

from   nbodykit.lab                  import  *
from   nbodykit.io.fits              import  FITSFile
from   nbodykit.io.csv               import  CSVFile
from   nbodykit.source.catalog.file  import  CSVCatalog, FITSCatalog
from   nbodykit.source.catalog.array import  ArrayCatalog
from   nbodykit.transform            import  StackColumns
from   nbodykit                      import  setup_logging
from   AbacusCosmos                  import  Halos


##  source activate nbodykit-env                                                                                                                                                                                                   
##                                                                                                                                                                                                                                   
##  export PYTHONPATH=/home/mjwilson/LBGSIMBA/:/home/mjwilson/LBGSIMBA/AbacusCosmos/:$PYTHONPATH                                                                                                                                     
##                                                                                                                                                                                                                                  
##  python3 nkit.py   

def get_simba(tracer, x, space=''):
    cat             = FITSCatalog('/home/mjwilson/LBGSIMBA/bigdat/simba_{}{}pos_{:.5f}.fits'.format(tracer, space, x))

    ##  Comoving Mpc/h.
    cat['Position'] = StackColumns(cat['x'], cat['y'], cat['z'])

    return  cat
    

compute    = True

t0         = time.time()

##  Simba. 
cols       = ['x', 'y', 'z']

Nmesh      = 1024

##  Real or redshift space, ['z', ''].
for space in ['']:
  for x in snaps.keys():
    if compute:
      setup_logging()

      ##  Simba.
      kmin            = 0.01
      boxsize         =  100.

      gcat            = get_simba( 'g', x)
      dmcat           = get_simba('dm', x)
      
      print('Solving for redshift:  {}, cross-spectra  and space: {}.'.format(x, space))
      
      ##  Box size in Mpc/h.
      gmesh           =  gcat.to_mesh(resampler='tsc', Nmesh=Nmesh, compensated=True, BoxSize=boxsize)
      dmmesh          = dmcat.to_mesh(resampler='tsc', Nmesh=Nmesh, compensated=True, BoxSize=boxsize)
      
      r               = FFTPower(first=gmesh, mode='2d', dk=0.1, kmin=0.01, Nmu=60, los=[0,0,1], poles=[0,2,4], BoxSize=boxsize, Nmesh=Nmesh, second=dmmesh)

      poles           = r.poles

      t1              = time.time()

      print(t1 - t0)

      print('Number of galaxies: {}'.format(len(gcat['Position'])))
      print('Shotnoise: {}'.format(poles.attrs['shotnoise']))

      for _ in r.power.attrs:
          print("%s = %s" %(_, str(r.power.attrs[_])))

      # The 2D power spectrum results.
      print("power = ", r.power)
      print("variables = ", r.power.variables)

      for name in r.power.variables:
          var = r.power[name]

          print("'%s' has shape %s and dtype %s" % (name, var.shape, var.dtype))
      
      # The multipoles. 
      for ell in [0]:    
        k      = poles['k']

        ##  https://arxiv.org/pdf/1309.5556.pdf:  eqn. A41.
        P      = 0.5 * (poles['power_%d' % ell] + np.conjugate(poles['power_%d' % ell]))
        P      = P.real
        
        np.savetxt('dat/{}{}pk_{:.5f}.txt'.format(space, 'c', x), np.c_[k, P, poles.attrs['shotnoise'] * np.ones_like(P)])            
        
