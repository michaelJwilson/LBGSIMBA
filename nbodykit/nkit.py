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

def get_simba(tracer, x):
    cat             = FITSCatalog('/home/mjwilson/LBGSIMBA/bigdat/simba_{}zpos_{:.5f}.fits'.format(tracer, x))

    ##  Comoving Mpc/h.
    cat['Position'] = StackColumns(cat['x'], cat['y'], cat['z'])

    return  cat
    
def get_abacus(bound=720.):
    fpath           = '/disk01/mjwilson/abacus/AbacusCosmos_720box_planck_products/AbacusCosmos_720box_planck_00-0_products/AbacusCosmos_720box_planck_00-0_FoF_halos/z1.500'

    # The path to the catalog will be different on your system                                                                                                                                                                      
    cat             = Halos.make_catalog_from_dir(dirname=fpath, load_subsamples=False, load_pids=False, halo_type='FoF')

    halos           = cat.halos

    ##  [-360., 360.] to [0., 720.]                                                                                                                                                                                                \
    halos['pos']   += 360.

    isin            = (halos['pos'][:,0] <= bound) & (halos['pos'][:,1] <= bound) & (halos['pos'][:,2] <= bound)
    sample          =  halos['pos'][isin]
    
    assert  sample.max() <= bound
    
    sample          = {'Position':  sample}
    
    ##  Comoving Mpc/h.
    return  ArrayCatalog(sample)


compute    = True

t0         = time.time()

##  Simba. 
cols       = ['x', 'y', 'z']
tracer     =  'g' # ['dm', 'g']

for x in snaps.keys():
    if compute:
      setup_logging()

      ##  Simba.
      kmin            = 0.01
      boxsize         = 100.
      cat             = get_simba(tracer, x)

      ##  Abacus
      ##  kmin        = 0.01
      ##  boxsize     = 100.  ##  Natural:  720.
      ##  cat         = get_abacus(bound=boxsize)
      
      ##  Box size in Mpc/h.
      mesh            = cat.to_mesh(resampler='tsc', Nmesh=1024, compensated=True, BoxSize=boxsize)
      
      r               = FFTPower(mesh, mode='2d', dk=0.05, kmin=0.01, Nmu=60, los=[0,0,1], poles=[0,2,4], BoxSize=boxsize, Nmesh=1024)

      poles           = r.poles

      t1              = time.time()

      print(t1 - t0)

      print('Number of galaxies: {}'.format(len(cat['Position'])))
      print('Shotnoise: {}'.format(poles.attrs['shotnoise']))

      for _ in r.power.attrs:
          print("%s = %s" %(_, str(r.power.attrs[_])))

      # The 2D power spectrum results
      print("power = ", r.power)
      print("variables = ", r.power.variables)

      for name in r.power.variables:
          var = r.power[name]

          print("'%s' has shape %s and dtype %s" % (name, var.shape, var.dtype))

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
        # np.savetxt('dat/abacus_pk_100.txt', np.c_[k, P]) 
        

        
