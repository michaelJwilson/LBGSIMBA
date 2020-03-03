import pickle
import fitsio
import pandas                        as      pd
import numpy                         as      np
import time
import astropy.io.fits               as      fits

from   snaps                         import  snaps

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

def get_simba(tracer, x, space='', insample=0):
    fpath           = '/home/mjwilson/LBGSIMBA/bigdat/simba_{}{}pos_{:.5f}.fits'.format(tracer, space, x)
    cat             = FITSCatalog(fpath)

    print('Retrieved {}.'.format(fpath))
    
    ##  Comoving Mpc/h.
    cat['Position'] = StackColumns(cat['x'], cat['y'], cat['z'])

    if (tracer == 'g') & insample:
      isin  = cat['INSAMPLE'].compute()
        
      print('Retaining {} of {} in sample.'.format(np.count_nonzero(isin), len(cat)))

      cat   = cat[isin.astype(bool)]
      
    return  cat

def get_abacus(bound=720., ArrayCatalog=True):
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

    if ArrayCatalog:
        ##  Comoving Mpc/h.
        return  ArrayCatalog(sample)

    else:
        return  sample

compute    = True

t0         = time.time()

##  Simba. 
cols       = ['x', 'y', 'z']

tracer     =  'g'  ##  ['dm', 'g'] 

insample   =  1

##  Real or redshift space. 
for space in ['']:
  for x in snaps.keys():
    if compute:
      setup_logging()

      ##  Simba.
      dk              =  0.1
      kmin            = 0.01
      boxsize         =  100.
      Nmesh           = 1024

      cat             = get_simba(tracer, x, insample=insample)

      ##  Abacus
      ##  kmin        = 0.01
      ##  boxsize     = 100.  ##  Natural:  720.
      ##  cat         = get_abacus(bound=boxsize)

      print('Solving for redshift:  {}, tracer:  {} and space: {}.'.format(x, tracer, space))
      
      ##  Box size in Mpc/h.
      mesh            = cat.to_mesh(resampler='tsc', Nmesh=Nmesh, compensated=True, BoxSize=boxsize)
      
      r               = FFTPower(mesh, mode='2d', dk=dk, kmin=kmin, Nmu=60, los=[0,0,1], poles=[0,2,4], BoxSize=boxsize, Nmesh=Nmesh)

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

      # np.save('dat/{}tdpk_{}_{:.5f}_k.npy'.format(space, tracer, x),  r.power['k'])
      # np.save('dat/{}tdpk_{}_{:.5f}_mu.npy'.format(space, tracer, x), r.power['mu'])
      # np.save('dat/{}tdpk_{}_{:.5f}_Pk.npy'.format(space, tracer, x), r.power['power'].real - poles.attrs['shotnoise'])
      
      # The multipoles. 
      for ell in [0]:    
        k      = poles['k']
        P      = poles['power_%d' % ell].real

        if ell == 0:
          P    = P - poles.attrs['shotnoise']

        np.savetxt('dat/{}{}pk_{:.5f}_insample_{}.txt'.format(space, tracer, x, insample), np.c_[k, P, poles.attrs['shotnoise'] * np.ones_like(P)])            
      
        

        
