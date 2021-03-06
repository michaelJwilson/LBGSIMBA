import  matplotlib; matplotlib.use('pdf')

import  sys
import  glob
import  h5py
import  Corrfunc 
import  fitsio
import  numpy               as      np
import  pylab               as      pl
import  matplotlib.pyplot   as      plt

from    scipy.spatial       import  KDTree 
from    itertools           import  product
from    hod                 import  get_data
from    sutils              import  latexify
from    get_data            import  get_caesar, snaps
from    Corrfunc.theory.DD  import  DD
from    Corrfunc.utils      import  convert_3d_counts_to_cf
from    simba_corrfunc      import  get_simba


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def get_abacus(bound=720.):
    from   AbacusCosmos      import Halos

    
    fpath           = '/disk01/mjwilson/abacus/AbacusCosmos_720box_planck_products/AbacusCosmos_720box_planck_00-0_products/AbacusCosmos_720box_planck_00-0_FoF_halos/z1.500'

    # The path to the catalog will be different on your system                                                                                                                                                                     
    cat             = Halos.make_catalog_from_dir(dirname=fpath, load_subsamples=False, load_pids=False, halo_type='FoF')

    halos           = cat.halos
    
    ##  [-360., 360.] to [0., 720.]
    halos['pos']    = halos['pos'] + 360.

    isin            = (halos['pos'][:,0] <= bound) & (halos['pos'][:,1] <= bound) & (halos['pos'][:,2] <= bound)
    sample          =  halos['pos'][isin]

    assert  sample.max() <= bound

    return  sample

def calc_xi(test, boxsize, redshift, tracer, space, insample=0):    
    # comoving Mpc/h.                                                                                                                                                                    
    pos         =  get_simba(tracer, redshift, space, insample=insample).astype(np.float32)
    
    if test:
      pos = pos[:10000]

    else:
      pos = pos[:5000000]
      
    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
    
    ##
    bins        = np.arange(0., 26., 1.)
    rs          = (bins[:-1] + bins[1:]) / 2.

    ##  Randoms                                                                                                                                                                                              
    nrand       = 3 * len(pos)

    np.random.seed(42)

    print('Picking {} and {} galaxies and randoms respectively.'.format(len(pos), nrand))
    
    randx       = np.random.uniform(0, boxsize, nrand).astype(np.float32)
    randy       = np.random.uniform(0, boxsize, nrand).astype(np.float32)
    randz       = np.random.uniform(0, boxsize, nrand).astype(np.float32)

    rpos        = np.c_[randx, randy, randz]
    
    # Distances along the :math:\pi direction are binned with unit depth. For instance, if pimax=40, then 40 bins will be created along the pi direction.
    DD          = Corrfunc.theory.DDrppi(X1= pos[:,0], Y1= pos[:,1], Z1= pos[:,2], periodic=True, boxsize=boxsize, nthreads=4, binfile=bins, autocorr=True,  pimax=25.0, output_rpavg=True)
    RR          = Corrfunc.theory.DDrppi(X1=rpos[:,0], Y1=rpos[:,1], Z1=rpos[:,2], periodic=True, boxsize=boxsize, nthreads=4, binfile=bins, autocorr=True,  pimax=25.0, output_rpavg=True)
    DR          = Corrfunc.theory.DDrppi(X1= pos[:,0], Y1= pos[:,1], Z1= pos[:,2],\
                                         X2=rpos[:,0], Y2=rpos[:,0], Z2=rpos[:,0],\
                                         pimax=25., periodic=True, boxsize=boxsize,\
                                         nthreads=4, binfile=bins, autocorr=False,\
                                         output_rpavg=True)
    
    xi          = convert_3d_counts_to_cf(len(pos), len(pos), nrand, nrand, DD, DR, DR, RR, estimator=u'LS')

    rps         = np.array([x['rpavg'] for x in RR])
    pis         = np.array([x['pimax'] for x in RR]) - 0.5
    
    ##  Save result:
    np.savetxt('dat/corrfunc3dxi_{}_{}_{:.3f}_insample_{}.txt'.format(tracer, space, redshift, insample), np.c_[rps, pis, xi], fmt='%.6le')

    ##  Abacus
    ##  np.savetxt('dat/corrfunc3dxi_abacus.txt', np.c_[rps, pis, xi], fmt='%.6le')   
        
if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    test        =  False

    insample    =  1
    
    #  [Mpc/h];  hubble      =  0.68  
    boxsize     =  100.     
    vol         =  boxsize ** 3.

    ##  5.0244
    redshifts   = [2.024621, 3.00307, 3.963392]

    tracer      = 'g'
    space       = ''
    
    for space in ['', 'z']:
      for redshift in redshifts:
        calc_xi(test, boxsize, redshift, tracer, space, insample=insample)
          
    print('\n\nDone.\n\n')
