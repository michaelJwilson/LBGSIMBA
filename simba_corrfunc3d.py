import  matplotlib; matplotlib.use('pdf')

import  glob
import  h5py
import  Corrfunc 
import  numpy               as      np
import  pylab               as      pl
import  matplotlib.pyplot   as      plt

from    scipy.spatial       import  KDTree 
from    itertools           import  product
from    hod                 import  get_data
from    utils               import  latexify
from    get_data            import  get_caesar, snaps
from    Corrfunc.theory.DD  import  DD
from    Corrfunc.utils      import  convert_3d_counts_to_cf


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

def calc_xi(test, boxsize, redshift):    
    caesar      =  get_caesar(boxsize, redshift)

    # comoving Mpc/h. 
    pos         =  np.array([list(x.pos.to('Mpc/h')) for x in caesar.galaxies])
    iscentral   =  np.array([x.central               for x in caesar.galaxies]).astype(np.int)
       
    pos         =  get_abacus()
    
    if test:
      pos = pos[:10000]
      
    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
                
    ##
    bins        = np.arange(0.1, 25., 0.5)
    rs          = (bins[:-1] + bins[1:]) / 2.

    ##  Randoms                                                                                                                                                                                              
    nrand       = 3 * len(pos)

    np.random.seed(42)

    randx       = np.random.uniform(0, boxsize, nrand)
    randy       = np.random.uniform(0, boxsize, nrand)
    randz       = np.random.uniform(0, boxsize, nrand)

    rpos        = np.c_[randx, randy, randz]

    print('Solving for redshift: {} {}'.format(len(pos), nrand))
    
    # Distances along the :math:\pi direction are binned with unit depth. For instance, if pimax=40, then 40 bins will be created along the pi direction.
    DD          = Corrfunc.theory.DDrppi(X1=pos[:,0],  Y1=pos[:,1],  Z1=pos[:,2], periodic=True, boxsize=boxsize, nthreads=4, binfile=bins, autocorr=True,  pimax=25.0, output_rpavg=True)
    RR          = Corrfunc.theory.DDrppi(X1=rpos[:,0], Y1=rpos[:,1], Z1=rpos[:,2], periodic=True, boxsize=boxsize, nthreads=4, binfile=bins, autocorr=True,  pimax=25.0, output_rpavg=True)
    DR          = Corrfunc.theory.DDrppi(X1=pos[:,0], Y1=pos[:,1], Z1=pos[:,2], X2=rpos[:,0], Y2=rpos[:,0], Z2=rpos[:,0], pimax=25., periodic=True,\
                                         boxsize=boxsize, nthreads=4, binfile=bins, autocorr=False, output_rpavg=True)

    xi          = convert_3d_counts_to_cf(len(pos), len(pos), nrand, nrand, DD, DR, DR, RR, estimator=u'LS')

    rps         = np.array([x['rpavg'] for x in RR])
    pis         = np.array([x['pimax'] for x in RR]) - 0.5
    
    ##  Save result:
    ##  np.savetxt('dat/corrfunc3dxi_{:.3f}.txt'.format(redshift), np.c_[rps, pis, xi], fmt='%.6le')

    ##  Abacus
    np.savetxt('dat/corrfunc3dxi_abacus.txt', np.c_[rps, pis, xi], fmt='%.6le')   
    
def plot_xi():
    colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, redshift in enumerate([2.024621, 3.00307, 3.963392, 5.0244]):
        rs, xi = np.loadtxt('dat/corrfuncxi_{:.3f}.txt'.format(redshift), unpack=True)

        pl.loglog(rs, xi, lw=0, alpha=0.8, label='{:.2f}'.format(redshift), c=colors[i], marker='^', markersize=3)

        # ZA
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

        b1, b2 = 20.0, 0.0

        # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2    
        cc     = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])
        cc     = np.array([0., 0., 1.0, b1, b2, 0.0,   0.0,     0.0])
        
        rs     = _[:,0]
        _      = _[:,0:8]
        result = np.dot(_, cc)
        
        pl.loglog(rs, result / rs / rs, lw=1, alpha=0.8, linestyle='-', c=colors[i])
        
    pl.xlim(0.5,  10.**1.47)
    pl.ylim(0.01, 5.e2)

    pl.xlabel(r'$r \ [h^{-1} \rm{Mpc}]$')
    pl.ylabel(r'$\xi(r)$')

    pl.legend(loc=3, frameon=False)
    
    plt.tight_layout()

    pl.savefig('plots/simba_corrfunc_xi.pdf')
    
    
if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    test        = False

    #  [Mpc/h];  hubble      =  0.68  
    boxsize     =  100.     
    vol         =  boxsize ** 3.

    redshifts   = [2.024621, 3.00307, 3.963392, 5.0244]
    
    for redshift in redshifts:
      calc_xi(test, boxsize, redshift)

      break
      
    ##  plot_3dxi()
      
    print('\n\nDone.\n\n')
