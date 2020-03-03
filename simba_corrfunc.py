import  glob
import  h5py
import  Corrfunc
import  fitsio
import  numpy                         as      np

from    astropy.table                 import  Table
from    scipy.spatial                 import  KDTree 
from    itertools                     import  product
from    snaps                         import  snaps


def get_simba(tracer, x, space='', insample=0):
    fpath  = '/home/mjwilson/LBGSIMBA/bigdat/simba_{}{}pos_{:.5f}.fits'.format(tracer, space, x)

    print('Retrieving {}.'.format(fpath))

    cat    = fitsio.read(fpath)

    if tracer == 'g':
      isin = cat['INSAMPLE']
      
    ##  Comoving Mpc/h.                                                                                                                                                                   
    cat    = np.c_[cat['x'], cat['y'], cat['z']]

    if tracer == 'dm':
        # sub-sample.
        cat   = cat[::10] 

    if (tracer == 'g') & insample:
        print('Retaining {} of {} in sample.'.format(np.count_nonzero(isin), len(cat)))

        cat   =  cat[isin.astype(bool)]
        
    return  cat
        
def calc_xi(test, boxsize, redshift, tracer, space, insample=0):    
    print('Solving for Test: {}, boxsize: {}, redshift: {}, tracer: {} and space: {}.'.format(test, boxsize, redshift, tracer, space))

    # Comoving Mpc/h. 
    pos         =  get_simba(tracer, redshift, space=space, insample=insample)
    
    if test:
      pos = pos[:1000]
    
    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
                
    ##
    bins        = np.logspace(-1., 1.47, 30)

    results     = Corrfunc.theory.xi(X=pos[:,0], Y=pos[:,1], Z=pos[:,2], boxsize=boxsize, nthreads=32, binfile=bins, output_ravg=True, verbose=True)

    print(results)
    
    ##  Save result:
    rs          = results['ravg']
    
    np.savetxt('dat/corrfuncxi_{}space_{}_{:.3f}_insample_{}.txt'.format(space, tracer, redshift, insample), np.c_[rs, results['xi']], fmt='%.6le')
    
def plot_xi(space='', insample=0):
    import  matplotlib; matplotlib.use('pdf')
    import  pylab                         as      pl
    import  matplotlib.pyplot             as      plt

    from    sutils                        import latexify
    from    mcfit                         import  P2xi

    
    latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

    bs     = [2.200000, 3.200000, 4.200000, 5.5000]
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, redshift in enumerate([2.024621, 3.00307, 3.963392, 5.0244]):
        fpath  = 'dat/corrfuncxi_{}space_{}_{:.3f}_insample_{}.txt'.format(space, 'g', redshift, insample)
                
        rs, xi = np.loadtxt(fpath, unpack=True)
        pl.loglog(rs, xi, lw=0, alpha=0.8, label='{:.2f}'.format(redshift), c=colors[i], marker='^', markersize=3)

        fpath  = 'dat/corrfuncxi_{}space_{}_{:.3f}.txt'.format(space, 'dm', redshift)

        rs, xi = np.loadtxt(fpath, unpack=True)
        pl.loglog(rs, xi, lw=0, alpha=0.8, c=colors[i], marker='s', markersize=3)
                
        # ZA.
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))
        
        # Lagrangian bias parameters: bE_1 = 1 + bL_1 (e.g. eqn. 213 of 1611.09787).
        b1, b2 = 0.0, 0.0
        
        # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2    
        cc     = np.array([0., 1.0, b1, b2, b1**2, b1*b2, b2**2])
        cc     = np.array([0., 1.0, b1, b2, 0.0,   0.0,     0.0])
        
        rs     = _[:,  0]
        _      = _[:,1:8] / rs[:,None] / rs[:,None]
        
        result = np.dot(_, cc)
        
        pl.loglog(rs, result, lw=1, alpha=0.5, linestyle='-', c=colors[i])
        
        # Linear (MCFIT).
        k, P, hf = np.loadtxt('dat/white/pklin_z{}.txt'.format(iz), unpack=True)
        # rs, xi = P2xi(k)(P)
        # pl.loglog(rs, (1. + b1) * (1. + b1) * xi, lw=1, alpha=0.8, linestyle='--', c=colors[i])

        # Halofit (MCFIT).                                                                                                                                                                 
        rs, xi   = P2xi(k)(hf)
        pl.loglog(rs, (1. + b1) * (1. + b1) * xi, lw=1, alpha=0.8, linestyle=':', c=colors[i])

        # Halofit (MCFIT).
        pl.loglog(rs, bs[i] * bs[i] * xi, lw=1, alpha=0.8, linestyle='-', c=colors[i])
        
    ##
    pl.xlim(0.5,  3.e1)
    pl.ylim(0.01, 2.e1)

    pl.xlabel(r'$r \ [h^{-1} \rm{Mpc}]$')
    pl.ylabel(r'$\xi(r)$')

    pl.legend(loc=3, frameon=False)
    
    plt.tight_layout()

    pl.savefig('plots/simba_corrfunc_xi_insample_{}.pdf'.format(insample))
    
    
if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    insample    =      1
    test        =  False

    #  [Mpc/h];  hubble      =  0.68  
    boxsize     =  100.     
    vol         =  boxsize ** 3.

    redshifts   = [2.024621, 3.00307, 3.963392, 5.0244]
    
    tracer      = 'g'   ##  ['g', 'dm']
    space       = ''    ##  ['z', '']
    '''
    for redshift in redshifts:          
      calc_xi(test, boxsize, redshift, tracer, space, insample=insample)
    '''
    plot_xi(space='', insample=1)
      
    print('\n\nDone.\n\n')
