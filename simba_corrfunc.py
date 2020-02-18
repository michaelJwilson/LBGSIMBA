import  matplotlib; matplotlib.use('pdf')

import  glob
import  h5py
import  Corrfunc
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree 
from    itertools         import  product
from    hod               import  get_data
from    utils             import  latexify
from    get_data          import  get_caesar, snaps
from    mcfit             import  P2xi


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def calc_xi(test, boxsize, redshift):    
    caesar      =  get_caesar(boxsize, redshift)

    # comoving Mpc/h. 
    pos         =  np.array([list(x.pos.to('Mpc/h')) for x in caesar.galaxies])
    iscentral   =  np.array([x.central               for x in caesar.galaxies]).astype(np.int)
    
    ##  Basic stats.
    print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))

    print('Number of centrals found: {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))
            
    if test:
      pos = pos[:100]

    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
                
    ##
    bins        = np.logspace(0., 1.47, 25)
    rs          = (bins[:-1] + bins[1:]) / 2.
    
    results     = Corrfunc.theory.xi(X=pos[:,0], Y=pos[:,1], Z=pos[:,2], boxsize=boxsize, nthreads=4, binfile=bins, periodic=True)

    ##  Save result:
    np.savetxt('dat/corrfuncxi_{:.3f}.txt'.format(redshift), np.c_[rs, results['xi']], fmt='%.6le')
    
def plot_xi():
    colors  = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, redshift in enumerate([2.024621, 3.00307, 3.963392, 5.0244]):
        rs, xi = np.loadtxt('dat/corrfuncxi_{:.3f}.txt'.format(redshift), unpack=True)
        pl.loglog(rs, xi, lw=0, alpha=0.8, label='{:.2f}'.format(redshift), c=colors[i], marker='^', markersize=3)
        
        # ZA.
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))
        
        # Lagrangian bias parameters: bE_1 = 1 + bL_1 (e.g. eqn. 213 of 1611.09787).
        b1, b2 = 10.0, 0.0
        '''
        # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2    
        cc     = np.array([0., 1.0, b1, b2, b1**2, b1*b2, b2**2])
        cc     = np.array([0., 1.0, b1, b2, 0.0,   0.0,     0.0])
        
        rs     = _[:,  0]
        _      = _[:,1:8] / rs[:,None] / rs[:,None]
        
        result = np.dot(_, cc)
        
        pl.loglog(rs, result, lw=1, alpha=0.8, linestyle='-', c=colors[i])
        '''
        # Linear (MCFIT).
        k, P, hf = np.loadtxt('dat/white/pklin_z{}.txt'.format(iz), unpack=True)
        rs, xi   = P2xi(k)(P)

        pl.loglog(rs, (1. + b1) * (1. + b1) * xi, lw=1, alpha=0.8, linestyle='--', c=colors[i])
        
    pl.xlim(1.e-3,  5.e3)
    pl.ylim(0.01,   5.e2)

    pl.xlabel(r'$r \ [h^{-1} \rm{Mpc}]$')
    pl.ylabel(r'$\xi(r)$')

    pl.legend(loc=3, frameon=False)
    
    plt.tight_layout()

    pl.savefig('plots/simba_corrfunc_xi.pdf')
    
    
if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    test        =  False

    #  [Mpc/h];  hubble      =  0.68  
    boxsize     =  100.     
    vol         =  boxsize ** 3.

    redshifts   = [2.024621, 3.00307, 3.963392, 5.0244]
    
    ##  for redshift in redshifts:
    ##    calc_xi(test, boxsize, redshift)

    plot_xi()
      
    print('\n\nDone.\n\n')
