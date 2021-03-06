import  matplotlib; matplotlib.use('pdf')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree 
from    itertools         import  product
from    hod               import  get_data
from    utils             import  latexify
from    get_data          import  get_caesar, snaps


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def calc_xi(test, reflect, boxsize, redshift):    
    caesar      =  get_caesar(boxsize, redshift)

    pos         =  [x.pos.to('Mpc/h') for x in caesar.galaxies]  # comoving Mpc/h. 
    iscentral   =  np.array([x.central for x in caesar.galaxies]).astype(np.int)
    
    ##  Basic stats.
    print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
    print('Number of centrals found: {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))
            
    if test:
      pos = pos[:1000]

    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
    
    catalogue   = np.copy(pos)
    
    if reflect:
      ##  Apply periodic reflections. 
      reflections = list(product('+0-', repeat=3))
      reflections.remove(tuple(['0', '0', '0']))

      convert     = {'+': +1., '0': 0., '-': -1.}

      for reflection in reflections:
        _copy       = np.copy(pos)
        signs       = [convert[x] for x in reflection]
    
        _copy[:,0] += signs[0] * boxsize
        _copy[:,1] += signs[1] * boxsize
        _copy[:,2] += signs[2] * boxsize

        catalogue   = np.vstack([catalogue, _copy])
        
      print(len(pos), 27 * len(pos), len(catalogue))
      print(catalogue.shape)
    
    ##  Grow tree. 
    PTree       = KDTree(pos)
    CTree       = KDTree(catalogue)
    
    ##  Note:  asymmetric catalogue-based Tree call and pos call - i.e. do not count pairs between
    ##         two reflections. 
    paired      =  CTree.query_ball_tree(PTree, 100.)
        
    dr          =  0.25
    bins        =  np.arange(-1.0, 100.0, dr)
        
    sep         =  []

    for i, row in enumerate(catalogue):
        for j, twin in enumerate(paired[i]):
            sep.append(np.sum((row - catalogue[twin])**2.))

    sep         =  np.sqrt(np.array(sep))
    bsep        =  np.digitize(sep, bins=bins)

    meanr       =  np.array([np.mean(sep[bsep == x]) for x in range(len(bins) - 1)]) 
    cnts        =  np.array([np.sum(bsep == x) for x in range(len(bins) - 1)])
                
    rr          =  (4. * np.pi / 3.) * nbar * nbar * vol * (( bins[:-1] + dr / 2. ) ** 3. - ( bins[:-1] - dr / 2. ) ** 3)
    xi          =  cnts / rr - 1.

    ##  Save result:
    np.savetxt('dat/xi_{:.3f}.txt'.format(redshift), np.c_[bins[:-1] + dr / 2., meanr, xi], fmt='%.6le')
    
def plot_xi():
    # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2

    for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
        midr, meanr, xi = np.loadtxt('dat/xi_{:.3f}.txt'.format(redshift), unpack=True)

        pl.loglog(meanr, xi, lw=1, alpha=0.8, label='{:.2f}'.format(redshift))

        # ZA
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

        b1, b2 = 1.0, 0.0

        cc     = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])

        rs     = _[:,0]
        _      = _[:,0:8]
        result = np.dot(_, cc)
        
        pl.loglog(rs, result, lw=1, alpha=0.8, linestyle='--')

        
    #pl.xlim(0.1,  60.5)
    #pl.ylim(0.1,  2.e3)

    pl.xlabel(r'$r \ [h^{-1} \rm{Mpc}]$')
    pl.ylabel(r'$\xi(r)$')

    pl.legend(loc=3, frameon=False)
    
    plt.tight_layout()

    pl.savefig('plots/xi.pdf')
    
    
if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    test        =  False
    reflect     =   True

    boxsize     =    100.           #  [Mpc/h];  hubble      =  0.68                                                                                                                                         
    vol         =  boxsize ** 3.

    redshifts   = [2.024621, 3.00307, 3.963392, 5.0244]
    redshifts   = [3.963392, 5.0244]
    
    for redshift in redshifts:
      calc_xi(test, reflect, boxsize, redshift)

    plot_xi()
      
    print('\n\nDone.\n\n')
