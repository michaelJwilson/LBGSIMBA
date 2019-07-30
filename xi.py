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


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.')

    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
    boxsize     = 100.
    getredshift = 3.00307

    f, p        =  get_data(boxsize, getredshift)

    ##
    gid         =  f['galaxy_data']['GroupID'][:]
    iscentral   =  f['galaxy_data']['central'][:]
    haloindex   =  f['galaxy_data']['parent_halo_index'][:]

    ##  print(p['COLOR_INFO'][:])

    ##  Basic stats.
    print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
    print('Number of centrals found: {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))
        
    ##  Positions in kpc.
    pos         = f['galaxy_data']['pos'][:]
    pos        /= 1.e3

    ngal        =  len(pos)
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol
    
    ##  Test.
    pos         = pos[:50]
    
    ##  Apply periodic reflections.
    catalogue   = np.copy(pos)

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

    print(catalogue)

    ##  Grow tree. 
    PTree       =  KDTree(pos)
    CTree       =  KDTree(catalogue)

    if compute:
        ##  Note:  asymmetric catalogue-based Tree call and pos call - i.e. do not count pairs between
        ##         two reflections. 
        paired      =  CTree.query_ball_tree(PTree, 25.)

        dr          =  0.1
        bins        =  np.arange(0.0, 12.5, dr)
        
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
        np.savetxt('dat/xi.txt', np.c_[meanr, xi], fmt='%.6le')
        
    else:
        meanr, xi   = np.loadtxt('dat/xi.txt')
        
    ##  rint(bins + dr / 2.) 
    print(meanr)
    ##  print(cnts)
    ##  print(rr)
    print(xi)

    pl.semilogx(meanr, xi, lw=1, c='darkcyan', alpha=0.8)

    pl.xlim(0.1, 12.5)
    pl.ylim(0.01, 10.)
    
    pl.xlabel(r'$r \ [h^{-1} \rm{Mpc}]$')
    pl.ylabel(r'$\xi(r)$')

    plt.tight_layout()
    
    pl.savefig('plots/xi.pdf')

    ##  Save result:
    np.savetxt('dat/xi.txt', np.c_[meanr, 1. + xi], fmt='%.6le')
    
    print('\n\nDone.\n\n')
