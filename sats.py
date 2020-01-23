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

def sats(test, boxsize, redshift):
    caesar      =  get_caesar(boxsize, redshift)

    pos         =  [x.pos.to('Mpc/h') for x in caesar.galaxies]  # comoving Mpc/h.                                                                                                               
    iscentral   =  np.array([x.central         for x in caesar.galaxies]).astype(np.int)

    ngal        =  len(iscentral)
    ncentral    =  np.sum(iscentral)
    nsat        =  np.sum(1 - iscentral)
    sfrac       =  1. * nsat / ngal
    
    print('\n\nFor redshift of {:.4f}'.format(redshift))
    print('\n\nNumber of galaxies found: {:d}'.format(ngal))
    print('Number of centrals found: {:d}'.format(ncentral))
    print('Number of satellites found: {:d}'.format(nsat))
    print('\nThe satellite fraction is {:.3f}'.format(sfrac))
    print('\n\n')
    
    return  caesar


if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.\n\n')

    test        =   True

    boxsize     =  100.           #  [Mpc/h];  hubble      =  0.68                                                                                                                                                                   
    vol         =  boxsize ** 3.

    for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
      caesar    = sats(test, boxsize, redshift)

    print('\n\nDone.\n\n')
