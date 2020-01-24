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
    caesar      =  get_caesar(boxsize, redshift, load_halo=True)

    pos         =  [x.pos.to('Mpc/h') for x in caesar.galaxies]  # comoving Mpc/h.                                                                                                               
    iscentral   =  np.array([x.central for x in caesar.galaxies]).astype(np.int)
    issatellite =  1 - iscentral
    
    ngal        =  len(iscentral)
    ncentral    =  np.sum(iscentral)
    nsat        =  np.sum(issatellite)
    sfrac       =  1. * nsat / ngal

    print('\n\nFor redshift of {:.4f}'.format(redshift))
    print('\n\nNumber of galaxies found: {:d}'.format(ngal))
    print('Number of centrals found: {:d}'.format(ncentral))
    print('Number of satellites found: {:d}'.format(nsat))
    print('\nThe satellite fraction is {:.3f}'.format(sfrac))
    print('\n\n')
    
    ##  Relative galaxy position from parent halo position (COM or?)
    vecs        =  [x.pos.to('Mpccm/h') - x.halo.pos.to('Mpccm/h') for x in caesar.galaxies]

    ##  vecs, but for only the satellites
    svecs       =  [x for x, y in zip(vecs, issatellite) if y == True]

    ##  Sanity check
    assert  len(svecs) == nsat

    ##  Virial radius for all the parent halos. 
    rs          =  [x.halo.radii['virial'].to('Mpccm/h') for x in caesar.galaxies]
    srs         =  [x for x, y in zip(rs, issatellite) if y == True]

    ##  Relative satellite position from parent halo position (COM or?) in units of
    ##  the halo virial radius.
    _           =  [] 
    
    for i, s in enumerate(svecs):
        _.append(s / srs[i])        
        
    svecs       =  _

    ss          = np.array([np.sqrt(np.sum(x**2.)).value for x in svecs]) 
    ss          = np.sort(ss)

    bins        = np.arange(0., 1.1, 0.05)
    cnts, _     = np.histogram(ss, bins=bins)

    pl.loglog(bins[:-1], cnts, 'k-', label='{:.2f}'.format(redshift))
    
    return  caesar


if __name__ == '__main__':
    print('\n\nWelcome to Simba xi.\n\n')

    test        =   True

    boxsize     =  100.           #  [Mpc/h];  hubble      =  0.68                                                                                                                                                                   
    vol         =  boxsize ** 3.

    for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
      caesar    = sats(test, boxsize, redshift)

      break
    
    pl.xlabel(r'$r / r_{\rm{vir}}$')
  
    pl.legend(frameon=False, loc=3)
  
    pl.savefig('plots/radial.pdf'.format(redshift))
    
    print('\n\nDone.\n\n')
