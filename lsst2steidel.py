import matplotlib;  matplotlib.use('PDF')

import numpy    as     np
import pylab    as     pl

from   insample import insample_steidel
from   insample	import insample

# Test of relations of App. A of LBGCMB.

def lsst2steidel():
    import pylab    as     pl
    from   get_data import get_pyloser


    nrows     =  -1
    boxsize   = 100.

    ##  2.024621,  3.963392                                                                                                                                                                                                          
    wave, two, ids = get_pyloser(boxsize, 2.024621, nrows=nrows, steidel=True)

    bxdrops        = insample_steidel('BX', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values)
    bmdrops        = insample_steidel('BM', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values)

    udrops         = insample('u',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)
    gdrops         = insample('g',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)
    rdrops         = insample('r',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)

    print(np.count_nonzero(bxdrops), np.count_nonzero(bmdrops), np.count_nonzero(udrops), np.count_nonzero(gdrops), np.count_nonzero(rdrops))

    return  two, bxdrops, bmdrops


if __name__ == '__main__':
    two, bxdrops, bmdrops = lsst2steidel()

    two['U-G']            = two['steidel_un']     - two['steidel_g']
    two['*U-G*']          = 0.97 * (two['LSST_u'] - two['LSST_g']) + 1.27
    
    two['G-R']            = two['steidel_g']      - two['steidel_rs']
    two['*G-R*']          = 0.32 * (two['LSST_r'] - two['LSST_i']) + 1.10
    
    pl.plot(two['U-G'][bmdrops], two['*U-G*'][bmdrops], marker='.', c='k', lw=0.0, markersize=2)

    pl.savefig('plots/lsst2steidel.pdf')
    
    print('\n\nDone.\n\n')


