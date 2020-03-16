import matplotlib;  matplotlib.use('PDF')

import numpy    as     np
import pylab    as     pl

from   insample import insample_steidel
from   insample	import insample

# Test of relations of App. A of LBGCMB.

def lsst2steidel():
    import pylab    as     pl
    from   get_data import get_pyloser


    nrows          =   -1
    boxsize        =  100.

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

    two['U-G']            = two['steidel_un'] - two['steidel_g']
    two['u-g']            = two['LSST_u']     - two['LSST_g'] 
    
    two['G-R']            = two['steidel_g']  - two['steidel_rs']
    two['g-r']            = two['LSST_g']     - two['LSST_r']
    
    pl.plot(two['U-G'],          two['u-g'],          marker='.', c='k',    lw=0.0, markersize=1, alpha=0.2)
    pl.plot(two['U-G'][bxdrops], two['u-g'][bxdrops], marker='.', c='gold', lw=0.0, markersize=1, alpha=0.2)
    pl.plot(two['U-G'][bmdrops], two['u-g'][bmdrops], marker='.', c='c',    lw=0.0, markersize=1, alpha=0.2)

    pl.xlabel(r'$U_n-G$')
    pl.ylabel(r'$u-g$')
    
    pl.savefig('plots/lsst2steidel_umg.pdf')

    ##
    pl.clf()
    
    pl.plot(two['G-R'],          two['g-r'],          marker='.', c='k',    lw=0.0, markersize=1, alpha=0.2)
    pl.plot(two['G-R'][bxdrops], two['g-r'][bxdrops], marker='.', c='gold', lw=0.0, markersize=1, alpha=0.2)
    pl.plot(two['G-R'][bmdrops], two['g-r'][bmdrops], marker='.', c='c',    lw=0.0, markersize=1, alpha=0.2)

    pl.xlabel(r'$G-R$')
    pl.ylabel(r'$g-r$')

    pl.savefig('plots/lsst2steidel_gmr.pdf')

    ##
    pl.clf()

    pl.plot(two['g-r'][bxdrops], two['u-g'][bxdrops], marker='.', c='gold', lw=0.0, markersize=1, alpha=0.2, label='BX')
    pl.plot(two['g-r'][bmdrops], two['u-g'][bmdrops], marker='.', c='c', lw=0.0, markersize=1, alpha=0.2, label='BM')

    pl.xlabel(r'$g-r$')
    pl.ylabel(r'$u-g$')
    pl.legend(frameon=False)
    
    pl.savefig('plots/lsst2steidel.pdf')

    ##                                                                                                                                                                                                                               
    pl.clf()

    pl.plot(two['G-R'][bxdrops], two['U-G'][bxdrops], marker='.', c='gold', lw=0.0, markersize=1, alpha=0.2, label='BX')
    pl.plot(two['G-R'][bmdrops], two['U-G'][bmdrops], marker='.', c='c', lw=0.0, markersize=1, alpha=0.2, label='BM')

    pl.xlabel(r'$G-R$')
    pl.ylabel(r'$Un-G$')
    
    pl.legend(frameon=False)
    
    pl.savefig('plots/lsst2steidel_steidel.pdf')

    print('\n\nDone.\n\n')


