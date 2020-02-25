import matplotlib;  matplotlib.use('PDF')

import numpy    as     np
import pylab    as     pl

from   get_data import get_pyloser


def insample(selection, u, g, r, i, z, y):
    if selection == 'u':
        isin = (u-g) > 1.5
        isin = isin & ((g-r) > -1.0) & ((g-r) < 1.2)
        isin = isin & (1.5 * (g-r) < (u-g) -0.75)
        
        return  isin

    elif selection == 'g':
        isin = (g-r) >	1.0 
        isin = isin & ((g-r) > 1.0)
        isin = isin & ((r-i) <	1.0)
        isin = isin & ((g-r) >	1.5 * (r-i) + 0.8)
        
        return isin

    elif selection == 'r':
        isin = (i-z) >  1.5
        isin = isin & ((z-y) < 0.5)
        isin = isin & ((i-z) >  2.0 * (z-y) + 1.1)

        return isin

    elif selection in ['BX', 'BM']:
        isin = (g-r) >= 0.2
        isin = isin & ((u-g) >= (g-r) + 0.2)
        isin = isin & ((u-g) <  (g-r) + 1.0)
        isin = isin & ((g-r) <= 0.2 * (u-g) + 0.4)
        
    else:
        raise  ValueWarning('Specific selection ({}) is not available.'.format(selection))


if __name__ == '__main__':
    import pylab as pl

    
    nrows     =  -1
    boxsize   = 100.

    ##  2.024621,  3.963392
    wave, two = get_pyloser(boxsize, 2.024621, nrows=nrows)

    bdrops    = insample('BX', two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)
    udrops    = insample('u',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)
    gdrops    = insample('g',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values) 
    rdrops    = insample('r',  two['LSST_u'].values, two['LSST_g'].values, two['LSST_r'].values, two['LSST_i'].values, two['LSST_z'].values, two['LSST_y'].values)

    print(len(two['LSST_u'].values), np.count_nonzero(bdrops), np.count_nonzero(udrops), np.count_nonzero(gdrops), np.count_nonzero(rdrops))
    
    umg       = two['LSST_u'].values - two['LSST_g'].values
    gmr       = two['LSST_g'].values - two['LSST_r'].values
    rmi       = two['LSST_r'].values - two['LSST_i'].values
    
    pl.plot(gmr,         umg,         marker='.', alpha=0.1, lw=0, c='k')
    pl.plot(gmr[bdrops], umg[bdrops], marker='.', lw=0, c='g')

    ##  pl.plot(rmi,         gmr,         marker='.', alpha=0.1, lw=0, c='k')
    ##  pl.plot(rmi[gdrops], gmr[gdrops], marker='.', lw=0, c='g') 

    pl.xlim(-2.0, 2.0)
    pl.ylim(-2.0, 3.5)
    
    pl.savefig('plots/insample.pdf')
    
    print('\n\nDone.\n\n')