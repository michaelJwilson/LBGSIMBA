import numpy as np
import pylab as pl


def insample(selection='u', u, g, r, i, z, y):
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
    
    else:
        raise  ValueWarning('Specific selection ({}) is not available.'.format(selection))


if __name__ == '__main__':
    from  get_data  import  get_pyloser

    print('\n\nDone.\n\n')
