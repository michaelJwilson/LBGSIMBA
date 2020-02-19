import matplotlib; matplotlib.use('PDF')

import numpy  as      np
import pylab  as      pl

from   mcfit             import  P2xi, TophatVar
from   scipy.interpolate import  CubicSpline, interp1d


def get_linearxi(redshift, sigma_8=False):
    iz        = int(100 * redshift + 0.001)
    
    k, P, hf  = np.loadtxt('dat/white/pklin_z{}.txt'.format(iz), unpack=True)    
    r, xi     = P2xi(k)(P)

    if sigma_8:
        R, var    = TophatVar(k)(P, extrap=True)
        varR      = CubicSpline(R, var)
        sigma8    = np.sqrt(varR(8))

    return  interp1d(r, xi, kind='cubic', copy=True, assume_sorted=False, bounds_error=False, fill_value=0.0)

def plot_linearxi():
    pl.loglog(r, xi, label=r'$\sigma_8: {:.3f}$'.format(sigma8))
    pl.legend(loc=2, frameon=False)
    pl.savefig('plots/test_hankel.pdf')

if __name__ == '__main__':
    print('\n\nWelcome to linear xi.\n\n')
    
    rs = np.arange(1.e-3, 5.e2, 0.1)
    
    for zz in np.arange(2, 5, 1):
        xi = get_linearxi(zz)

        pl.loglog(rs, xi(rs), '-', label=r'$z={:.2f}$'.format(zz))

    pl.legend(frameon=False)

    pl.savefig('plots/linear_xi.pdf')
    
    print('\n\nDone.\n\n')
