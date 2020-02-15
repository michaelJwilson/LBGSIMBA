import numpy  as      np
import pylab  as      pl

from   mcfit             import  P2xi, TophatVar
from   scipy.interpolate import  CubicSpline


k, P, hf  = np.loadtxt('dat/white/pklin_z300.txt', unpack=True)
r, xi     = P2xi(k)(P)

R, var    = TophatVar(k, lowring=True)(P, extrap=True)
varR      = CubicSpline(R, var)
sigma8    = np.sqrt(varR(8))

pl.loglog(r, xi, label=r'$\sigma_8: {:.3f}$'.format(sigma8))
pl.legend(loc=2, frameon=False)
pl.savefig('plots/test_hankel.pdf')
