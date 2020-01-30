import matplotlib; matplotlib.use('PDF')

import numpy as     np
import pylab as     pl

from   scipy         import signal
from   scipy.special import gamma  as gammaf

from   utils         import latexify


def p(z, z0, sigma):    
    # Normalised Gaussian in z. 
    norm = 1. / sigma / np.sqrt(2. * sigma)
    exp  = (z - z0) / sigma
    exp  = exp * exp
    exp /= -2.
    
    return  np.exp(exp) / norm

def xi(r, r0=5., gamma=1.8):
    # r0 in Mpc/h
    
    return (r / r0) ** -gamma

def Aw(gamma, rc, dr, r0): 
    # narrow slice, eqn. (19) of Simon.
    _     = rc ** (1. - gamma) / 2. / dr

    norm  = np.sqrt(np.pi) * r0 ** gamma
    norm *= gammaf(gamma / 2. - 1. / 2.)
    norm /= gammaf(gamma / 2.)

    return  _ * norm
        
def pow_wtheta(theta, gamma, rc, dr, r0):
    ##  theta in radians. 
    return Aw(gamma, rc, dr, r0) * theta ** (1. - gamma)

def limber():
    # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2
    
    for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
        # ZA                                                                                                                                                                                                                                                                              
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

	b1, b2 = 1.0, 0.0

        cc     = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])

        rs     = _[:,0]
        _      = _[:,0:8]

        # Zeldovich correlation fn. 
        result = np.dot(_, cc)

    
if __name__ == '__main__':    
    # https://arxiv.org/pdf/astro-ph/0609165.pdf

    # zs = np.arange(0, 10, 0.025)
    # pl.plot(zs, p(zs, 2.0, 0.5))
    # pl.show()

    ts  = np.arange(0.1, 10., 0.01)  ## degs.
    ts *= np.pi / 180.              ## radians. 

    wt  = pow_wtheta(ts, 1.8, 2000., 100., 5.) 

    pl.loglog(ts, wt)
    
    pl.savefig('plots/wtheta.pdf')
    
    print('\n\nDone.\n\n')
