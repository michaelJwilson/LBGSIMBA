import matplotlib; matplotlib.use('PDF')

import scipy
import numpy              as     np
import pylab              as     pl
import matplotlib.pyplot  as     plt
import astropy.units      as     u
import astropy.constants  as     const

from   scipy              import signal
from   scipy.special      import gamma  as gammaf
from   utils              import latexify
from   cosmo              import cosmo
from   astropy.cosmology  import z_at_value as _z_at_value
from   scipy.interpolate  import interp1d
from   scipy.integrate    import quad


latexify(fontsize=6)

def p(z, z0, sigma):    
    # Normalised Gaussian in z. 
    norm = 1. / sigma / np.sqrt(2. * sigma)
    exp  = (z - z0) / sigma
    exp  = exp * exp
    exp /= -2.
    
    return  np.exp(exp) / norm

def z_at_value(vals, func):
    vals = np.array([vals])    
    return np.array([_z_at_value(func, x * u.Mpc / cosmo.h, zmax=10.) for x in vals])

def pc(chi, z0=2., sigma=0.5):
    zee    = z_at_value(chi, cosmo.comoving_distance)    

    return p(zee, z0, sigma) * 100. * cosmo.efunc(zee) / const.c.to('km/s').value
    
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

def inner(theta, rbar, xi):
    def _(dr):
        R = np.sqrt(rbar**2. * theta**2. + dr**2.)

        return  xi(R)

    # NOTE:  Ranges on dr of inner integral. 
    return  scipy.integrate.quad(_, 1.0, 250., limit=500)[0]

@np.vectorize
def limber_wtheta(theta, p1, p2, xi):
    def p1p2(x):
        return p1(x) * p2(x) * inner(theta, x, xi)

    # NOTE:  Integral over \bar r.  Range hard coded.
    return  scipy.integrate.quad(p1p2, 10., 6000., limit=500)[0]

def zel(index):
    redshifts = [2.024621, 3.00307, 3.963392, 5.0244]
    redshift  = redshifts[index]
    
    # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2
    iz        = int(100 * redshift + 0.001)
    _         = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

    b1, b2    = 1.0, 0.0

    cc        = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])

    rs        = _[:,0]
    _         = _[:,0:8]

    result    = np.dot(_, cc)

    return  redshift, rs, interp1d(rs, result, fill_value=0.0, bounds_error=False)


if __name__ == '__main__':    
    # https://arxiv.org/pdf/astro-ph/0609165.pdf
    
    ts   = np.arange(0.1, 10., 1.0)  # degs.
    cs   = ts * np.pi / 180.         # radians. 

    # wt = pow_wtheta(cs, 1.8, 2000., 100., 5.)  
    
    for i in np.arange(4):
      zz, rs, xi = zel(i)      
      wt         = limber_wtheta(cs, pc, pc, xi)

      pl.loglog(ts, wt, label=zz)

    pl.xlabel(r'$\theta$ [deg.]')
    pl.ylabel(r'$\omega(\theta$)')

    plt.tight_layout()
    
    pl.savefig('plots/wtheta.pdf')
    
    print('\n\nDone.\n\n')
