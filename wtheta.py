import matplotlib;                matplotlib.use('PDF')

import time
import scipy
import numpy              as      np
import pylab              as      pl
import matplotlib.pyplot  as      plt
import astropy.units      as      u
import astropy.constants  as      const

from   scipy              import  signal
from   scipy.special      import  gamma  as gammaf
from   utils              import  latexify
from   cosmo              import  cosmo
from   astropy.cosmology  import  z_at_value as _z_at_value
from   scipy.interpolate  import  interp1d
from   scipy.integrate    import  quad, quadrature
from   pz.hildebrandt_09  import  getpz_H09
from   xi_hankel          import  get_linearxi


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

MAXITER  = 500
colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']

def p(z, z0=2., sigma=0.5):    
    # Normalised Gaussian in z. 
    norm = 1. / sigma / np.sqrt(2. * sigma)
    exp  = (z - z0) / sigma
    exp  = exp * exp
    exp /= -2.
    
    return  np.exp(exp) / norm

def z_at_value(vals, func):
    vals = np.array([vals])    
    return np.array([_z_at_value(func, x * u.Mpc / cosmo.h, zmax=10.) for x in vals])

def pc(chi, pz):
    zee  = z_at_value(chi, cosmo.comoving_distance)    

    return  pz(zee) * 100. * cosmo.efunc(zee) / const.c.to('km/s').value

def hildebrandt_pz(chi, sample='u'):
    return  pc(chi, getpz_H09(sample=sample, interp=True))

def tophatc(chi, rc=2.e3, dr=1.e2):
    norm = 1. / (2. * dr)
    
    on   = np.abs(chi - rc) <= dr  
    
    return  norm * on.astype(np.float)
    
def xi(r, r0=5., gamma=1.8):
    # r0 in Mpc/h
    return  (r / r0) ** -gamma

def Aw(rc, dr, gamma, r0): 
    # narrow slice, eqn. (19) of Simon.
    _     = rc ** (1. - gamma) / 2. / dr

    norm  = np.sqrt(np.pi) * r0 ** gamma
    norm *= gammaf(gamma / 2. - 1. / 2.)
    norm /= gammaf(gamma / 2.)

    return  _ * norm
        
def pow_wtheta(theta, rc, dr, gamma=1.8, r0=5.):
    ##  theta in radians. 
    return  Aw(rc, dr, gamma, r0) * theta ** (1. - gamma)

def inner(theta, rbar, xi):
    Rmin = rbar	* theta

    def _(dlnr):
        dr = np.exp(dlnr)

        ##  Eqn. (5) of www.aanda.org/articles/aa/pdf/2007/39/aa6352-06.pdf
        R  = np.sqrt(rbar**2. * theta**2. + dr**2.)

        return  xi(R) / dr
   
    ##  scipy.integrate.quad(_, np.log(Rmin), np.log(30.), limit=MAXITER)[0] 
    I = scipy.integrate.quadrature(_, np.log(Rmin), np.log(30.), tol=1.5e-05, rtol=1.5e-05, maxiter=MAXITER, vec_func=True)[0]
    
    return  2. * I

@np.vectorize
def limber_wtheta(theta, p1, p2, xi, zmin=1.0, zmax=5.0):
    def _(lnx):
      x = np.exp(lnx)
        
      return  p1(x) * p2(x) * inner(theta, x, xi)

    print('Solving for theta: {} arcsec.'.format(3600. * theta))
    
    lower = cosmo.h * cosmo.comoving_distance(zmin).value  # Mpc/h 
    upper = cosmo.h * cosmo.comoving_distance(zmax).value  # Mpc/h

    ##  scipy.integrate.quad(_, np.log(lower), np.log(upper), limit=MAXITER)[0]
    I     = scipy.integrate.quadrature(_, np.log(lower), np.log(upper), tol=1.5e-08, rtol=1.5e-08, maxiter=MAXITER, vec_func=False)[0]
    
    return  I

def zel(index):
    redshifts = [2.024621, 3.00307, 3.963392, 5.0244]
    redshift  = redshifts[index]
    
    #  r [Mpc/h]    r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2
    iz        = int(100 * redshift + 0.001)
    _         = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

    b1, b2    = 1.0, 0.0

    cc        = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])

    rs        = _[:,0]
    _         = _[:, 0:8]

    result    = np.dot(_, cc)

    return  redshift, rs, interp1d(rs, result, fill_value=0.0, bounds_error=False)

def plot_wtheta():    
    for i, zz in enumerate([2.024621, 3.00307, 3.963392, 5.0244]):
        ts, lwt, wt = np.loadtxt('dat/wtheta_{:.3f}'.format(zz).replace('.', 'p') + '.txt', unpack=True)

        prefac = 60. * 60.
        
        if i == 0:
            pl.loglog(prefac * ts, lwt, c=colors[i],  label='Limber', alpha=0.5)

        else:
            pl.loglog(prefac * ts, lwt, c=colors[i],  label='', alpha=0.5)

        break

    pl.legend(frameon=False)
    
    pl.xlabel(r'$\theta$ [arcsec.]')                                                                                                                                                                                         
    pl.ylabel(r'$\omega(\theta$)')                                                                                                                                                                                             

    plt.tight_layout()                                                                                                                                                                                                   

    pl.savefig('plots/wtheta.pdf')   


if __name__ == '__main__':    
    #  https://arxiv.org/pdf/astro-ph/0609165.pdf    
    ts = np.arange(0.001, 0.02, 0.0001)  # degs.
    cs = ts * np.pi / 180.               # radians. 
    
    for i, redshift in enumerate([2.024621, 3.00307, 3.963392, 5.0244]):
      zz, rs, _  = zel(i)

      xi         = get_linearxi(redshift)

      ##  lwt    = limber_wtheta(cs, tophatc, tophatc, xi)
      lwt        = limber_wtheta(cs, hildebrandt_pz, hildebrandt_pz, xi)
      
      wt         =    pow_wtheta(cs, 2.e3, 1.e2)

      np.savetxt('dat/wtheta_{:.3f}'.format(zz).replace('.', 'p') + '.txt', np.c_[ts, lwt, wt], fmt='%.6le')

      break
      
    plot_wtheta()
    
    print('\n\nDone.\n\n')
