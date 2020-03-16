import  matplotlib;                matplotlib.use('PDF')

import  os
import  time
import  scipy
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt
import  astropy.units      as      u
import  astropy.constants  as      const

from    scipy              import  signal
from    scipy.special      import  gamma  as gammaf
from    sutils             import  latexify
from    cosmo              import  cosmo
from    astropy.cosmology  import  z_at_value as _z_at_value
from    scipy.interpolate  import  interp1d
from    scipy.integrate    import  quad, quadrature
from    pz.hildebrandt_09  import  getpz_H09
from    xi_hankel          import  get_linearxi 


latexify(columns=1, equal=True, fontsize=12, ggplot=True, usetex=True)

MAXITER  = 500
colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']

@np.vectorize
def Gaussian_pz(z, z0=4., sigma=0.5):    
    # Normalised Gaussian in z. 
    norm = 1. / sigma / np.sqrt(2. * sigma)
    exp  = (z - z0) / sigma
    exp  = exp * exp
    exp /= -2.
    
    return  np.exp(exp) / norm

def pc(chi, pz):
    zee  = z_at_value(chi, cosmo.comoving_distance)    

    return  pz(zee) * 100. * cosmo.efunc(zee) / const.c.to('km/s').value

def tophatc(chi, rc=2.e3, dr=1.e2):
    norm = 1. / (2. * dr)
    
    on   = np.abs(chi - rc) <= dr  
    
    return  norm * on.astype(np.float)
    
def xi(r, r0=5.0, gamma=1.8):
    # r0 in Mpc/h
    return  (r / r0) ** -gamma

def Aw(rc, dr, gamma, r0): 
    # narrow slice, eqn. (19) of Simon.
    _     = rc ** (1. - gamma) / 2. / dr

    norm  = np.sqrt(np.pi) * r0 ** gamma
    norm *= gammaf(gamma / 2. - 1. / 2.)
    norm /= gammaf(gamma / 2.)

    return  _ * norm

def beta2gamma(beta):
    ##  -beta  = 1. - gamma
    return  1. + beta

def pow_wtheta(theta, rc=None, dr=None, gamma=1.8, r0=5., Aw=None):
    ##  theta in radians. 
    if Aw is None:
      Aw = Aw(rc, dr, gamma, r0)

    if (rc is None) & (dr is None) & (Aw is None):
        raise ValueError('Provide either Aw or rc and dr.')
      
    return  Aw * theta ** (1. - gamma)

@np.vectorize
def inner(theta, rbar, xi):
    def _(lndr):
        dr = np.exp(lndr)

        ##  Eqn. (5) of www.aanda.org/articles/aa/pdf/2007/39/aa6352-06.pdf
        R  = np.sqrt(rbar**2. * theta**2. + dr**2.)

        return  dr * xi(R)
   
    ##  scipy.integrate.quad(_, np.log(Rmin), np.log(30.), limit=MAXITER)[0] 
    I = scipy.integrate.quadrature(_, np.log(1.e-8), np.log(1.e4), tol=1.0e-05, rtol=1.0e-04, maxiter=MAXITER, vec_func=True)[0]
    
    return  2. * I

@np.vectorize
def limber_wtheta(theta, p1, p2, xi, zmin=1.0, zmax=5.0):
    print('Solving for theta: {} arcsec.'.format(3600. * theta * 180. / np.pi))

    def _(z):
      x  = cosmo.h * cosmo.comoving_distance(z).value           # [Mpc/h].
      Hz =    100. * cosmo.efunc(z) / const.c.to('km/s').value  # [(Mpc/h)^-1].

      return  p1(z) * p2(z) * inner(theta, x, xi)  
    
    I    = scipy.integrate.quadrature(_, zmin, zmax, tol=1.0e-05, rtol=1.0e-04, maxiter=MAXITER, vec_func=True)[0]
    
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

def plot_wtheta(park=False, qiu=False):    
    zs = [2.024621, 3.00307, 3.963392, 5.0244]
    
    for i, zz in enumerate(zs):
      ts, lwt, wt = np.loadtxt('dat/wtheta_{:.3f}'.format(zz).replace('.', 'p') + '.txt', unpack=True)

      prefac = 60. * 60.
        
      pl.loglog(prefac * ts, lwt, c=colors[i],  label='$z$={:.2f}'.format(zz), alpha=1., linestyle='-', lw=0.3)

    if park:
        Aw, beta = np.loadtxt('dat/park16_wtheta.dat', unpack=True, usecols=(0,1))
        gamma    = beta2gamma(beta)

        npark    = len(Aw)
        
        for j, (A, g) in enumerate(zip(Aw, gamma)):
          wt     = pow_wtheta(3600. * ts, rc=None, dr=None, gamma=g, r0=5., Aw=A)
          
          if j == 0:
              label = 'Park16: $z\simeq 4.0$'

          else:
              label = ''
              
          pl.loglog(prefac * ts, wt, linestyle='dashdot', c=colors[2], label=label)

    if qiu:
        Aw, beta = np.loadtxt('dat/qiu18_wtheta.dat', unpack=True, usecols=(0,1))
        gamma    = beta2gamma(beta)

        nqiu     = len(Aw)
        
        for j, (A, g) in enumerate(zip(Aw, gamma)):
          wt     = pow_wtheta(3600. * ts, rc=None, dr=None, gamma=g, r0=5., Aw=A)

          if j == 0:
              label = 'Qiu18: $z\simeq 3.8$'

          else:
              label = ''

          pl.loglog(prefac * ts, wt, linestyle='dotted', c=colors[2], label=label)

    ##
    axes    = pl.gca()
    lines   = axes.get_lines()

    legend1 = pl.legend([lines[i] for i in [0,1,2,3]], ['$z$={:.2f}'.format(zz) for zz in zs], loc=1, frameon=False)
    legend2 = pl.legend([lines[i] for i in [4, 4 + npark]], ['Park16: $z\simeq 4.0$', 'Qiu18: $z\simeq 3.8$'], loc=3, frameon=False)
    
    
    axes.add_artist(legend1)
    axes.add_artist(legend2)
            
    ##  pl.legend(frameon=False, loc=3)
    
    pl.xlabel(r'$\theta$ [arcsec.]')                                                                                                                                                                                         
    pl.ylabel(r'$\omega(\theta$)')                                                                                                                                                                                             
    plt.tight_layout()                                                                                                                                                                                                   
    pl.savefig('plots/wtheta.pdf')   

def test_xi():
    zs  = np.arange(0.1, 0.2, 0.1)

    for i, redshift in enumerate(zs):
      xi = get_linearxi(redshift, root=os.environ['LBGSIMBA'] + '/white/dat/', halofit=True)

      rs = np.logspace(-3., 3., 500.)
      xs = xi(rs)
      
      pl.loglog(rs, np.abs(xs), '-', label=redshift)

    pl.legend(frameon=False)
    pl.savefig('plots/test_xi.pdf')

def test_innera():
  def _(lndr, floor):
      dr = np.exp(lndr)
      
      ##  Eqn. (5) of www.aanda.org/articles/aa/pdf/2007/39/aa6352-06.pdf
      R  = np.sqrt(floor + dr**2.)

      return  dr * xi(R)

  xi   = get_linearxi(3.00307, root=os.environ['LBGSIMBA'] + '/white/dat/', halofit=True)

  lndr = np.log(np.logspace(np.log(1.e-8), np.log(1.e4), 1000., base=np.exp(1)))
  
  for floor in [0.0, 10., 100., 1000.]:
    Ps   = _(lndr, floor)
    pl.loglog(np.exp(lndr), np.abs(Ps), '-', label=floor)

  pl.xlim(1.e-8, 1.e4)
    
  pl.xlabel(r'$\Delta$')
  pl.title('I={}'.format(I))
  pl.legend(frameon=False, loc=4, ncol=2)
  pl.savefig('plots/test_inner_wtheta.pdf')

def test_innerb():
  ts  = np.logspace(np.log10(6.), np.log10(600.), 25)         # arcseconds                                                                                                                                                         
  ts /= 3600.                                                 # degs.                                                                                                                                                              
  cs  = ts * np.pi / 180.                                     # radians.      

  xi  = get_linearxi(3.00307, root=os.environ['LBGSIMBA'] + '/white/dat/', halofit=True)
  
  for z in np.arange(0.1, 1.0, 0.1):
    chi = cosmo.h * cosmo.comoving_distance(z).value          # [Mpc/h].
    I   = inner(cs, chi, xi)
  
    print('{} \t {}'.format(z, I))
    
    
if __name__ == '__main__':    
    #  https://arxiv.org/pdf/astro-ph/0609165.pdf    
    ts  = np.logspace(np.log10(6.), np.log10(3600.), 25)      # arcseconds
    ts /= 3600.                                               # degs. 
    cs  = ts * np.pi / 180.                                   # radians. 

    ##  hildebrandt_pz = getpz_H09(sample='u', interp=True)

    zs  = [2.024621, 3.00307,  3.963392, 5.0244]
    bs  = [2.200000, 3.200000, 4.200000, 5.5000]
    
    for i, redshift in enumerate(zs):
      ##  zz, rs, _  = zel(i)

      xi             = get_linearxi(redshift, root=os.environ['LBGSIMBA'] + '/white/dat/', halofit=True)

      def bxi(r):
          return  bs[i] * bs[i] * xi(r)

      def pp(z):
          return  Gaussian_pz(z, z0=redshift, sigma=0.5)
      
      ##  Ts         = test_inner(cs, redshift, xi)

      ##  lwt        = limber_wtheta(cs, tophatc, tophatc, xi)
      lwt            = limber_wtheta(cs, pp, pp, bxi)
      
      ##  wt         = pow_wtheta(cs, 2.e3, 1.e2)
      wt             = np.zeros_like(lwt)

      np.savetxt('dat/wtheta_{:.3f}'.format(redshift).replace('.', 'p') + '.txt', np.c_[ts, lwt, wt], fmt='%.6le')
    
    plot_wtheta(park=True, qiu=True)
    
    print('\n\nDone.\n\n')
