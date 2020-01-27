import numpy as np

from spherical_bessel_transform import SphericalBesselTransform
from loginterp import loginterp

from cleft_fftw import CLEFT


class CorrelationFunction(CLEFT):
    '''
    Class to calculate the real space correlation function given a linear power spectrum and some bias parameters.
    
    Inherits the the CLEFT class.
    
    '''

    def __init__(self, *args, kmin=3e-4, kmax=3, nk=200, **kw):
        '''
         Same keywords and arguments as the other two classes for now.
        '''
        
        # Setup ffts etc.
        CLEFT.__init__(self, *args, **kw)

        self.kmin, self.kmax, self.nk = kmin, kmax, nk
        self.kv = np.logspace(np.log10(kmin), np.log10(kmax), nk); self.nk = nk
        
        self.kint = np.logspace(-5,3,4000)
        self.sph_gsm  = SphericalBesselTransform(self.kint,L=3,fourier=True)
        self.rint = np.logspace(-3,5,4000)
        self.rint = self.rint[(self.rint>0.1)*(self.rint<600)] #actual range of integration
        
        self.setup_power_spectrum()
        self.setup_config()

    def setup_power_spectrum(self):
        self.make_ptable(kmin = self.kmin, kmax = self.kmax, nk = self.nk)

    def setup_config(self):
        # Fourier transform power spectrum template
        
        # the correlation function
        self.xitable = np.zeros((len(self.rint),12))

        for ii in range(12):
            _integrand = loginterp(self.pktable[:,0], self.pktable[:,1+ii])(self.kint)
            qs, xs = self.sph_gsm.sph(0,_integrand)
            self.xitable[:,ii] = np.interp(self.rint, qs, xs)
            
        _integrand = loginterp(self.pktable[:,0], self.pktable[:,1]+self.pktable[:,2]+self.pktable[:,3])(self.kint)
        qint, ximatter = self.sph_gsm.sph(0,_integrand)
        self.ximatter = np.interp(self.rint, qint, ximatter)

        _integrand = loginterp(self.pktable[:,0], self.pktable[:,0]**2 * self.pktable[:,1])(self.kint)
        qint, xict = self.sph_gsm.sph(0,_integrand)
        self.xict = np.interp(self.rint, qint, xict)
        

    def combine_bias_terms(self, b1, b2, bs, alpha):
        '''
        Calculate velocity moments and turn into cumulants.
        The bvec format is [b1, b2, bs, alpha, alpha_v, alpha_s0, alpha_s2]
        
        '''
        # Compute each moment
        self.xieft = self.ximatter + b1*self.xitable[:,3] + b1**2*self.xitable[:,4]\
        + b2*self.xitable[:,5] + b1*b2*self.xitable[:,6] + b2**2 * self.xitable[:,7]\
        + bs*self.xitable[:,8] + b1*bs*self.xitable[:,9] + b2*bs*self.xitable[:,10]\
        + bs**2*self.xitable[:,11] + alpha*self.xict
        
        return self.rint, self.xieft
