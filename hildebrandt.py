import  matplotlib;  matplotlib.use('PDF')
import  numpy as     np

from    utils import latexify


@np.vectorize
def merr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22):
  ##  Magnitude error.
  if m < mstar:
    return  estar * 10. ** (0.4 * (alphab + 1.) * (m - mstar))

  else:
    return  estar * np.exp(10. ** (alphaf * (m - mstar))) / 2.72 

def ferr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None):
  flux = 10. ** -((m + 48.60) / 2.5)  ##  Assumes AB mag.  [erg/s/cm2/Hz]. 
  
  sigm = merr(m, mstar, estar=estar, alphab=alphab, alphaf=alphaf)

  ##  E.g. eqn. (6) of https://arxiv.org/pdf/1509.00870.pdf
  sigf = flux * sigm * np.log(10.) / 2.5

  if lim_snr is not None:
    sigf = np.sqrt(sigf**2. + (flux / lim_snr) ** 2.)

  return  sigf

def snr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None):
  flux   = 10. ** -((m + 48.60) / 2.5)  ##  Assumes AB mag.  [erg/s/cm2/Hz].

  ##  Depth map.
  sigf   = ferr(m, mstar, estar=estar, alphab=alphab, alphaf=alphaf)
  
  if lim_snr is not None:
    sigf = np.sqrt(sigf**2. + (flux / lim_snr) ** 2.)

  snr = flux / sigf

  return  snr

def onesig_mag(mstar, alphaf=0.22, estar=0.2):
  #  Return the magnitude at which a source would have S/N of unity in
  #  the PHAT model.
  x  = np.log(1.086 * 2.72 / estar)
  x  = np.log10(x) / alphaf

  return  mstar + x
  
  
if __name__ == '__main__':
  import  matplotlib;  matplotlib.use('PDF')
  import  pylab as pl


  latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)

  ##
  ms   = np.arange(16., 31., 0.1)
    
  for mstar in [28., 27., 26., 25., 24.]:
    Fs =  10. ** -((ms    + 48.60) / 2.5)
    F0 =  10. ** -((mstar + 48.60) / 2.5)
    
    pl.semilogy(ms, 5. * (Fs / F0)**0.5, 'gold')

    pl.semilogy(ms, snr(ms, mstar),  label=r'$m_*$ = ' + str(mstar))
    pl.semilogy(ms, snr(ms, mstar,  lim_snr=100.), '--')

    pl.axvline(x=mstar, ymin=0., ymax=1., c='k', alpha=0.3)
    
  pl.axhline(y=5.0, xmin=0, xmax=1, c='k', alpha=0.3)
  pl.axhline(y=1.0, xmin=0, xmax=1, c='k', alpha=0.3)
    
  pl.legend(loc=3, frameon=False)

  pl.ylim(0.1, 1.e3)

  pl.xlabel('$m$')
  pl.ylabel(r'$S/N$')

  pl.savefig('plots/Hildebrandt.pdf')

    
