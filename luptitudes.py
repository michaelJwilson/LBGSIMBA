import matplotlib; matplotlib.use('PDF')
import numpy as np
import pylab as pl

from   hildebrandt  import  onesig_mag


@np.vectorize
def luptitude(F, SigF):
  ##  Luptitudes, https://arxiv.org/pdf/astro-ph/9903081.pdf
  ##  Note:  x = F / F0, where F0 is the flux of a zero mag. object.
  a  = 2.5 * np.log10(np.exp(1))
  F0 = 10. ** (-48.60 / 2.5)

  x  = F / F0
  
  if np.isinf(SigF):
    return  np.inf

  # b   = 1.042 * (SigF / F0)

  # Better over / under flow.
  lnb   = np.log(1.042) + np.log(SigF) - np.log(F0)

  xonb  = F / SigF / 1.042 

  return  -a * (np.arcsinh(xonb / 2.) + lnb) 

def lup_lim(SigF):
  ##  x = 0 Luptitude limit.   
  a     = 2.500 * np.log10(np.exp(1))
  F0    = 10. ** (-48.60 / 2.5)
  
  lnb   = np.log(1.042) + np.log(SigF) - np.log(F0)

  return  -a * lnb


if __name__ == '__main__':
  from  hildebrandt import ferr


  print('\n\nWelcome to Luptitudes.\n\n')
  
  fivesig = 25.30
  onesig  = onesig_mag(fivesig)           ##  AB mags. 
  onesig  = 10. ** -((onesig + 48.60) / 2.5)
  
  print(onesig)
  
  ##
  ms      = np.arange(20., 45., 0.2)

  Fs      = np.array([10. ** (-(x + 48.60) / 2.5) for x in ms])
  sigFs   = np.array([onesig] * len(Fs))

  ##  luptitude(Fs,   sigFs)
  Lups    = luptitude(Fs, sigFs)
  lims    = lup_lim(sigFs)
  
  print(Fs)
  print(sigFs)
  print(Lups)
  
  # lim     = lup_lim(SigF)

  pl.plot(ms, ms,   c='gold', label=r'')
  pl.plot(ms, Lups,     'k-', label=r'$5\sigma depth: %.2f$' % fivesig)
  pl.plot(ms, lims,     'r-')
  
  pl.axvline(x=fivesig, ymin=0., ymax=1.)
  # pl.axhline(y=lim, xmin=0., xmax=1.)

  pl.xlabel('Source AB magnitude')
  pl.ylabel('Source AB Luptitude')
  
  pl.xlim( 20., 46.)
  # pl.ylim(-20.,  29.)

  pl.savefig('plots/luptitudes.pdf')

  print('\n\nDone.\n\n')
