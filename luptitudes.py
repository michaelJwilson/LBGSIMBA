import numpy as np
import pylab as pl


@np.vectorize
def luptitude(F, SigF):
  ##  Luptitudes, https://arxiv.org/pdf/astro-ph/9903081.pdf
  ##  Note:  x = F / F0, where F0 is the flux of a zero mag. object.
  a     = 2.500 * np.log10(np.exp(1))
  F0    = 10. ** (-48.60 / 2.5)

  x     = F / F0
  b     = 1.042 * (SigF / F0)

  return  -a * (np.arcsinh(x / 2. / b) + np.log(b)) 

def lup_lim(SigF):
  a     = 2.500 * np.log10(np.exp(1))
  F0    = 10. ** (-48.60 / 2.5)
  
  b     = 1.042 * (SigF / F0)

  return  -a * np.log(b)


if __name__ == '__main__':
  print('\n\nWelcome to Luptitudes.\n\n')

  mag25  = 10. ** (-(25. + 48.60) / 2.5)
  SigF   = mag25 / 5.

  ms     = np.arange(20., 35., 0.1)

  Fs     = np.array([10. ** (-(x + 48.60) / 2.5) for x in ms])
  Lups   = luptitude(Fs, SigF)

  lim    = lup_lim(SigF)

  pl.plot(ms, Lups, 'k-')

  pl.axvline(x=25., ymin=0., ymax=1.)
  pl.axhline(y=lim, xmin=0., xmax=1.)

  pl.xlim(24., 32.)
  pl.ylim(24., 27.)

  pl.show()

  print('\n\nDone.\n\n')
