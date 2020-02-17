import numpy as np
import pylab as pl
import os


def getpz_H09(sample='u'):
  # Samples at supplied zs (every dz = 0.1)
  if sample == 'u':
    data  = np.loadtxt('dat/udrops_r23p0t25p5_Masim.dat')
    
  elif sample == 'g':
    data  = np.loadtxt('dat/gdrops_i23p5t26p0_Masim.dat')

  else:
    raise ValueError('Only u and g drops are aviailable for CARS currently.')

  dzs   = np.diff(data[:,0]).astype(np.float32)

  assert  np.all(dzs == dzs[0])

  dz    = data[1,0] - data[0,0]
  norm  = dz * np.sum(data[:,1])

  pz    = data[:,1] / norm
  midz  = data[:,0] + 0.5 * dz

  assert  np.isclose(dz * np.sum(pz), 1.00)

  return  midz, pz


if __name__ == "__main__":
  for sample in ['u', 'g']:
    midz, pz = getpz_H09(sample=sample)
    pl.plot(midz, pz, label=r'${}$-dropouts'.format(sample))

  pl.legend(frameon=False)
  pl.show()
