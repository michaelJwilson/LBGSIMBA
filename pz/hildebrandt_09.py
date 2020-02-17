import  numpy  as  np
import  pylab  as  pl
import  os

import  matplotlib.pyplot  as      plt

from    scipy.interpolate  import  interp1d


colors    = plt.rcParams['axes.prop_cycle'].by_key()['color']

def getpz_H09(sample='u', interp=True):
  # Samples at supplied zs (every dz = 0.1)
  if sample == 'u':
    root  = os.environ['LBGSIMBA']
    data  = np.loadtxt(root + '/pz/dat/udrops_r23p0t25p5_Masim.dat')
    
  elif sample == 'g':
    root  = os.environ['LBGSIMBA']
    data  = np.loadtxt(root + '/pz/dat/gdrops_i23p5t26p0_Masim.dat')

  else:
    raise ValueError('Only u and g drops are aviailable for CARS currently.')

  dzs   = np.diff(data[:,0]).astype(np.float32)

  assert  np.all(dzs == dzs[0])

  dz    = data[1,0] - data[0,0]
  norm  = dz * np.sum(data[:,1])

  pz    = data[:,1] / norm
  midz  = data[:,0] + 0.5 * dz

  assert  np.isclose(dz * np.sum(pz), 1.00)

  if interp:
    return  interp1d(midz, pz, copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False, kind='cubic')
    
  else:   
    return  midz, pz


if __name__ == "__main__":
  for i, sample in enumerate(['u', 'g']):
    midz, pz = getpz_H09(sample=sample, interp=False)
    pl.plot(midz, pz, label=r'${}$-dropouts'.format(sample), marker='^', c=colors[i], markersize=3)

    zs       = np.arange(0.0, 10.0, 0.01)
    ps       = getpz_H09(sample=sample, interp=True)(zs) 

    pl.plot(midz, pz, label=r'${}$-dropouts'.format(sample), alpha=0.5, c=colors[i])
    
  pl.legend(frameon=False)
  pl.show()
