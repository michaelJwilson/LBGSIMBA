import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  pandas                as      pd
import  numpy                 as      np
import  pylab                 as      pl
import  matplotlib.pyplot     as      plt

from    scipy.spatial         import  KDTree 
from    itertools             import  product
from    get_data              import  get_data, print_keys, get_phys
from    utils                 import  latexify


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)


print('\n\nWelcome to Simba stellar mass.')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

##  Closest redshifts:  [2.024621, 3.00307, 3.963392, 5.0244].
for boxsize, linestyle, min_smass, label in zip([100., 50.], ['-', '--'], [1.82e7, 2.82e6], [1, 0]):
  vol             =  boxsize ** 3.

  for i, (getredshift, name) in enumerate(zip([2.024621, 3.00307, 3.963392, 5.0244], ['two', 'three', 'four', 'five'])):
    color         =  colors[i]

    stellarmass   =  get_phys(boxsize, getredshift, printit=False, nrows=-1)['stellarmass']
    ssfr          =  get_phys(boxsize, getredshift, printit=False, nrows=-1)['ssfr']
    
    print('\n\nRange of stellar mass: {} to {} solar masses.'.format(stellarmass.min() / 1e10, stellarmass.max() / 1e10))
    print('\n\nRange of SSFR:         {} to {} per Gyr.'.format(ssfr.min(), ssfr.max()))
    
    ##  Read sample selection.                                                                                                                                                        
    lsst_sample   =  pd.read_pickle("bigdat/{}.pkl".format(name))
    isin          =  lsst_sample['INSAMPLE'].values

    for x, y, alpha, label in zip([stellarmass, stellarmass[isin]], [ssfr, ssfr[isin]], [0.5, 1.0], ['', r'$z$ = %.2lf' % getredshift]):
      ##  Stellar mass function.
      bins         =  np.arange(12.25, 8.75, -0.5)

      ssfr         =  np.log10(y)
      stellarmass  =  np.log10(x)

      ssfr         =         ssfr[stellarmass <= bins[0]]
      stellarmass  =  stellarmass[stellarmass <= bins[0]]
      
      ## 
      bsmass       =  np.digitize(stellarmass, bins=bins, right=True)
      
      mean_ssfr    =  np.array([np.mean(       ssfr[bsmass == _bin]) for _bin in np.arange(len(bins))])
      mean_smass   =  np.array([np.mean(stellarmass[bsmass == _bin]) for _bin in np.arange(len(bins))])

      if not label:
        label = ''
      
      pl.plot(mean_smass, mean_ssfr, lw=1, alpha=alpha, label=label, linestyle=linestyle, color=color)

  break

##
pl.xlim(8.50, 12.5)
pl.ylim(-.5,   1.5)

pl.xlabel(r'$\log_{10}|M_*|$')
pl.ylabel(r'$\log_{10}|\dot M_* / M_* / \rm{Gyr}^{-1}|$')

plt.tight_layout()

pl.legend(loc=1, frameon=False, handlelength=1, borderaxespad=1.)

pl.savefig('plots/ssfr.pdf')

pl.clf()

print('\n\nDone.\n\n')
    
