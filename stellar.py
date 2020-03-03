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
from    sutils                import  latexify
from    insample              import  read_insample


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

print('\n\nWelcome to Simba stellar mass.')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
##  pl.axvline(np.log10(100. * 1.82e7), ymin=0., ymax=1., c='k', alpha=1.0, lw=1)
##  pl.axvline(np.log10(100. * 2.28e6), ymin=0., ymax=1., c='k', alpha=1.0, lw=1)

##  pl.axhline(np.log10(1.e-2),  xmin=0., xmax=1., c='k', linestyle='--', alpha=1.0, lw=1)
##  pl.axhline(np.log10(1.e-3),  xmin=0., xmax=1., c='k', linestyle='--', alpha=1.0, lw=1)

for boxsize, linestyle, min_smass, label in zip([100., 50.], ['-', '--'], [1.82e7, 2.82e6], [1, 0]):
  vol           =  boxsize ** 3.

  for i, (getredshift, name) in enumerate(zip([2.024621, 3.00307, 3.963392, 5.0244], ['two', 'three', 'four', 'five'])):
    color       =  colors[i]

    stellarmass =  get_phys(boxsize, getredshift, printit=False, nrows=-1)['stellarmass']
    
    print('\n\nRange of stellar mass: {} to {} solar masses.'.format(stellarmass.min() / 1e10, stellarmass.max() / 1e10))

    ##  Read sample selection.                                                                                                                                                        
    lsst_sample   =  read_insample(getredshift)
    isin          =  lsst_sample['INSAMPLE'].values

    for x, alpha, label in zip([stellarmass, stellarmass[isin]], [0.5, 1.0], ['', r'$z$ = %.2lf' % getredshift]):
      ##  Stellar mass function.
      bins         =  np.arange(10.25, 8.75, -0.1)

      stellarmass  =  np.log10(x)
    
      stellarmass  =  stellarmass[stellarmass <= bins[0]]

      bsmass       =  np.digitize(stellarmass, bins=bins, right=True)

      ubins, cnts  =  np.unique(bsmass, return_counts = True)

      cnts         =   cnts[ubins < len(bins)]
      ubins        =  ubins[ubins < len(bins)]

      missing      =  np.setdiff1d(np.arange(len(bins)), ubins)

      ubins        =  np.concatenate([ubins, missing])
      cnts         =  np.concatenate([cnts,  np.zeros_like(missing)])

      sortind      =  np.argsort(ubins)

      ubins        =  ubins[sortind]
      cnts         =   cnts[sortind]

      mean_smass   =  np.array([np.mean(stellarmass[bsmass == _bin]) for _bin in np.arange(len(bins))])

      assert len(ubins) == len(bins)

      if not label:
        label = ''
      
      pl.plot(mean_smass, np.log10(np.cumsum(cnts) / vol), lw=1, alpha=alpha, label=label, linestyle=linestyle, color=color)

  break

##
pl.xlim(8.50, 10.5)
pl.ylim( -5., -1.5)

pl.xlabel(r'$\log_{10}|M_*|$')
pl.ylabel(r'$\log_{10}|\bar n(> M_*) / (h^{-1} \rm{Mpc})^{-3}|$')

plt.tight_layout()

pl.legend(loc=1, frameon=False, handlelength=1, borderaxespad=1.)

pl.savefig('plots/smf.pdf')

pl.clf()

print('\n\nDone.\n\n')
    
