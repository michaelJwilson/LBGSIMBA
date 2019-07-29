import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy                 as      np
import  pylab                 as      pl
import  matplotlib.pyplot     as      plt

from    scipy.spatial         import  KDTree 
from    itertools             import  product
from    main                  import  get_data, print_keys
from    utils                 import  latexify


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

print('\n\nWelcome to Simba stellar mass.')
    
##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize     =  100.
getredshift =  3.00307

f, p        =  get_data(boxsize, getredshift)

##
gid         =  f['galaxy_data']['GroupID'][:]
iscentral   =  f['galaxy_data']['central'][:]
haloindex   =  f['galaxy_data']['parent_halo_index'][:]
    
##  Basic stats.
print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
print('Number of centrals found: {}'.format(np.sum(iscentral)))
print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))

##  
sfr         =  f['galaxy_data']['sfr'][:]            ##  [Solar mass per year].
bhmdot      =  f['galaxy_data']['bhmdot'][:]
gfrac       =  f['galaxy_data']['gas_fraction'][:]

sfr        *=  1.e9                                  ##  [Solar mass per giga year].

smass       =  1.82e7                                ##  Stellar mass (particle) resolution [Solar masses].
stellarmass =  f['galaxy_data']['nstar'][:] * smass
ssfr        =  sfr / stellarmass


print('\n\nRange of stellar mass: {} to {} solar masses.'.format(stellarmass.min() / 1e10, stellarmass.max() / 1e10))
print('Range of SFR: {} to {}.'.format(sfr.min(), sfr.max()))
print('Range of SSFR: {} to {}.'.format(ssfr.min(), ssfr.max()))

##  Stellar mass function.
bins          =  np.logspace(9., 13., 10, endpoint=True, base=10.0)
bsmass        =  np.digitize(stellarmass, bins=bins)

ninds, counts =  np.unique(bsmass, return_counts = True)
mean_smass    =	 np.array([np.mean(stellarmass[bsmass == _bin]) for _bin in np.arange(len(bins))])
mean_ssfr     =  np.array([np.mean(ssfr[bsmass == _bin]) for _bin in np.arange(len(bins))])

pl.plot(np.log10(mean_smass), np.log10(mean_ssfr), c='darkcyan', lw=1)

pl.ylim(-1., 1.)

pl.xlabel(r'$\log_{10}|M_*|$')
pl.ylabel(r'$\log_{10}|\dot M_* / M_* / \rm{Gyr}^{-1} |$')

plt.tight_layout()

pl.savefig('plots/ssfr.pdf')

print('\n\nDone.\n\n')
    
