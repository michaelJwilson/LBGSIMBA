import  matplotlib

matplotlib.use('PDF')

import  glob
import  h5py
import  numpy          as      np
import  pylab          as      pl

from    scipy.spatial  import  KDTree 
from    itertools      import  product
from    main           import  get_data, print_keys
from    utils          import  latexify


if __name__ == '__main__':
    print('\n\nWelcome to Simba HOD.')

    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244                                                                                                                                                                                                                                                                                                                                        
    boxsize     = 100.
    getredshift = 3.00307

    f, p        =  get_data(boxsize, getredshift)

    ##
    gid         =  f['galaxy_data']['GroupID'][:]
    iscentral   =  f['galaxy_data']['central'][:]
    haloindex   =  f['galaxy_data']['parent_halo_index'][:]
    
    ##  Basic stats.
    print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
    print('Number of centrals found: {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1 - iscentral)))

    ##
    hid         =  f['halo_data']['GroupID'][:]
    ndm         =  f['halo_data']['ndm'][:]

    ##  DM particle mass: 9.6e7 [Solar Mass]
    pmass       =  9.6e7
    hmass       =  ndm[haloindex] * pmass

    ##  Basic halo stats.                                                                                                                                                                                           
    print('\n\nNumber of halos found: {} millions.'.format(len(hid) / 1e6))
    print('Range of halo masses: {} to {} [10^9]'.format(hmass.min() / 1e9, hmass.max() / 1e9))

    ##  Now bin galaxies and halos by halo mass.  
    bins        =  np.logspace(10., 14., 10, endpoint=True, base=10.0)
    bmass       =  np.digitize(hmass, bins=bins)

    ##
    result      =  {}

    print('\n\n')

    for i, _bin in enumerate(bins[:-1]):
      nhalos         = np.sum(bmass == i)
      mean_mass      = np.mean(hmass[bmass == i])
        
      sample         = iscentral[bmass == i].astype(np.float)
        
      uhalos, counts = np.unique(haloindex[bmass == i], return_counts = True)
  
      label          = str(_bin / 1e10)
      result[label]  = {'mean_hmass': mean_mass / 1e10, 'cen': np.sum(sample), 'sat': np.sum(1 - sample), 'tot': len(sample), 'satfrac': 100. * np.sum(1. - sample) / len(sample),\
                        'nhalos': nhalos, 'nuhalos': len(uhalos)}
  
      print(result[label])


    masses = np.array([result[key]['mean_hmass'] for key in result.keys()])
    expcen = np.array([result[key]['cen'] / result[key]['nuhalos'] for key in result.keys()])
    
    print(masses)
    print(expcen)

    pl.semilogy(np.log10(masses), expcen, label='Centrals')

    ##pl.xlim(11., 14.5)
    ##  pl.ylim(0.1, 300.)

    pl.xlabel(r'$\log_{10} | M_h / M_\odot|$')
    pl.ylabel(r'$\langle N_g \rangle$')
    
    pl.savefig('plots/hod.pdf')
    
    print('\n\nDone.\n\n')
