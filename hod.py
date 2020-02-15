import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as   np
import  pylab             as   pl
import  matplotlib.pyplot as   plt

from    scipy.spatial  import  KDTree 
from    itertools      import  product
from    get_data       import  get_data, print_keys, get_caesar
from    utils          import  latexify
from    fithod         import  cen_model, sat_model


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def run_hod(boxsize=100., getredshift=3.00307):
    ##  https://caesar.readthedocs.io/en/latest/usage.html  

    f           =  get_caesar(boxsize, getredshift, load_halo=True)

    no_central  = np.array([x.central_galaxy == None for x in links.halos], dtype=np.int)

    ##  List, one for each halo. 
    Nc          = np.array([(1-x) for x in no_central], dtype=float)
    Ns          = np.array([np.float(len(x.satellite_galaxies)) for x in links.halos])
    
    nhalo       =  len(links.halos)
    ncentral    =  np.count_nonzero(central)
    nsats       =  np.sum(Ns)

    hmass       = [x.masses['dm'] for x in links.halos]
    hmass       = np.array([np.atleast_1d(x.value)[0] for x in hmass])
    
    ##  Now bin galaxies and halos by halo mass.  
    bins        =  np.logspace(10., 14., 10, endpoint=True, base=10.0)
    bmass       =  np.digitize(hmass, bins=bins)
    
    ##
    result      =  {}

    for i, _bin in enumerate(bins):
      print('\n\n')  

      ##  Number of haloes in this mass bin.   
      nhalos     = np.count_nonzero(bmass == i)

      ##  Mean halo mass of this sample.
      mean_mass  = np.mean(hmass[bmass == i])

      result[i]  = {'mean_hmass': mean_mass,\
                    'cen':   np.mean(Nc[bmass == i]), 'sat':   np.mean(Ns[bmass == i]),\
                    'cen2':   np.var(Nc[bmass == i]), 'sat2':   np.var(Ns[bmass == i]),\
                    'nhalo':  nhalos}  

      for x in result.keys():
          print(x, result[x])

    ##
    masses  = np.array([result[key]['mean_hmass'] for key in range(len(bins))])

    expcen  = np.array([result[key]['cen'] for key in range(len(bins))])
    expsat  = np.array([result[key]['sat'] for key in range(len(bins))])

    stdcen  = np.array([np.sqrt(result[key]['cen2']) for key in range(len(bins))])
    stdsat  = np.array([np.sqrt(result[key]['sat2']) for key in range(len(bins))])
    
    print(masses)
    
    print(expcen)
    print(expsat)
    
    print(stdcen)
    print(stdsat)

    ##  Write results.
    np.savetxt('dat/hod_{}.txt'.format(str(getredshift).replace('.', 'p')), np.c_[masses, expcen, stdcen, expsat, stdsat], fmt='%.6le')

def plot_hod(boxsize=100., getredshift=3.00307, plot_model=False):
    masses, expcen, stdcen, expsat, stdsat = np.loadtxt('dat/hod_{}.txt'.format(str(getredshift).replace('.', 'p')), unpack=True)

    pl.errorbar(np.log10(masses), expcen, stdcen, label='Centrals',   c='k',        alpha=0.8, lw=0, marker='^', markersize=2)
    pl.errorbar(np.log10(masses), expsat, stdsat, label='Satellites', c='darkcyan', alpha=0.8, lw=0, marker='^', markersize=2)

    if plot_model:
      ##  Plot best-fit models.                                                                                                                                                                                                   
      ordinate  = np.logspace(10., 15., num=200)                                                                                                                                                                                                                                                                                                                                                                                                                       
      cenparams = np.loadtxt('dat/hod-nc-params.txt')                                                                                                                                                                        
      pl.semilogy(np.log10(ordinate), cen_model(ordinate, cenparams), c='k', alpha=0.8, lw=1)                                                                                                                                                                                                                                                                                                                                                               
      satparams = np.loadtxt('dat/hod-ns-params.txt')                                                                                                                                                                         
      pl.semilogy(np.log10(ordinate), sat_model(ordinate, satparams), c='darkcyan', alpha=0.8, lw=1)
      
    ## 
    pl.xscale('log')
    pl.yscale('log')

    pl.xlim(10.5, 14.5)
    pl.ylim(0.01, 100.)

    pl.xlabel(r'$\log_{10} | M_h / M_\odot|$')
    pl.ylabel(r'$\langle N_g \rangle$')

    pl.legend(loc=2, frameon=False)

    plt.tight_layout()

    pl.savefig('plots/hod_{:.3f}.pdf'.format(getredshift).replace('.', 'p'))

    
if __name__ == '__main__':
    print('\n\nWelcome to Simba HOD.')

    redshifts = [2.024621, 3.00307, 3.963392, 5.0244]

    for redshift in redshifts:
        run_hod(100., getredshift=redshift)

        plot_hod(100., getredshift=redshift)
        
    print('\n\nDone.\n\n')
    
