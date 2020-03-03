import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  pandas            as   pd
import  numpy             as   np
import  pylab             as   pl
import  matplotlib.pyplot as   plt

from    scipy.spatial  import  KDTree 
from    itertools      import  product
from    get_data       import  get_data, print_keys, get_caesar
from    sutils         import  latexify
from    fithod         import  cen_model, sat_model
from    insample       import  read_insample


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def run_hod(boxsize=100., getredshift=3.00307, set_insample=0):
    ##  https://caesar.readthedocs.io/en/latest/usage.html  
    f           =  get_caesar(boxsize, getredshift, load_halo=True)
        
    no_central  =  np.array([x.central_galaxy == None for x in f.halos], dtype=np.int)

    ##  List, one for each halo. 
    Nc          =  np.array([(1-x) for x in no_central], dtype=float)
    Ns          =  np.array([np.float(len(x.satellite_galaxies)) for x in f.halos])
    
    nhalo       =  len(f.halos)

    ncentral    =  np.sum(Nc)
    nsats       =  np.sum(Ns)

    ngal        =  np.int(ncentral + nsats)

    if set_insample:
      ##  Get in sample.
      insample  =  read_insample(getredshift)
      insample  =  insample.head(n=ngal)

      print(insample)
          
      insample  =  insample['INSAMPLE']
      
    else:
      insample  =  np.ones_like(Nc)

    ##
    hmass       = [x.masses['dm'] for x in f.halos]
    hmass       =  np.array([np.atleast_1d(x.value)[0] for x in hmass])
    
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
                    'cen':   np.mean(Nc[(bmass == i) & insample]),\
                    'sat':   np.mean(Ns[(bmass == i) & insample]),\
                    'cen2':   np.var(Nc[(bmass == i) & insample]),\
                    'sat2':   np.var(Ns[(bmass == i) & insample]),\
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
    np.savetxt('dat/hod_{}_insample_{}.txt'.format(str(getredshift).replace('.', 'p'), set_insample), np.c_[masses, expcen, stdcen, expsat, stdsat], fmt='%.6le')

def plot_hod(boxsize=100., plot_model=False):
    ##                                                                                                                                                                                                                              
    latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

    ##                                                                                                                                                                                                                              
    fig, axes    = plt.subplots(nrows=1, ncols=3, sharey=False)
    
    plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)
    
    for i, getredshift in enumerate([2.024621, 3.00307, 3.963392]):
        masses, expcen, stdcen, expsat, stdsat = np.loadtxt('dat/hod_{}.txt'.format(str(getredshift).replace('.', 'p')), unpack=True)

        print(masses)

        print(expcen)
        print(expsat)

        print(stdcen)
        print(stdsat)
        
        axes[i].errorbar(np.log10(masses), expcen, stdcen, label='Centrals',   c='k',        alpha=0.8, marker='^', markersize=3, lw=0, elinewidth=1)
        axes[i].errorbar(np.log10(masses), expsat, stdsat, label='Satellites', c='darkcyan', alpha=0.8, marker='^', markersize=3, lw=0, elinewidth=1)

        axes[i].set_xticks(ticks=np.logspace(10.5, 14.5, 20))
        axes[i].set_xticklabels(labels=['{:.1f}'.format(x) for x in np.logspace(10.5, 14.5, 20)])
        
        if plot_model:
            ##  Plot best-fit models.                                                                                                                                                                                             
            ordinate  = np.logspace(10., 15., num=200)

            cenparams = np.loadtxt('dat/hod-nc-params_{}.txt'.format(str(getredshift).replace('.', 'p')))
            axes[i].semilogy(np.log10(masses), cen_model(masses, cenparams), c='k', alpha=0.8, lw=1)                                                                                                                                                                                                                                                                                                                                       
            satparams = np.loadtxt('dat/hod-ns-params_{}.txt'.format(str(getredshift).replace('.', 'p')))                                                                                                                        
            axes[i].semilogy(np.log10(masses), sat_model(masses, satparams), c='darkcyan', alpha=0.8, lw=1)
            
        ##
        axes[i].set_xscale('linear')
        axes[i].set_yscale('log')

        axes[i].set_xlim(10.5, 14.5)
        axes[i].set_ylim(0.01, 100.)

        axes[i].set_xlabel(r'$\log_{10} | M_h / M_\odot|$', labelpad=15)

    axes[0].set_ylabel(r'$\langle N_g \rangle$')
    axes[0].legend(loc=2, frameon=False)
    
    for ax in axes:
        ax.set_axis_on()

        ax.spines['bottom'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['right'].set_color('black')

        ax.tick_params(axis='x', colors='black')
        
    plt.tight_layout()
    
    pl.savefig('plots/hod.pdf')

    
if __name__ == '__main__':
    print('\n\nWelcome to Simba HOD.')

    test          =  True
    insample      =  1
    redshifts     = [2.024621, 3.00307, 3.963392, 5.0244]
    
    if test:
        boxsize   = 50.
    
    else:
        boxsize   = 100.
        
    for redshift in redshifts:
        run_hod(boxsize, getredshift=redshift, set_insample=insample)
        
        break
        
    # plot_hod(boxsize=100., plot_model=True)
        
    print('\n\nDone.\n\n')
    
