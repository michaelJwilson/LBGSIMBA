import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  corner
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    scipy.spatial      import  KDTree 
from    itertools          import  product
from    get_data           import  get_data
from    utils              import  latexify


matplotlib.rc('text', usetex = True)
##  plt.style.use('ggplot')

if __name__ == '__main__':
    print('\n\nWelcome to Simba cornerplot.')
    
    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
    boxsize     =  100.
    getredshift =  3.00307

    f, p        =  get_data(boxsize, getredshift)

    ##
    gid         =  f['galaxy_data']['GroupID'][:]
    iscentral   =  f['galaxy_data']['central'][:]
    haloindex   =  f['galaxy_data']['parent_halo_index'][:]
    cid         =  p['CAESAR_ID'][:]

    ##  1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'   
    LUMUV       =  p['absmag_19'][:]

    ##  Assert on the CAESAR ordering:  i.e. equating GroupID to CAESAR_ID for galaxy catalogue. 
    assert  np.all(gid == cid)
    
    ##  Basic stats.
    print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
    print('Number of centrals found: {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))
    
    ##  
    sfr         =  f['galaxy_data']['sfr'][:]
    bhmdot      =  f['galaxy_data']['bhmdot'][:]
    gfrac       =  f['galaxy_data']['gas_fraction'][:]
    
    smass       =  1.82e7                                ##  [Solar masses].
    stellarmass =  f['galaxy_data']['nstar'][:] * smass
    ssfr        =  sfr / stellarmass
    
    toplot      =  np.c_[sfr[sfr < 20.], np.log10(stellarmass[sfr < 20.])]

    ##  Quantiles=[0.16, 0.5, 0.84], hist2d_kwargs=hist2d_kwargs
    figure      = corner.corner(toplot, labels=[r"$SFR$", r"$\log_{10}(M_*)$"], show_titles=False,\
                                title_kwargs={"fontsize": 12}, use_math_text=True, plot_contours=False)
    '''
    axes        =  figure.axes

    for ax in axes: 
      ax.set_axis_on()

      ax.spines['bottom'].set_color('black')
      ax.spines['top'].set_color('black')
      ax.spines['left'].set_color('black')
      ax.spines['right'].set_color('black')
    '''
    pl.savefig('plots/cornerplot.pdf')
    
    print('\n\nDone.\n\n')
