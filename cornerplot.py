import  matplotlib;  matplotlib.use('Agg')

import  glob
import  h5py 
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt
import  corner    	   as 	   corner

from    scipy.spatial      import  KDTree 
from    itertools          import  product
from    get_data           import  get_data
from    utils              import  latexify


##  matplotlib.rc('text', usetex = True)
plt.style.use('ggplot')

if __name__ == '__main__':
    print('\n\nWelcome to Simba cornerplot.')
    '''
    K        =   5
    factor   = 2.0           # size of one side of one panel                                                                                                                                                                    

    lbdim   = 0.5 * factor   # size of left/bottom margin                                                                                                                                                                         
    trdim   = 0.2 * factor   # size of top/right margin                                                                                                                                                                            

    whspace = 0.1 * 2.      # w/hspace size                                                                                                                                                                                       
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim     = lbdim + plotdim + trdim
    
    figure, axes = pl.subplots(K, K, figsize=(dim, dim))
    axes         = figure.axes
    '''
    for getredshift, reverse, color in zip([2.024621, 3.00307, 3.963392, 5.0244], [False] * 3, ['yellow', 'dodgerblue', 'green', 'brickred']):
        pl.clf()

        print('Solving for z={}.'.format(getredshift))
        
        ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
        boxsize     =  100.
        f, p        =  get_data(boxsize, getredshift)

        ##
        gid         =  f['galaxy_data']['GroupID'][:]
        iscentral   =  f['galaxy_data']['central'][:]
        haloindex   =  f['galaxy_data']['parent_halo_index'][:]
        cid         =  p['CAESAR_ID'][:]
        
        ##  1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'   
        MUV         =  p['absmag_19'][:]
        
        ##  
        lsstu       =  p['appmag_28'][:]
        lsstg       =  p['appmag_29'][:]
        lsstr       =  p['appmag_30'][:]
        lssti       =  p['appmag_31'][:]
        
        ##  Assert on the CAESAR ordering:  i.e. equating GroupID to CAESAR_ID for galaxy catalogue. 
        assert  np.all(gid == cid)
        
        ##  Basic stats.
        print('\n\nNumber of galaxies found: {}'.format(len(iscentral)))
        print('Number of centrals found: {}'.format(np.sum(iscentral)))
        print('Number of satellites found: {}'.format(np.sum(1. - iscentral)))
        
        ##  
        sfr         =  f['galaxy_data']['sfr'][:]
    
        smass       =  1.82e7                                ##  [Solar masses].
        stellarmass =  f['galaxy_data']['nstar'][:] * smass
        ssfr        =  sfr / stellarmass    

        ##
        hid         =  f['halo_data']['GroupID'][:]
        ndm         =  f['halo_data']['ndm'][:]

        ##  DM particle mass: 9.6e7 [Solar Mass]                                                                                                                                            
        pmass       =  9.6e7
        _hmass      =  ndm * pmass
        hmass       =  np.array([_hmass[hid == x] for x in haloindex])
        
        toplot      =  np.c_[np.log10(hmass[sfr > 0.]), np.log10(stellarmass[sfr > 0.]), MUV[sfr > 0.], np.log10(sfr[sfr > 0.]), lssti[sfr > 0.]]
        
        hist_kwargs = {'color': color, 'alpha': 0.5}
    
        ##  Quantiles=[0.16, 0.5, 0.84];  **hist_kwargs
        figure      = corner.corner(toplot, bins=10, labels=[r"$\log_{10}(M_h)$", r"$\log_{10}(M_*)$", r"$M_{UV}$", r"$\log_{10}$(SFR)", r"$i_{AB}$"], show_titles=False,\
                                    title_kwargs={"fontsize": 12}, use_math_text=False, plot_contours=False, smooth=None, plot_density=False, reverse=reverse, hist_kwargs=hist_kwargs, **hist_kwargs)
        
        pl.savefig('plots/cornerplot_z{}.pdf'.format(np.round(getredshift).astype(np.int)))
    
    print('\n\nDone.\n\n')
