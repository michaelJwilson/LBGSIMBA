import  matplotlib;  matplotlib.use('PDF')

import  copy
import  glob
import  h5py
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    scipy.spatial      import  KDTree 
from    itertools          import  product
from    get_data           import  get_data, print_keys
from    utils              import  latexify
from    fithod             import  cen_model, sat_model


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

def plot_tree(shuffled, seed):
  for key in shuffled.keys():
    halo_pos  = shuffled[key][0][0]

    cen_pos   = halo_pos + shuffled[key][1][0]
    
    for sat_pos in shuffled[key][2]:
      sat_pos = halo_pos + sat_pos[0]

      pl.plot(sat_pos[1], sat_pos[2], 'c.', markersize=7)

    pl.plot( cen_pos[1],  cen_pos[2], 'k.', markersize=3)
    pl.plot(halo_pos[1], halo_pos[2], 'w.', markersize=3, alpha=0.3)
    
  pl.savefig('plots/shuffled_{}.pdf'.format(seed))


if __name__ == '__main__':
    print('\n\nWelcome to Simba shuffled HOD.')

    ##  Greater than zero for shuffling. 
    seed           =  1

    np.random.seed(seed)
    
    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
    boxsize        =  100.
    getredshift    =  3.00307

    f, p           =  get_data(boxsize, getredshift)
    
    ##  ID for each galaxy. 

    _gid           =  f['galaxy_data']['GroupID'][:]

    ##  Is this galaxy a central. 
    _iscentral     =  f['galaxy_data']['central'][:]
    
    ##  Positions in kpc.
    _pos           =  f['galaxy_data']['pos'][:]
    _pos          /=  1.e3  ##  [Mpc].                                                                                                                                                                                                                                                                                                      
    _pos          *=  0.68  ##  [Mpc/h].
    
    ##  What is the corresponding halo ID. 
    _haloindex     =  f['galaxy_data']['parent_halo_index'][:]

    ##  Keep those galaxies within 5 Mpc/h of sim. wall.
    lim            =                        1.0
 
    isin           =            _pos[:,0] < lim
    pos            =       _pos[_pos[:,0] < lim]
    
    gid            =       _gid[_pos[:,0] < lim]
    iscentral      = _iscentral[_pos[:,0] < lim]
    haloindex      = _haloindex[_pos[:,0] < lim]
    
    ##
    hid            =  f['halo_data']['GroupID'][:]

    ##  Number of dark matter particles in this halo. 
    ndm            =  f['halo_data']['ndm'][:]

    ##  DM particle mass: 9.6e7 [Solar Mass]                                                                                                                                                                                        
    pmass          =  9.6e7
    hmass          =  ndm * pmass

    ##
    hpos           =  f['halo_data']['pos'][:]
    hpos          /=  1.e3  ##  [Mpc].
    hpos          *=  0.68  ##  [Mpc/h].        

    ##  Basic galaxy stats.
    print('\n\n')
    print('Number of galaxies found:   {}'.format(len(iscentral)))
    print('Number of centrals found:   {}'.format(np.sum(iscentral)))
    print('Number of satellites found: {}'.format(np.sum(1 - iscentral)))

    ##  Basic halo stats.                                                                                                                                                                                           
    print('\n\nNumber of halos found: {}M.'.format(len(hid) / 1e6))
    print('Range of halo masses: {} [10^9] to {} [10^13].'.format(hmass.min() / 1e9, hmass.max() / 1e13))

    ##  Now bin galaxies and halos by halo mass.  
    bins        =  np.logspace(10., 14., 10, endpoint=True, base=10.0)

    ##  Binned mass. 
    bmass       =  np.digitize(hmass, bins=bins)
    
    ##
    result      =  {}

    for i, _bin in enumerate(bins):
      ##  Number of haloes in this mass bin.   
      nhalos       = np.sum(bmass == i)

      ##  Mean halo mass of this sample.
      mean_mass    = np.mean(hmass[bmass == i])

      ##  Halo IDs that make the sample. 
      hsample      = hid[bmass == i]

      uhalos, cnts = np.unique(hsample, return_counts=True)

      ##  Each halo is included only once. 
      assert np.all(cnts == 1)

      ##  All the galaxies in this mass range. 
      ids_in      = gid[[x in uhalos for x in haloindex]] 

      ##  IDs of centrals in this halo mass range. 
      cids_in     = ids_in[iscentral[[x in uhalos for x in haloindex]]]

      ##  Sort these IDs. 
      cids_in     = np.sort(cids_in)
      
      ##  Loop through haloes in this mass range, to pair centrals to satellites based on the same halo (ID).
      for u in uhalos:
        ##  Of the galaxies in this halo, which is the central. 
        _central  = iscentral[haloindex == u]

        ##  ID of the central.
        _cid      = gid[haloindex == u][_central]
        
        if len(_cid) > 0:
            ##  One central per halo.
            assert (len(_cid) == 1) & (_cid in cids_in)

            _cid      = _cid[0]

            ##  Position of parent halo.
            _hpos     = hpos[hid == u]
            
            ##  Position of central.                                                                                                                                                                                                    
            _cpos     = pos[haloindex == u][_central]
            
            ##  IDs of the satellites.
            _sid      = gid[haloindex == u][~_central]

            ##  Save Position of satellites relative to their central.
            result[_cid] = []

            result[_cid].append(_hpos)
            result[_cid].append(_cpos - _hpos)
            result[_cid].append(pos[haloindex == u][~_central] - _hpos)
      
    ##  Sorted keys.
    cids_srtd   = np.sort(result.keys())
    
    ##  Shuffled version of the central IDs.
    cids_in_shf = np.copy(cids_srtd)

    shuffled    = copy.copy(result)

    processed   = []
    
    if seed == 0:    
      ##  No shuffling.
      pass
      
    else:  
      np.random.shuffle(cids_in_shf)  
      
      for i, cid in enumerate(cids_srtd):
        if cid not in processed:
          shuffled[cid][1]            = result[cids_in_shf[i]][1]
          shuffled[cids_in_shf[i]][1] = result[cid][1]
          
          shuffled[cid][2]            = result[cids_in_shf[i]][2]
          shuffled[cids_in_shf[i]][2] = result[cid][2]

          processed.append(cid)
          processed.append(cids_in_shf[i])
                        
    plot_tree(shuffled, seed)
                
    print('\n\nDone.\n\n')
