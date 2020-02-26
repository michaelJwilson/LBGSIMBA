import  matplotlib;  matplotlib.use('PDF')

import  copy
import  glob
import  h5py
import  fitsio
import  numpy               as      np
import  pylab               as      pl
import  matplotlib.pyplot   as      plt

from    snaps               import  snaps
from    scipy.spatial       import  KDTree 
from    itertools           import  product
## from    utils               import  latexify
from    fithod              import  cen_model, sat_model
from    caesar.quick_loader import  quick_load

## latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

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


def gen_shuffled(seed, boxsize, getredshift):    
    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
    snap           =  snaps[getredshift]
    
    cc             =  quick_load('/home/rad/data/m50n1024/s50/Groups/m50n1024_026.hdf5')

    halos          =  cc.halos

    ##  Now bin galaxies and halos by halo mass.                                                                                                                                                                    
    bins           =  np.logspace(10., 14., 10, endpoint=True, base=10.0)

    ## 
    tree           =  {}
    newtree        =  {}
    
    for i in np.arange(len(bins)):
      tree[i]      = []
      newtree[i]   = []

    ##  Conservation of sats and centrals with shuffling. 
    _ncentral      = 0
    _nsatellite    = 0 
    
    for hh in halos:
      hmass        = hh.masses['dm']                                          # Solar. ['baryon', 'dm', 'gas', 'stellar', 'total']
      digmass      = np.digitize(hmass, bins=bins)
      
      central      = hh.central_galaxy
      satellites   = hh.satellite_galaxies     

      hpos         = hh.pos.to('Mpccm/h').value                               # comoving Mpc/h. 

      if central is not None:
        cpos         = central.pos.to('Mpccm/h').value - hpos                 # comoving Mpc/h. 
        _ncentral   += 1

      else:
        cpos         = None
        
      if len(satellites) > 0:
        spos         = [x.pos.to('Mpccm/h').value -hpos for x in satellites]    # comoving Mpc/h. 
        _nsatellite += len(satellites)
        
      else:
        spos         = []

      if central is not None:
        tree[digmass].append([hmass, hpos, cpos, spos])

    ##  Seed the shuffled tree.    
    np.random.seed(seed=seed)
        
    for i in np.arange(len(bins)):
      print('\n\n---------------------------------------------------------------\n\n')

      for x in tree[i]:
        print(x[0], x[1], x[2])
      
      ncentral     = len(tree[i])

      order        = np.arange(ncentral)  

      np.random.shuffle(order)
      
      for j in np.arange(ncentral):
        entry      = tree[i][j] 
        swap       = tree[i][order[j]]
        
        # Append new halos masses. 
        newtree[i].append([entry[0], entry[1], swap[2], swap[3]])
      
      print('\n\n')

      for x in newtree[i]:
        print(x[0], x[1], x[2])

    ##
    centrals   = []
    satellites = []

    for i in np.arange(len(bins)):
       ncentral      = len(tree[i])    

       for j in np.arange(ncentral):
         entry       = tree[i][j]

         centrals.append(list(entry[1] + entry[2]))
         
         for sat in entry[3]:
           satellites.append(list(entry[1] + sat))

    centrals   = np.array(centrals) 
    satellites = np.array(satellites)

    assert  len(centrals)   == _ncentral
    assert  len(satellites) == _nsatellite

    result     = np.vstack((centrals, satellites))
    
    ##  Write shuffled catalogue.                                                                                                                                                                                                 
    fpath      = '/home/mjwilson/LBGSIMBA/bigdat/simba_gpos_{}_shuffled{:d}.fits'.format(getredshift, seed)
        
    names      = ['x', 'y', 'z']
    towrite    = {'x': result[:,0], 'y': result[:,1], 'z': result[:,2]}
        
    fitsio.write(fpath, data=towrite, names=names, clobber=True)
    

if __name__ == '__main__':
    print('\n\nWelcome to Simba shuffled HOD.')

    ##  Greater than zero for shuffling.                                                                                                                                                                                         
    seed    = 42

    np.random.seed(seed)

    boxsize = 100.

    zs      = snaps.keys()
    
    for getredshift in zs:
      print('\n\nSolving for redshift: {}.'.format(getredshift))

      gen_shuffled(seed, boxsize, getredshift)
    
    print('\n\nDone.\n\n')
