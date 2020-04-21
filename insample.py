import matplotlib;  matplotlib.use('PDF')

import os
import numpy    as     np
import pylab    as     pl
import pandas   as     pd

from   snaps    import snaps  as snaps_c
from   snaps    import snaps_lhalf, snaps_hhalf


def insample(selection, u, g, r, i, z, y, maglim=None, default=True):
    if selection == 'u':
        isin = (u-g) > 1.5
        isin = isin & ((g-r) > -1.0) & ((g-r) < 1.2)
        isin = isin & (1.5 * (g-r) < (u-g) -0.75)

        if maglim is not None:
          if default:
            maglim = 24.6
            
          isin = isin & (i <= maglim)
        
        return  isin

    elif selection == 'g':
        isin = (g-r) >	1.0 
        isin = isin & ((r-i) <	1.0)
        isin = isin & ((g-r) >	1.5 * (r-i) + 0.8)

        if maglim is not None:
          if default:
            maglim = 25.8  

          isin = isin &	(i < maglim)
        
        return isin

    elif selection == 'r':
        isin = (i-z) >  1.5
        isin = isin & ((z-y) < 0.5)
        isin = isin & ((i-z) >  2.0 * (z-y) + 1.1)

        if maglim is not None:
          if default:
            maglim = 25.8 
            
          isin = isin &	(z < maglim)
        
        return isin

def insample_steidel(selection, Un, G, R, maglim=None, default=True):
    if selection == 'BX':
        ##  Steidel filters.
        ##  https://arxiv.org/pdf/astro-ph/0401445.pdf
        isin = (G-R)          >= -0.2
        isin = isin & ((Un-G) >= (G-R) + 0.2)
        isin = isin & ((Un-G) <  (G-R) + 1.0)
        isin = isin & (( G-R) <= 0.2 * (Un-G) + 0.4)

        if maglim is not None:
          if default:
            maglim = 25.5
            
          isin = isin & (R >= 23.5) & (R <= maglim)
        
        return isin

    elif selection == 'BM':
	##  Steidel filters; 1.4 < z < 3.3, R > 23.5
        isin = (G-R)          >= -0.2
        isin = isin & ((Un-G) >= (G-R) - 0.1)
        isin = isin & ((Un-G) <  (G-R) + 0.2)
        isin = isin & (( G-R) <= 0.2 * (Un-G) + 0.4)

        if maglim is not None:
          if default:
            maglim = 25.5
            
          isin = isin & (R >= 23.5) & (R <= maglim)
        
        return isin

    else:
        raise  ValueWarning('Specific selection ({}) is not available.'.format(selection))

def read_insample(getredshift):
    ##  Read insample and luptitudes.
    root   = os.environ['LBGSIMBA']
    fpaths = {2.024621: '/bigdat/insample_two.h5', 3.00307: '/bigdat/insample_three.h5', 3.963392: '/bigdat/insample_four.h5', 5.0244: '/bigdat/insample_five.h5'}

    fpath  = root + fpaths[getredshift] 

    print('\n\nReading in sample: {}.'.format(fpath))

    frame  = pd.read_hdf(fpath, 'df')

    return  frame

def test_insample(root=None):
    import pylab    as     pl
    from   get_data import get_pyloser

    
    nrows              =   -1
    boxsize            =  100.

    snaps              =  snaps_c
    zs                 =  np.sort(list(snaps.keys()))

    # Overwrite with custom simba run.  Gaussian EBV of magnitude 0.3 
    root               = '/home/mjwilson/LBGSIMBA/pylosers/m100n1024/s50/run/'

    if root is None:
      root             = '/home/mjwilson/LBGSIMBA/'
    
    mlims              =  np.arange(22.5, 30.0, 0.1)
    results            =  {}
    
    for nodust, name in zip([True, False], ['nodust', 'dust']):
      results[name]    = []
        
      if snaps == snaps_c:
        # wave, two,   ids2  =  get_pyloser(boxsize, zs[0], nrows=nrows, snaps=snaps, nodust=nodust, steidel=True, root=root)
        wave, three, ids3    =  get_pyloser(boxsize, zs[1], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False, root=root)
        # wave, four,  ids4  =  get_pyloser(boxsize, zs[2], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)
        # wave, five,  ids5  =  get_pyloser(boxsize, zs[3], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)

      else:
          # No dust magnitudes not available at displaced redshifts, e.g. 2.0 -> 1.5 etc. 
          # wave, two,   ids2 =  get_pyloser(boxsize, zs[0], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False, root=root)
          wave, three, ids3   =  get_pyloser(boxsize, zs[1], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False, root=root)
          # wave, four,  ids4 =  get_pyloser(boxsize, zs[2], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)
	  # wave, five,  ids5 =  get_pyloser(boxsize, zs[3], nrows=nrows, snaps=snaps, nodust=nodust, steidel=False)
          
      for maglim in mlims:
        '''
        if snaps == snaps_c:
          bxdrops            =  insample_steidel('BX', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=maglim, default=False)
          bmdrops            =  insample_steidel('BM', two['steidel_un'].values, two['steidel_g'].values, two['steidel_rs'].values, maglim=maglim, default=False)

        else:
          bxdrops            =  np.zeros_like(two['LSST_u'].values).astype(bool)
          bmdrops            =  np.zeros_like(two['LSST_u'].values).astype(bool)
        '''
        
        ##  
        udrops               =  insample('u',  three['LSST_u'].values, three['LSST_g'].values, three['LSST_r'].values, three['LSST_i'].values, three['LSST_z'].values, three['LSST_y'].values, maglim=maglim, default=False)
        # gdrops             =  insample('g',   four['LSST_u'].values,  four['LSST_g'].values,  four['LSST_r'].values,  four['LSST_i'].values,  four['LSST_z'].values,  four['LSST_y'].values, maglim=maglim, default=False)
        # rdrops             =  insample('r',   five['LSST_u'].values,  five['LSST_g'].values,  five['LSST_r'].values,  five['LSST_i'].values,  five['LSST_z'].values,  five['LSST_y'].values, maglim=maglim, default=False)

        # Color selected.
        tmp                  =  [np.count_nonzero(udrops)]
        # tmp                =  [np.count_nonzero(bxdrops), np.count_nonzero(bmdrops), np.count_nonzero(udrops), np.count_nonzero(gdrops), np.count_nonzero(rdrops)]

        '''
        # Pure magnitude limits.
        tmp2                 =  [np.count_nonzero((two['steidel_rs'] >= 23.5) & (two['steidel_rs'] <= maglim)),\
                                 np.count_nonzero(three['LSST_i'] < maglim), np.count_nonzero(four['LSST_i'] < maglim),\
                                 np.count_nonzero(five['LSST_z'] < maglim)]
        '''

        tmp2                 =  [np.count_nonzero(three['LSST_i'] < maglim)]
        
        results[name].append(tmp + tmp2)
    
    results['nodust'] = pd.DataFrame(np.c_[mlims, np.array(results['nodust'])], columns=['mlim', 'u', 'three'])
    results['dust']   = pd.DataFrame(np.c_[mlims, np.array(results['dust'])],   columns=['mlim', 'u', 'three'])
                                 
    # results['nodust'] = pd.DataFrame(np.c_[mlims, np.array(results['nodust'])], columns=['mlim', 'BX', 'BM', 'u', 'g', 'r', 'two', 'three', 'four', 'five'])
    # results['dust']   = pd.DataFrame(np.c_[mlims, np.array(results['dust'])],   columns=['mlim', 'BX', 'BM', 'u', 'g', 'r', 'two', 'three', 'four', 'five'])

    print('\n\n ----  No dust  ----')
    print(results['nodust'])
    print('\n\n ----  Dust  ----')
    print(results['dust'])                                      

    results['nodust'].to_csv(root + 'no_dust_Ns_{}.csv'.format(boxsize))
    results['dust'].to_csv(root + 'dust_Ns_{}.csv'.format(boxsize))
    
    
if __name__ == '__main__':
    snaps = {2.024621: '078', 3.00307: '062', 3.963392: '051', 5.0244: '042'}
    '''
    for key in snaps.keys():
      frame = read_insample(key)

      print(frame)
    '''

    test_insample()
    
    print('\n\nDone.\n\n')

