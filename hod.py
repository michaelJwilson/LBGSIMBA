import  glob
import  h5py
import  numpy         as      np
import  pylab         as      pl

from    scipy.spatial import  KDTree 
from    itertools     import  product


def print_keys(arg):
    print('\n\nAvailable keys for {}:\n{}'.format(arg, arg.keys()))

def get_data(boxsize, getredshift):
    ##  Find the snapshot file closest to the desired redshift, getredshift. 
    files      = glob.glob('/home/mjwilson/sim/m100n1024/s50/Groups/m100n1024_*.hdf5')
    snapshots  = [int(x.split('/')[-1].split('_')[-1].split('.hdf5')[0]) for x in files]
    sortind    = np.argsort(snapshots)
    sortfiles  = [files[i] for i in sortind]
    snapshots  = [int(x.split('/')[-1].split('_')[-1].split('.hdf5')[0]) for x in sortfiles]

    redshifts  = np.loadtxt('/home/mjwilson/snapshot_redshifts.txt', dtype={'names': ['a', 'z', 'id'], 'formats': ['float', 'float', 'int']})
    redshifts  = np.array([x[1] for x in redshifts[snapshots]])

    zdiff      = np.abs(redshifts - 3.00307)
    index      = np.where(zdiff == zdiff.min())[0]

    getfile    = sortfiles[index[0]]

    parts      = getfile.split('/')
    file       = parts[-1]
    rpath      = '/'.join(parts[:-1])
    photfile   = rpath + '/phot_' + file

    print('\n\nGetting snapshot {} for redshift {}.'.format(snapshots[index[0]], getredshift))
    print('\n\nRetrieving:  {}'.format(getfile))
    print('\n\nRetrieving:  {}'.format(photfile))

    ##  Extract desired data from this redshift snapshot. 
    f         = h5py.File(getfile)
    p         = h5py.File(photfile)

    print_keys(f)
    print_keys(f['galaxy_data'])
    print_keys(f['halo_data'])
    
    print_keys(p)

    return f, p 


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

    ##  print(p['COLOR_INFO'][:])

    ##
    cid         =  p['CAESAR_ID'][:]
    LUMUV       =  p['absmag_19'][:]  ## 1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'

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
    
    print('\n\nRange of stellar mass: {} to {} solar masses.'.format(stellarmass.min() / 1e10, stellarmass.max() / 1e10))

    ##  Stellar mass function.
    bins          =  np.logspace(9., 13., 10, endpoint=True, base=10.0)
    bsmass        =  np.digitize(stellarmass, bins=bins)

    ninds, counts =  np.unique(bsmass, return_counts = True)
    mean_smass    =	 np.array([np.mean(stellarmass[bsmass == _bin]) for _bin in np.arange(len(bins))])

    ##  Check that all bins are populated by at least one galaxy. 
    assert  len(ninds) == (len(bins) - 1)

    ##  print(mean_smass)
    ##  print(counts)

    ##  Limit for star-forming to quenched segregation based on specific star formation rate.  Pg. 11 of 1901.10203
    qlimit      = 10. ** (-1.8 + 0.3 * getredshift)      ##  [Per Gyr].

    ##  print('\n\nNumber of star-forming galaxies found: {}'.format(np.sum(ssfr > qlimit)))
    ##  print('Number of quenched galaxies found: {}'.format(np.sum(ssfr < qlimit)))
    
    ##  Positions in kpc.
    pos         = f['galaxy_data']['pos'][:]
    pos        /= 1.e3
    
    ##  Test.
    ##  pos     = pos[:50]
    
    ##  Apply periodic reflections.
    catalogue   = np.copy(pos)

    reflections = list(product('+0-', repeat=3))
    reflections.remove(tuple(['0', '0', '0']))

    convert     = {'+': +1., '0': 0., '-': -1.}

    for reflection in reflections:
        _copy       = np.copy(pos)
        signs       = [convert[x] for x in reflection]
    
        _copy[:,0] += signs[0] * boxsize
        _copy[:,1] += signs[1] * boxsize
        _copy[:,2] += signs[2] * boxsize

        catalogue   = np.vstack([catalogue, _copy])
        
    print(len(pos), 27 * len(pos), len(catalogue))
    print(catalogue.shape)

    print(catalogue)

    ##  Grow tree. 
    PTree       =  KDTree(pos)
    CTree       =  KDTree(catalogue)

    ##  Note:  asymmetric catalogue-based Tree call and pos call - i.e. do not count pairs between
    ##         two reflections. 
    paired      =  CTree.query_ball_tree(PTree, 25.)

    dr          =  2.5
    bins        =  np.arange(0.0, 30.0, dr)

    sep         =  []

    for i, row in enumerate(catalogue):
        for j, twin in enumerate(paired[i]):
            sep.append(np.sum((row - catalogue[twin])**2.))

    sep         =  np.sqrt(np.array(sep))
    bsep        =  np.digitize(sep, bins=bins)

    meanr       =  np.array([np.mean(sep[bsep == x]) for x in range(len(bins) - 1)]) 
    cnts        =  np.array([np.sum(bsep == x) for x in range(len(bins) - 1)])

    ngal        =  catalogue.shape[0]
    vol         =  boxsize ** 3.
    nbar        =  ngal / vol

    rr          =  (4. * np.pi / 3.) * nbar * nbar * vol * (( bins[:-1] + dr / 2. ) ** 3. - ( bins[:-1] - dr / 2. ) ** 3.)
    
    xi          =  cnts / rr - 1.

    print(bins + dr / 2.) 
    print(meanr)
    print(cnts)
    print(rr)
    print(xi)

    ##
    hid         =  f['halo_data']['GroupID'][:]
    ndm         =  f['halo_data']['ndm'][:]

    ##  allthere    =   any([ind in hid for ind in haloindex])
    ##  print('\n\nAll parent halos (to galaxies) present?  {}'.format(allthere))

    ##  DM particle mass: 9.6e7 [Solar Mass]
    pmass       =  9.6e7
    hmass       =  ndm[haloindex] * pmass

    ##  Basic halo stats.                                                                                                                                                                                           
    print('\n\nNumber of halos found: {}M'.format(len(hid) / 1e6))
    print('Range of halo masses: {} to {} [10^9]'.format(hmass.min() / 1e9, hmass.max() / 1e9))

    ##  Now bin galaxies and halos by halo mass.  
    bins        =  np.logspace(10., 14., 10, endpoint=True, base=10.0)
    bmass       =  np.digitize(hmass, bins=bins)

    ##
    result      =  {}

    print('\n\n')

    for i, _bin in enumerate(bins):
        nhalos         = np.sum(bmass == i)
        mean_mass      = np.mean(hmass[bmass == i])
        
        sample         = iscentral[bmass == i]
        
        uhalos, counts = np.unique(haloindex[bmass == i], return_counts = True)
  
        label          = str(_bin / 1e10)
        result[label]  = {'mean_hmass': mean_mass / 1e10, 'cen': np.sum(sample), 'sat': np.sum(1 - sample), 'tot': len(sample), 'satfrac': 100. * np.sum(1. - sample) / len(sample),\
                          'nhalos': nhalos, 'uhalos': len(uhalos)}
  
        print(result[label])

    print('\n\nDone.\n\n')
