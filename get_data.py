import  glob
import  h5py
import  numpy           as      np
import  pylab           as      pl
import  astropy.io.fits as      fits

from    scipy.spatial   import  KDTree 
from    itertools       import  product
from    astropy.table   import  Table


def print_keys(arg):
    print('\n\nAvailable keys for {}:\n{}'.format(arg, arg.keys()))

def get_data(boxsize, getredshift):
    ##  Find the snapshot file closest to the desired redshift, getredshift. 
    files      = glob.glob('/home/mjwilson/LBGSIMBA/100MPC/m100n1024_*.hdf5')
    snapshots  = [int(x.split('/')[-1].split('_')[-1].split('.hdf5')[0]) for x in files]
    sortind    = np.argsort(snapshots)
    sortfiles  = [files[i] for i in sortind]
    snapshots  = [int(x.split('/')[-1].split('_')[-1].split('.hdf5')[0]) for x in sortfiles]

    redshifts  = np.loadtxt('/home/mjwilson/LBGSIMBA/snapshot_redshifts.txt', dtype={'names': ['a', 'z', 'id'], 'formats': ['float', 'float', 'int']})
    redshifts  = np.array([x[1] for x in redshifts[snapshots]])

    zdiff      = np.abs(redshifts - getredshift)
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
    f          = h5py.File(getfile)
    p          = h5py.File(photfile)

    print_keys(f)
    print_keys(f['galaxy_data'])
    print_keys(f['halo_data'])
    
    print_keys(p)

    return f, p 


if __name__ == '__main__':
    print('\n\nWelcome to Simba get_data.')

    boxsize     =  50.

    ##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
    f2, p2      =  get_data(boxsize, 2.024621)
    f3, p3      =  get_data(boxsize, 3.00307)
    f4, p4      =  get_data(boxsize, 3.963392)
    f5, p5      =  get_data(boxsize, 5.024400)

    print
    print(p2['HEADER_INFO'][:])
    print
    print(p3['HEADER_INFO'][:])
    print
    print(p4['HEADER_INFO'][:])
    print
    print(p5['HEADER_INFO'][:])
    print
    print(p2['COLOR_INFO'][:])
    
    print('\n\nDone.\n\n')
