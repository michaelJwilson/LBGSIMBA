import  numpy             as      np
import  astropy.io.fits   as      fits

from    readgadget        import  *
from    astropy.table     import  Table
from    get_data          import  snaps, get_caesar
from    astropy.cosmology import  FlatLambdaCDM


cosmo      = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)
tracer     = 'dm'                                                                        

def write_all(boxsize=100.):
  for x in snaps.keys():
    snap   = snaps[x]    
    fpath  = '/home/rad/data/m100n1024/s50/snap_m100n1024_{}.hdf5'.format(snap)
    
    # https://bitbucket.org/rthompson/pygadgetreader/src/default/#markdown-header-readsnap
    # readsnap(snapfile,'pos','star',units=1,suppress=1)/(1+redshift)/h # pkpc
    pos    = readsnap(fpath, 'pos', tracer, units=1, nth=32, debug=True)                   # in comoving kpc/h.
    pos   /= 1.e3                                                                          # in comoving Mpc/h.   

    pos    = Table(pos, names=('x', 'y', 'z'))

    for _ in ['x', 'y', 'z']: 
      assert  np.isclose(boxsize, np.max(pos[_]), atol=1.e-4, rtol=0.0)
    
    pos['x'].unit = 'Mpc/h'
    pos['y'].unit = 'Mpc/h'
    pos['z'].unit = 'Mpc/h'
    
    print('\n\nz={}\n\n'.format(x))
    print(pos)
    
    # Velocity
    vel   =  readsnap(fpath, 'vel', tracer, units=1, suppress=1, nth=32, debug=True)       # vel. in physical km/s; peculiar?
    vel  *= (1. + x)
    vel  /=  100. * cosmo.efunc(x)                                                         # [Mpc/h].

    vel   =  Table(vel, names=('vx', 'vy', 'vz'))
    
    vel['vx'].unit = 'Mpc/h'
    vel['vy'].unit = 'Mpc/h'
    vel['vz'].unit = 'Mpc/h'

    # vel.sort('vz')

    print(vel)
    
    # Redshift-space position. 
    zpos      = Table(pos, copy=True)
    zpos['z'] = pos['z'] + vel['vz']

    # Wrap redshift-space position.
    zpos['z'] = np.mod(zpos['z'], boxsize)
    
    pos.write('../dat/simba_{}pos_{:.5f}.fits'.format(tracer, x),   format='fits', overwrite=True)
    zpos.write('../dat/simba_{}zpos_{:.5f}.fits'.format(tracer, x), format='fits', overwrite=True)


if __name__ == '__main__':
    print('\n\nWriting Simba dark matter.\n\n')
    
    pos = write_all()
    
    print('\n\nDone.\n\n')
