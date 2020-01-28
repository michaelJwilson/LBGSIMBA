import  numpy             as      np
import  astropy.io.fits   as      fits

from    readgadget        import  *
from    astropy.table     import  Table
from    get_data          import  snaps 
from    astropy.cosmology import  FlatLambdaCDM


cosmo      = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)

for x in snaps.keys():
    snap   = snaps[x]    
    fpath  = '/home/rad/data/m100n1024/s50/snap_m100n1024_{}.hdf5'.format(snap)
    
    # https://bitbucket.org/rthompson/pygadgetreader/src/default/#markdown-header-readsnap
    # readsnap(snapfile,'pos','star',units=1,suppress=1)/(1+redshift)/h # pkpc
    dmpos  = readsnap(fpath, 'pos', 'dm', units=1, nth=32, debug=True)                # in comoving kpc/h.
    dmpos /= 1.e3                                                                     # in comoving Mpc/h.   

    dmpos  = Table(dmpos, names=('x', 'y', 'z'))

    dmpos['x'].unit = 'Mpc/h'
    dmpos['y'].unit = 'Mpc/h'
    dmpos['z'].unit = 'Mpc/h'
    
    print('\n\nz={}\n\n'.format(x))
    print(dmpos)
    
    # Velocity
    dmvel   =  readsnap(fpath, 'vel', 'dm', units=1, suppress=1, nth=32, debug=True)  # vel. in physical km/s; peculiar?
    dmvel  *= (1. + x)
    dmvel  /=  100. * cosmo.efunc(x)                                                  # [Mpc/h].

    dmvel   = Table(dmvel, names=('vx', 'vy', 'vz'))
    
    dmvel['vx'].unit = 'Mpc/h'
    dmvel['vy'].unit = 'Mpc/h'
    dmvel['vz'].unit = 'Mpc/h'

    # dmvel.sort('vz')

    print(dmvel)
    
    # Redshift-space position. 
    dmzpos      = Table(dmpos, copy=True)
    dmzpos['z'] = dmpos['z'] + dmvel['vz']
    
    dmpos.write('../dat/simba_dmpos_{:.5f}.fits'.format(x),   format='fits', overwrite=True)
    dmzpos.write('../dat/simba_dmzpos_{:.5f}.fits'.format(x), format='fits', overwrite=True)
    
print('\n\nDone.\n\n')
