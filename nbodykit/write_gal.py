import  numpy             as      np
import  astropy.io.fits   as      fits

from    readgadget        import  *
from    astropy.table     import  Table
from    get_data          import  snaps, get_caesar
from    astropy.cosmology import  FlatLambdaCDM


cosmo = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)

def read_zpos(boxsize, redshift):
    caesar         =  get_caesar(boxsize, redshift)

    pos            =  [list(x.pos.to('Mpccm/h').value) for x in caesar.galaxies]  # comoving Mpc/h.         
    
    pos            =  Table(np.array(pos), names=('x', 'y', 'z'))

    pos['x'].unit  = 'Mpc/h'
    pos['y'].unit  = 'Mpc/h'
    pos['z'].unit  = 'Mpc/h'
        
    # Velocity
    vel            =  [list(x.vel.to('km/s').value * (1. + redshift) / 100. / cosmo.efunc(redshift)) for x in caesar.galaxies]  # km/s.  
    vel            =  Table(np.array(vel), names=('vx', 'vy', 'vz'))

    vel['vx'].unit = 'Mpc/h'
    vel['vy'].unit = 'Mpc/h'
    vel['vz'].unit = 'Mpc/h'

    # vel.sort('vz')                                                                                                                                                                                                                                                                        
    print(vel)

    # Redshift-space position.
    zpos           = Table(pos, copy=True)
    zpos['z']      = pos['z'] + vel['vz']

    # Wrap redshift-space position.                                                                                                                                                                                                    
    zpos['z']      = np.mod(zpos['z'], boxsize)
    
    # print(vel)

    pos.write('../dat/simba_gpos_{:.5f}.fits'.format(redshift),   format='fits', overwrite=True)
    zpos.write('../dat/simba_gzpos_{:.5f}.fits'.format(redshift), format='fits', overwrite=True)

    return  0 


if __name__ == '__main__':
  for x in snaps.keys():
    read_zpos(100., x)
    
  print('\n\nDone.\n\n')
