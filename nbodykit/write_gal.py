import  numpy             as      np
import  astropy.io.fits   as      fits

from    readgadget        import  *
from    astropy.table     import  Table
from    snaps             import  snaps
from    get_data          import  get_caesar
from    astropy.cosmology import  FlatLambdaCDM
from    insample          import  read_insample


cosmo              =  FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.725)

def read_zpos(boxsize, redshift, set_insample=0):
    print('\n\nSolving for redshift {}.'.format(redshift))

    caesar         =  get_caesar(boxsize, redshift)

    pos            =  [list([x.GroupID]) + list(x.pos.to('Mpccm/h').value) for x in caesar.galaxies]  # comoving Mpc/h.             
    pos            =  Table(np.array(pos), names=('id', 'x', 'y', 'z'))
    
    pos['x'].unit  = 'Mpc/h'
    pos['y'].unit  = 'Mpc/h'
    pos['z'].unit  = 'Mpc/h'

    print(boxsize, np.max(pos['x']), np.max(pos['y']), np.max(pos['z']))
    
    for _ in ['x', 'y', 'z']:
      assert  np.isclose(boxsize, np.max(pos[_]), atol=1.e-1, rtol=0.0)
    
    # Velocity
    vel            =  [list(x.vel.to('km/s').value * (1. + redshift) / 100. / cosmo.efunc(redshift)) for x in caesar.galaxies]  # km/s.  
    vel            =  Table(np.array(vel), names=('vx', 'vy', 'vz'))

    vel['vx'].unit = 'Mpc/h'
    vel['vy'].unit = 'Mpc/h'
    vel['vz'].unit = 'Mpc/h'

    # vel.sort('vz')
    
    # Redshift-space position.
    zpos           = Table(pos, copy=True)
    zpos['z']      = pos['z'] + vel['vz']

    # Wrap redshift-space position.
    zpos['z']      = np.mod(zpos['z'], boxsize)
    
    print(zpos)

    if set_insample:
      insample     = read_insample(redshift)

      print(insample)

      targetids    = insample['TARGETIDS'].to_numpy()
      insample     = insample['INSAMPLE'].to_numpy()

      print(len(pos))
      print(len(targetids))
      
      assert  np.all(np.array(pos['id']).astype(np.int) == targetids.astype(np.float))

      pos['INSAMPLE']  = insample
      zpos['INSAMPLE'] = insample
      
    pos.write('../bigdat/simba_gpos_{:.5f}.fits'.format(redshift),   format='fits', overwrite=True)
    zpos.write('../bigdat/simba_gzpos_{:.5f}.fits'.format(redshift), format='fits', overwrite=True)

    return  0 


if __name__ == '__main__':
  set_insample = 1

  for x in snaps.keys():
    read_zpos(100., x, set_insample=set_insample)
    
  print('\n\nDone.\n\n')
