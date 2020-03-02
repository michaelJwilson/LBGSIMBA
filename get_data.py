import  glob
import  h5py
import  caesar
import  pandas          as      pd
import  numpy           as      np
import  astropy.io.fits as      fits

from    scipy.spatial   import  KDTree 
from    itertools       import  product
from    astropy.table   import  Table
from    hildebrandt     import  ferr
from    snaps           import  snaps


def print_keys(arg):
    print('\n\nAvailable keys for {}:\n{}'.format(arg, arg.keys()))

def get_data(boxsize, getredshift, printit=False):
    ##  Find the CAESAR snapshot file closest to the desired redshift, getredshift. 
    files      = glob.glob('/home/mjwilson/LBGSIMBA/100/m100n1024_*.hdf5')
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

    ##  DEPRECATED:  use pyloser instead of phot. 
    photfile   = rpath + '/phot_' + file

    print('\n\nGetting snapshot {} for redshift {}.'.format(snapshots[index[0]], getredshift))
    print('\n\nRetrieving:  {}'.format(getfile))
    print('\n\nRetrieving:  {}'.format(photfile))

    ##  Extract desired data from this redshift snapshot. 
    f          = h5py.File(getfile)
    _          = h5py.File(photfile)

    if printit:
        print_keys(f)
        print_keys(f['galaxy_data'])
        print_keys(f['halo_data'])
    
        print_keys(_)

    return f, _

def get_phys(boxsize, getredshift, printit=False, nrows=-1):
    f, p                     = get_data(boxsize, getredshift, printit=printit)

    to_return                = {}
    
    ##                                                                                                                                                                                                                              
    to_return['gid']         =  f['galaxy_data']['GroupID'][:]
    to_return['iscentral']   =  f['galaxy_data']['central'][:]
    to_return['haloindex']   =  f['galaxy_data']['parent_halo_index'][:]
    to_return['cid']         =  p['CAESAR_ID'][:]
    
    ##  1300 1700 1510 Idealized 1500A bandpass: rounded tophat centered'                                                                                                                                                           
    to_return['MUV']         =  p['absmag_19'][:]

    ## Ignore photometry that's not in the files.
    ## to_return['lsstu']    =  p['appmag_28'][:]
    ## to_return['lsstg']    =  p['appmag_29'][:]
    ## to_return['lsstr']    =  p['appmag_30'][:]
    ## to_return['lssti']    =  p['appmag_31'][:]

    ##                                                                                                                                                                                                                              
    to_return['sfr']         =  f['galaxy_data']['sfr'][:] * 1.e9     ##  [Solar masses / Gyr].

    to_return['smass_res']   =  1.82e7                                ##  [Solar masses].

    to_return['stellarmass'] =  f['galaxy_data']['nstar'][:] * to_return['smass_res']
    to_return['ssfr']        =  to_return['sfr'] / to_return['stellarmass'] 
    
    ##                                                                                                                                                                                                         
    to_return['hid']         =  f['halo_data']['GroupID'][:]
    to_return['ndm']         =  f['halo_data']['ndm'][:]

    ##  DM particle mass: 9.6e7 [Solar Mass]                                                                                                                                                                                        
    to_return['dm_pmass']    =  9.6e7
    to_return['_hmass']      =  to_return['ndm'] * to_return['dm_pmass']
    to_return['hmass']       =  np.array([to_return['_hmass'][to_return['hid'] == x][0] for x in to_return['haloindex']])

    for key in to_return.keys():
      if key not in ['smass_res', 'dm_pmass']:    
        to_return[key]       =  to_return[key][:nrows] 
    
    return  to_return  
    

def get_caesar(boxsize, redshift, load_halo=False):
    #  https://caesar.readthedocs.io/en/latest/usage.html
    #
    #  0.330620   2.024621  078
    #  0.249808   3.003070  062
    #  0.201475   3.963392  051
    #  0.165992   5.024400  042
    
    root  = '/home/rad/data/m{}n1024/s50/Groups/'.format(np.int(np.floor(boxsize)))
    
    snap  = snaps[redshift]
    fpath = root + 'm{}n1024_{}.hdf5'.format(np.int(np.floor(boxsize)), snap)

    return  caesar.load(fpath, LoadHalo=np.int(load_halo))

def get_pyloser(boxsize, redshift, printit=False, nrows=-1, magtype='app'):
    #  Currently, photometry only. 
    root   = '/home/mjwilson/LBGSIMBA/100/'
    snap   = snaps[redshift]

    fpath  = root + 'pyloser_m100n1024_{}.hdf5'.format(snap)
    
    links  = h5py.File(fpath, 'r')
    
    attrs  = links.attrs.items()

    bands  = list(links.attrs['bands'])

    for x in bands:
        print(x)
    
    frame  = pd.DataFrame(data=links['{}mag'.format(magtype)][:], columns=bands) 
    
    retain = ['LSST_u', 'LSST_g', 'LSST_r', 'LSST_i', 'LSST_y', 'LSST_z']

    if magtype == 'abs':
        retain += ['i1500']
    
    frame  = frame[retain]
    frame  = frame[:nrows]
    
    if printit:
        for x in attrs:
            print(x)
        
        for x in bands:
            print(x)

    sel    = [x in retain for x in bands]
    bands  = np.array(bands)[sel]

    waves  = np.array(links['mag_wavelengths'])[sel]
    waves  = dict(zip(bands, waves))
            
    return  waves, frame

def get_pyloser_fluxes(boxsize, redshift, printit=False, nrows=-1):
    from  depths  import  get_depths


    depths       = get_depths()
    
    wave, frame  = get_pyloser(boxsize, redshift, printit=printit)
    
    retain       = ['LSST_u', 'LSST_g', 'LSST_r', 'LSST_i', 'LSST_z', 'LSST_y']
    frame        = frame[retain]
    frame        = frame[:nrows]

    wave         = dict(zip(retain, [wave[x] for x in retain]))
    
    for band in retain:
      print('Solving for {}.'.format(band))

      b                               = band.split('_')[-1]
      
      frame['FLAMBDA_' + b]           = np.zeros_like(frame[band])          ##  Fv [ergs/s/cm2/Hz].
      frame['FLAMBDAERR_' + b]        = np.zeros_like(frame[band])
      
      for i, y in enumerate(frame[band]):  
        frame['FLAMBDA_' + b][i]      = 10. ** (-(y + 48.60) / 2.5)
        frame['FLAMBDAERR_' + b][i]   = ferr(y, depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)

        # Fv to Fl conversion.
        frame['FLAMBDA_' + b][i]     *= 2.9979e8 / wave[band] / wave[band]  ##  Fl [ergs/s/cm2/A].
        frame['FLAMBDAERR_' + b][i]  *= 2.9979e8 / wave[band] / wave[band]

        frame['FLAMBDA_' + b][i]     *= 1.0e18                              ##  1.e-18  Fl [ergs/s/cm2/A].
        frame['FLAMBDAERR_' + b][i]  *= 1.0e18 
        
      del frame[band]
        
    return wave, frame

      
if __name__ == '__main__':
    print('\n\nWelcome to Simba get_data.')

    boxsize     =  100.
    '''
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
    '''

    # links       = get_caesar(boxsize, 2.024621)

    waves, frame  = get_pyloser(boxsize, 2.024621, printit=True, magtype='abs')
    # wave, links = get_pyloser_fluxes(boxsize, 2.024621, printit=True, nrows=10)   

    # result      = get_phys(boxsize, 2.024621, printit=False)
        
    print('\n\nDone.\n\n')
