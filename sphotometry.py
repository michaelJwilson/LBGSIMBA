# Reads in the photometry files output by Loser (https://bitbucket.org/romeeld/closer/src/default/)
# These are ASCII files, with suffix .app for apparent mags, .abs for absolute mags.

import pylab    as plt
import numpy    as np
import sys
import os
import function as fu
import caesar

from   get_data import snaps


###########################################################################
def read_mags(SNAP, infile=None, magcols=None, SUFF='app'):
    if infile == None:
      infile   = '/home/rad/data/m100n1024/s50/Groups/loser_m100n1024_{}.{}'.format(SNAP, SUFF)

    if magcols == None:
        # ugrizyYJH
        magcols = [28, 29, 30, 31, 33, 32, 25, 26, 27]

    print('Reading file: {}'.format(infile))
        
    redshift, boxsize                                                     = np.genfromtxt(infile,usecols=(2,10),unpack=True,comments=None,max_rows=1)
    galid, _, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, nband = np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10),unpack=True)

    nbands  = int(nband[0])
    ngal    = len(sfr)

    Lmag    = []
    Lmag_nd = []

    for i in range(len(magcols)):
        Lmag.append(np.loadtxt(infile,usecols=(11+int(magcols[i])),unpack=True))
        Lmag_nd.append(np.loadtxt(infile,usecols=(12+int(magcols[i])+nbands),unpack=True))

    # magnitudes of galaxies in each desired band
    mag    = np.asarray(Lmag, np.dtype([('u',float), ('g',float), ('r',float), ('i',float), ('z',float), ('y',float), ('Y',float), ('J',float), ('H',float)]))

    # no-dust magnitudes
    mag_nd = np.asarray(Lmag_nd)

    return redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, mag, mag_nd


if __name__ == '__main__':
    SNAP = '062'
    
    print('\n\nWelcome to SIMBA LBG.\n\n')

    print('Available redshifts: {}'.format(snaps.keys()))
    print('Available snapshots: {}'.format(snaps.values()))

    redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, mag, mag_nd = read_mags(SNAP)
          
    print('Snapshot %s: at redshift %.6lf'   % (SNAP, redshift))
    print('With boxsize %.6lf Mpc/h'         % boxsize)
    print('Number of available filters:  %s' % nbands)
    print('Number of available galaxies: %d' % ngal)

    ## Transpose magnitudes.  nd ->  no intrinsic dust extinction. 
    mag, mag_nd = mag.T, mag_nd.T

    print('\n\nExample magnitudes:\n\n')
    
    ##  Magnitudes in all requested bands of first 10 galaxies; 
    print(mag['u'][:10])
    print(mag['g'][:10])
    print(mag['r'][:10])
    
    print('\n\nDone.\n\n')
