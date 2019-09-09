# Reads in the photometry files output by Loser (https://bitbucket.org/romeeld/closer/src/default/)
# These are ASCII files, with suffix .app for apparent mags, .abs for absolute mags.

import pylab as plt
import numpy as np
import sys
import os
import function as fu
import caesar

# choose simulation and snapshot to look at
MODEL = sys.argv[1]  # e.g. m50n512 - size and particle number
WIND  = sys.argv[2]  # e.g. s50 for Simba
SNAP  = sys.argv[3]  # snapshot number
SUFFIX = sys.argv[4]  # "abs" for absolute mags, "app" for apparent mags
magcols = [0]  # default colors
if len(sys.argv) > 5:
    magcols = sys.argv[5:] # the desired color numbers from the header of the *.app or *.abs files.  e.g. 0 for Johnson V, and so on.

###########################################################################

def read_mags(infile,magcols):
    redshift,boxsize = np.genfromtxt(infile,usecols=(2,10),unpack=True,comments=None,max_rows=1)
    galid,junk,sfr,LyC,mformed,mstar,L_FIR,meanage,Zstar,A_V,nband = np.loadtxt(infile,usecols=(0,1,2,3,4,5,6,7,8,9,10),unpack=True)
    nbands = int(nband[0])
    ngal = len(sfr)
    Lmag = []
    Lmag_nd = []
    for i in range(len(magcols)):
        Lmag.append(np.loadtxt(infile,usecols=(11+int(magcols[i])),unpack=True))
        Lmag_nd.append(np.loadtxt(infile,usecols=(12+int(magcols[i])+nbands),unpack=True))
    Lmag = np.asarray(Lmag)        # magnitudes of galaxies in each desired band
    Lmag_nd = np.asarray(Lmag_nd)  # no-dust magnitudes
    return redshift,boxsize,nbands,ngal,sfr,LyC,mformed,mstar,L_FIR,meanage,Zstar,A_V,Lmag,Lmag_nd


if __name__ == '__main__':
    lfile = '/home/rad/data/%s/%s/Groups/loser_%s_%03d.%s' % (MODEL,WIND,MODEL,int(SNAP),SUFFIX)
    
    redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, Lmag, Lmag_nd = read_mags(lfile, magcols)

    print('\n\nWelcome to SIMBA LBG.\n\n')

    print('Reading: %s\n\n' % lfile)
    
    print('Snapshot %s: at redshift %.6lf' % (SNAP, redshift))
    print('With boxsize %.6lf Mpc/h' % boxsize)
    print('Number of filters request:  %s' % nbands)
    print('Number of available galaxies: %d' % ngal)

    print('\n\nM* for first 10 galaxies:')

    for Ms in mstar[:10]:
        print('%.6le' % Ms)

    ## Transpose magnitudes.  nd ->  no intrinsic dust extinction. 
    Lmag, Lmag_nd = Lmag.T, Lmag_nd.T

    print('\n\nExample magnitudes:\n\n')
    
    ##  Magnitudes in all requested bands of first 10 galaxies; 
    for mags in Lmag[:10]:
        print(mags)
    
    print('\n\nDone.\n\n')
