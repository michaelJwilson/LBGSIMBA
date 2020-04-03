import numpy as np
import pylab as pl


def color_box(ax, dband = 'g'):    
    ##  Hildebrandt (2018), https://arxiv.org/pdf/0903.3951.pdf;  Otherwise, GoldRush (Ono, 2018).  
    hsc        = {'u': {'bcol': 'u-g', 'rcol': 'g-r', 'minbcol': 1.5, 'maxrcol': 1.2, 'gradient': 1.5, 'intercept': 0.75, 'hiz': 4.5},\
                  'g': {'bcol': 'g-r', 'rcol': 'r-i', 'minbcol': 1.0, 'maxrcol': 1.0, 'gradient': 1.0, 'intercept': 1.00, 'hiz': 6.0},\
                  'r': {'bcol': 'r-i', 'rcol': 'i-z', 'minbcol': 1.2, 'maxrcol': 0.7, 'gradient': 1.5, 'intercept': 1.00, 'hiz': 7.5},\
                  'z': {'bcol': 'i-z', 'rcol': 'z-y', 'minbcol': 1.5, 'maxrcol': 0.5, 'gradient': 2.0, 'intercept': 1.10, 'hiz': 7.5}}
    
    bcol       =  hsc[dband]['bcol']
    rcol       =  hsc[dband]['rcol']

    minbcol    =  hsc[dband]['minbcol']        ## Detect the break
    maxrcol    =  hsc[dband]['maxrcol']        ## Flat spectra above break. 
    
    gradient   =  hsc[dband]['gradient']
    intercept  =  hsc[dband]['intercept']

    rcols      =  np.arange(-1.5, 4.0, 0.001)
    bcols      =  gradient * rcols + intercept

    minrcol    = (minbcol - intercept) / gradient
    bcollim    =  gradient * maxrcol + intercept

    ##  Plot selection box perimeter. 
    ax.plot(rcols[rcols < minrcol], minbcol * np.ones_like(rcols[rcols < minrcol]), 'k', lw=0.4)
    ax.plot(maxrcol * np.ones_like(bcols[bcols > bcollim]), bcols[bcols > bcollim], 'k', lw=0.4)
    
    ##  Gradient
    ax.plot(rcols[(rcols > minrcol) & (bcols < bcollim)], bcols[(rcols > minrcol) & (bcols < bcollim)], c='k', linestyle='-', lw=0.4)

    return  0

def color_box_steidel(ax, ttype = 'BX'):
    ##  https://arxiv.org/pdf/astro-ph/0401445.pdf
    steidel    = {'bx': {'bcol': 'Un-G', 'rcol': 'G-R', 'minrcol': -0.2, 'gradient': 1.0, 'intercept': 0.2},\
                  'bm': {'bcol': 'Un-G', 'rcol': 'G-R', 'minrcol': -0.2, 'gradient': 1.0, 'intercept': 0.2}}

    bcol       =  steidel[ttype]['bcol']
    rcol       =  steidel[ttype]['rcol']

    minrcol    =  steidel[ttype]['minrcol']                                                                                                                                                 

    gradient   =  steidel[ttype]['gradient']
    intercept  =  steidel[ttype]['intercept']

    rcols      =  np.arange(-1.5, 4.0, 0.001)

    bcollims   = []
    bcolmaxs   = []
    
    for i, intercept in enumerate([-1.0, 0.2, 1.0]):
      bcols    =  gradient * rcols + intercept

      rcollim  =  rcols < 0.2 * bcols + 0.4

      bcollims.append(bcols[(rcols > minrcol) & rcollim].min())
      bcolmaxs.append(bcols[(rcols > minrcol) & rcollim].max())
      
      ax.plot(rcols[(rcols > minrcol) & rcollim], bcols[(rcols > minrcol) & rcollim], c='k', linestyle='-', lw=0.4)

    bcollims   =  np.array(bcollims) 
    bcolmaxs   =  np.array(bcolmaxs)
    
    ax.plot(minrcol * np.ones_like(np.arange(bcollims[0], bcollims[-1], 0.001)), np.arange(bcollims[0], bcollims[-1], 0.001), lw=0.4, c='k')
    ax.plot(rcols[((5.0 * (rcols - 0.4)) > bcolmaxs[0]) & ((5.0 * (rcols - 0.4)) < bcolmaxs[-1])], (5.0 * (rcols - 0.4))[((5.0 * (rcols - 0.4)) > bcolmaxs[0]) & ((5.0 * (rcols - 0.4)) < bcolmaxs[-1])], lw=0.4, c='k')

    print((5.0 * (rcols - 0.4))[(bcols > bcolmaxs[0]) & (bcols < bcolmaxs[-1])])
    
    return  0


if __name__ == '__main__':
    colourtrack('g')
    
    print('\n\nDone.\n\n')
