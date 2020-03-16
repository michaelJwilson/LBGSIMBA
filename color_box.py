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


if __name__ == '__main__':
    colourtrack('g')
    
    print('\n\nDone.\n\n')
