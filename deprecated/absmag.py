import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy             as      np
import  pylab             as      pl
import  matplotlib.pyplot as      plt

from    scipy.spatial     import  KDTree
from    itertools         import  product
from    get_data          import  get_data
from    utils             import  latexify


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.  ##  Mpc/h.
vol          =  boxsize ** 3.

##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244  
f2, p2       =  get_data(boxsize, 2.024621)
f3, p3       =  get_data(boxsize, 3.00307)
f4, p4       =  get_data(boxsize, 3.963392) 
f5, p5       =  get_data(boxsize, 5.024400)

# Color 28 3053 4086  3571  LSST u band 305.30 - 408.60
# Color 29 3863 5670  4768  LSST g band 386.30 - 567.00
# Color 30 5369 7060  6216  LSST r band 536.90 - 706.00
# Color 31 6759 8330  7546  LSST i band 675.90 - 833.00
# Color 32 9083 10996 10042 LSST y band 908.30 - 1099.60
# Color 33 8029 9386  8709  LSST z band 802.90 - 938.60

lsstu2       =  p2['absmag_28'][:]
lsstg2       =  p2['absmag_29'][:]
lsstr2       =  p2['absmag_30'][:]

lsstu3       =  p3['absmag_28'][:]
lsstg3       =  p3['absmag_29'][:]
lsstr3       =  p3['absmag_30'][:]

lsstg4       =  p4['absmag_29'][:]
lsstr4       =  p4['absmag_30'][:]
lssti4       =  p4['absmag_31'][:]

lsstr5       =  p5['absmag_30'][:]
lssti5       =  p5['absmag_31'][:]
lsstz5       =  p5['absmag_33'][:]

umg2         =  lsstu2 - lsstg2
gmr2         =  lsstg2 - lsstr2

umg3         =  lsstu3 - lsstg3
gmr3         =  lsstg3 - lsstr3

gmr4         =  lsstg4 - lsstr4
rmi4         =  lsstr4 - lssti4

imz5         =  lssti5 - lsstz5
rmi5         =  lsstr5 - lssti5

##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

##  plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)
  
axes[0].plot(gmr2, umg2, '.', markersize=.25, c='y', label=r'$z={:.1f}$'.format(2.0246), alpha=0.6)  
axes[1].plot(gmr3, umg3, '.', markersize=.25, c='b', label=r'$z={:.1f}$'.format(3.0031), alpha=0.6)
axes[2].plot(rmi4, gmr4, '.', markersize=.25, c='g', label=r'$z={:.1f}$'.format(3.9633), alpha=0.6)
axes[3].plot(imz5, rmi5, '.', markersize=.25, c='r', label=r'$z={:.1f}$'.format(5.0244), alpha=0.6)

axes[0].set_xlabel(r'$g-r$')
axes[0].set_ylabel(r'$u-g$')

axes[1].set_xlabel(r'$g-r$')
axes[1].set_ylabel(r'$u-g$')

axes[2].set_xlabel(r'$r-i$')
axes[2].set_ylabel(r'$g-r$')

axes[3].set_xlabel(r'$i-z$')
axes[3].set_ylabel(r'$r-i$')

for ax in axes:
  ax.set_axis_on()

  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['left'].set_color('black')
  ax.spines['right'].set_color('black')

  ax.set_xlim(-.5,   .5)
  ax.set_ylim(-.5, 1.75)

  ax.legend(frameon=False, loc=1)
  
plt.tight_layout()

pl.savefig('plots/colorcolor_abs.pdf')
