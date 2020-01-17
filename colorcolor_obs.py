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
from    luptitudes        import  luptitude
from    hildebrandt       import  ferr
from    depths            import  get_depths
from    scipy.stats       import  norm       as normal_rand


##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244
boxsize      =  100.
vol          =  boxsize ** 3.

depths       =  get_depths()

##  Closest redshifts:  2.024621, 3.00307, 3.963392, 5.0244  
f2, _p2      =  get_data(boxsize, 2.024621)
f3, _p3      =  get_data(boxsize, 3.00307)
f4, _p4      =  get_data(boxsize, 3.963392) 
f5, _p5      =  get_data(boxsize, 5.024400)

# Color 28 3053 4086  3571  LSST u band 305.30 - 408.60
# Color 29 3863 5670  4768  LSST g band 386.30 - 567.00
# Color 30 5369 7060  6216  LSST r band 536.90 - 706.00
# Color 31 6759 8330  7546  LSST i band 675.90 - 833.00
# Color 32 9083 10996 10042 LSST y band 908.30 - 1099.60
# Color 33 8029 9386  8709  LSST z band 802.90 - 938.60

##  u, g, r.
p2           =  np.hstack((_p2['appmag_28'][:], _p2['appmag_29'][:], _p2['appmag_30'][:]))
p2           =  np.array(p2, dtype = [('u', 'float32'), ('g','float32'), ('r','float32')])

##  u, g, r.
p3           =  np.hstack((_p3['appmag_28'][:], _p3['appmag_29'][:], _p3['appmag_30'][:]))
p3           =  np.array(p3, dtype = [('u', 'float32'), ('g','float32'), ('r','float32')])

##  g, r, i.
p4           =  np.hstack((_p4['appmag_29'][:], _p4['appmag_30'][:], _p4['appmag_31'][:]))
p4           =  np.array(p4, dtype = [('g', 'float32'), ('r', 'float32'), ('i', 'float32')])

##  r, i, z. 
p5           =  np.hstack((_p5['appmag_30'][:], _p5['appmag_31'][:], _p5['appmag_33'][:]))
p5           =  np.array(p5, dtype = [('r', 'float32'), ('i', 'float32'), ('z', 'float32')])

for band in ['u', 'g', 'r']:
  for i, y in enumerate(x['appmag_{:d}'.format(col)][:]):    
      Flux      =  10. ** (-(y + 48.60) / 2.5)
      SigF      =  ferr(y, depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)

      Noise     =  normal_rand(loc=0.0, scale=SigF).rvs(size=1)[0]

      ##  See also eqn. (10.2) of Chromey, Introduction to Observational Astronomy.                                                                                                                                               
      lup       =  luptitude(Flux + Noise, SigF)

      x['appmag_{:d}'.format(col)][i] = lup

exit(1)

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

pl.savefig('plots/colorcolor_obs.pdf')
