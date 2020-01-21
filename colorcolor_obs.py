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
from    sphotometry       import  read_mags
from    fast_scatter      import  fast_scatter


boxsize      =  100.

depths       =  get_depths()

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]                                                                                                                                                                        
##  Available snapshots: ['062',   '078',    '051',    '042']                                                                                                                                                                          
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, two,   two_nd    =  read_mags('078', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, three, three_nd  =  read_mags('062', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, four,  four_nd   =  read_mags('051', infile=None, magcols=None, SUFF='app')
redshift, boxsize, nbands, ngal, sfr, LyC, mformed, mstar, L_FIR, meanage, Zstar, A_V, five,  five_nd   =  read_mags('042', infile=None, magcols=None, SUFF='app')

nruns        =  len(two['u']) * len(two.dtype.names) + len(three['u']) * len(three.dtype.names) + len(four['u']) * len(four.dtype.names) + len(five['u']) * len(five.dtype.names)
count        =  0

for x in [two, three, four, five]:
  for band in x.dtype.names:
    for i, y in enumerate(x[band]):    
      Flux       =  10. ** (-(y + 48.60) / 2.5)  ##  Nanomaggies. 
      SigF       =  ferr(y, depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)

      Noise      =  normal_rand(loc=0.0, scale=SigF).rvs(size=1)[0]

      ##  See also eqn. (10.2) of Chromey, Introduction to Observational Astronomy.                                                                                                                                               
      x[band][i] = luptitude(Flux + Noise, SigF)

      print(100. * count / nruns)

      count += 1

exit(1)
      
umg2         =  p2['u'] - p2['g']
gmr2         =  p2['g'] - p2['r']

umg3         =  p3['u'] - p3['g']
gmr3         =  p3['g'] - p3['r']

gmr4         =  p4['g'] - p4['r']
rmi4         =  p4['r'] - p4['i']

imz5         =  p5['i'] - p5['z']
rmi5         =  p5['r'] - p5['i']

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
