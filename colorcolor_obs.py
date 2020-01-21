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
from    get_data          import  get_pyloser


boxsize      =  100.

depths       =  get_depths()

two          =  get_pyloser(boxsize, 2.024621)
three        =  get_pyloser(boxsize, 3.003070)
four         =  get_pyloser(boxsize, 3.963392)
five         =  get_pyloser(boxsize, 5.024400)

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]                                                                                                                                                                        
##  Available snapshots: ['062',   '078',    '051',    '042']                                                                                                                                                                          

nruns        =  5 * (len(two['LSST_u'])  + len(three['LSST_u']) + len(four['LSST_u']) + len(five['LSST_u']))
count        =  0

for x in [two, three, four, five]:
  for band in ['LSST_u', 'LSST_g', 'LSST_r', 'LSST_i', 'LSST_z']:
    for i, y in enumerate(x[band]):    
      Flux       =  10. ** (-(y + 48.60) / 2.5)  ##  Nanomaggies. 
      SigF       =  ferr(y, depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)

      Noise      =  normal_rand(loc=0.0, scale=SigF).rvs(size=1)[0]

      ##  See also eqn. (10.2) of Chromey, Introduction to Observational Astronomy.                                                                                                                                               
      x[band][i] = luptitude(Flux + Noise, SigF)

      print(100. * count / nruns)

      count += 1

umg2         =  two['LSST_u'].values   - two['LSST_g'].values
gmr2         =  two['LSST_g'].values   - two['LSST_r'].values

umg3         =  three['LSST_u'].values - three['LSST_g'].values
gmr3         =  three['LSST_g'].values - three['LSST_r'].values

gmr4         =  four['LSST_g'].values  - four['LSST_r'].values
rmi4         =  four['LSST_r'].values  - four['LSST_i'].values

imz5         =  five['LSST_i'].values  - five['LSST_z'].values
rmi5         =  five['LSST_r'].values  - five['LSST_i'].values
      
##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

fast_scatter(axes[0], gmr2, umg2, np.ones_like(gmr2), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[1], gmr3, umg3, np.ones_like(gmr3), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[2], rmi4, gmr4, np.ones_like(rmi4), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)
fast_scatter(axes[3], imz5, rmi5, np.ones_like(imz5), 0.9, 1.1, 10, markersize=0.1, cmap='autumn_r', printit=False, alpha=1.0)

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
  ax.set_ylim(-.5, 2.75)

  ax.legend(frameon=False, loc=1)

plt.tight_layout()

pl.savefig('plots/colorcolor_obs.png')
