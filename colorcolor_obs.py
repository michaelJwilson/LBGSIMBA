import  matplotlib;  matplotlib.use('PDF')

import  glob
import  h5py
import  numpy                                 as      np
import  pylab                                 as      pl
import  matplotlib.pyplot                     as      plt

from    scipy.spatial                         import  KDTree
from    itertools                             import  product
from    get_data                              import  get_data
from    sutils                                import  latexify
from    luptitudes                            import  luptitude, lup_lim
from    hildebrandt                           import  ferr
from    depths                                import  get_depths, get_depths27
from    scipy.stats                           import  norm       as normal_rand
from    sphotometry                           import  read_mags
from    insample                              import  insample
from    fast_scatter                          import  fast_scatter
from    get_data                              import  get_pyloser, get_phys
from    color_box                             import  color_box
from    beta_uv                               import  det_bands
from    hildebrandt                           import  onesig_mag
from    mpl_toolkits.axes_grid1.inset_locator import  inset_axes


boxsize            =  100.

depths             =  get_depths()

nrows              =  -1
prop               =  'hmass'  ## 'smass'     

wave, two,   ids2  =  get_pyloser(boxsize, 2.024621, nrows=nrows)

wave, three, ids3  =  get_pyloser(boxsize, 3.003070, nrows=nrows)
wave, four,  ids4  =  get_pyloser(boxsize, 3.963392, nrows=nrows)
wave, five,  ids5  =  get_pyloser(boxsize, 5.024400, nrows=nrows)

phys2              =  get_phys(boxsize, 2.024621, nrows=nrows)['hmass']
phys3              =  get_phys(boxsize, 3.003070, nrows=nrows)['hmass']
phys4              =  get_phys(boxsize, 3.963392, nrows=nrows)['hmass']
phys5              =  get_phys(boxsize, 5.024400, nrows=nrows)['hmass']

##  Available redshifts: [3.00307, 2.024621, 3.963392, 5.0244]                                                                                                                                                                        
##  Available snapshots: ['062',   '078',    '051',    '042']
nruns              =  5 * (len(two['LSST_u'])  + len(three['LSST_u']) + len(four['LSST_u']) + len(five['LSST_u']))
count              =  0

bands              =  ['LSST_u', 'LSST_g', 'LSST_r', 'LSST_i', 'LSST_z', 'LSST_y']

for redshift, x, ids in zip([2.024621, 3.003070, 3.963392, 5.0244], [two, three, four, five], [ids2, ids3, ids4, ids5]):
  _, dbands        = det_bands(redshift, wave, bands, lim=1500.)

  x['TARGETIDS']   = ids

  x['DBANDS']      = ''.join(x.split('_')[-1] for x in dbands)
  x['ISIN']        = np.array([''] * len(x['LSST_u']))
  x['SNR2']        = np.zeros_like(x['LSST_u'])
  
  for band in bands:
    x['LUP_' + band]  = np.zeros_like(x[band])
    
    onesig            = onesig_mag(depths[band])                                                 ##  AB mags.

    print('Solving for {} depth in {} ({}).'.format(onesig, band, ''.join(x.split('_')[-1] for x in dbands)))

    ##  Magnitudes.
    x['FIVSIG_' + band] = depths[band]
    x['ONESIG_' + band] = onesig
    
    ##  Fvs.
    onesig              = 10. ** -((onesig + 48.60) / 2.5)
    x['LUPLIM_' + band] = lup_lim(onesig)
  
    for i, y in enumerate(x[band]):    
      Flux       =  10. ** (-(y + 48.60) / 2.5)                                                 ##  Fv.

      ##  estar  = 0.2 sets 5 sigma depth provided. 
      SigF       =  ferr(y, depths[band], estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None)

      Noise      =  normal_rand(loc=0.0, scale=SigF).rvs(size=1)[0]
            
      ##  See also eqn. (10.2) of Chromey, Introduction to Observational Astronomy.
      ##  x.at[i, 'LUP_' + band] = luptitude(Flux,         onesig)

      if Flux > SigF:
        x.at[i, 'LUP_' + band]   = luptitude(Flux + Noise, onesig)

      else:
        x.at[i, 'LUP_' + band]   = lup_lim(onesig)

      
      x.at[i, 'SNR2']           += (Flux / SigF)**2.

      if band in dbands:
        #  https://www.dataquest.io/blog/settingwithcopywarning/
        x.at[i, 'ISIN']         +=  band.split('_')[-1] if (x.at[i, band] <= depths[band]).astype(np.int) else ''
      
      # print(100. * count / nruns)
      
      count += 1
    

physs  = {}
frames = {}

##  Sufficiently detected. 
for name, frame, phys, selection in zip(['two', 'three', 'four', 'five'], [two, three, four, five], [phys2, phys3, phys4, phys5], ['BM', 'u', 'g', 'r']):
  isin              = [len(x) >= 1 for x in frame['ISIN']]
  isin              = isin & insample(selection, frame['LUP_LSST_u'].values, frame['LUP_LSST_g'].values, frame['LUP_LSST_r'].values, frame['LUP_LSST_i'].values, frame['LUP_LSST_z'].values, frame['LUP_LSST_y'].values)

  frame['INSAMPLE'] = isin
  framae            = frame.drop(['ISIN'], axis=1)
  
  ##  frame.to_pickle('bigdat/insample_{}.pkl'.format(name))
  frame.to_hdf('bigdat/insample_{}.h5'.format(name), key='df', mode='w')  
  
  print('\n\n')
  print(frame.loc[isin, :])

  physs[name]       = phys
  frames[name]      = frame 
  
##
two          =  frames['two']
three        =  frames['three']
four         =  frames['four']
five         =  frames['five']

phys2        =  physs['two']
phys3        =	physs['three']
phys4        =	physs['four']
phys5        =	physs['five']


##
umg2         =    two['LUP_LSST_u'].values - two['LUP_LSST_g'].values
gmr2         =    two['LUP_LSST_g'].values - two['LUP_LSST_r'].values

umg3         =  three['LUP_LSST_u'].values - three['LUP_LSST_g'].values
gmr3         =  three['LUP_LSST_g'].values - three['LUP_LSST_r'].values

gmr4         =   four['LUP_LSST_g'].values - four['LUP_LSST_r'].values
rmi4         =   four['LUP_LSST_r'].values - four['LUP_LSST_i'].values

imz5         =   five['LUP_LSST_i'].values - five['LUP_LSST_z'].values
rmi5         =   five['LUP_LSST_r'].values - five['LUP_LSST_i'].values


##
latexify(columns=2, equal=False, fontsize=8, ggplot=True, usetex=True, ratio=0.35)

##
fig, axes    = plt.subplots(nrows=1, ncols=4, sharey=False)

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.75, hspace=None)

cmap         = 'viridis'

isin         = two['INSAMPLE'].values
scatter      = axes[0].scatter(gmr2[ isin], umg2[ isin], marker='.', c=np.log10(phys2[ isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14.)
scatter      = axes[0].scatter(gmr2[~isin], umg2[~isin], marker='x', c=np.log10(phys2[~isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14., alpha=0.5)

isin         = three['INSAMPLE'].values
scatter      = axes[1].scatter(gmr3[ isin], umg3[ isin], marker='.', c=np.log10(phys3[ isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14.)
scatter      = axes[1].scatter(gmr3[~isin], umg3[~isin], marker='x', c=np.log10(phys3[~isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14., alpha=0.7)

isin         = four['INSAMPLE'].values
scatter      = axes[2].scatter(rmi4[ isin], gmr4[ isin], marker='.', c=np.log10(phys4[ isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14.)
scatter      = axes[2].scatter(rmi4[~isin], gmr4[~isin], marker='x', c=np.log10(phys4[~isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14., alpha=0.7)

isin         = five['INSAMPLE'].values
scatter      = axes[3].scatter(imz5[ isin], rmi5[ isin], marker='.', c=np.log10(phys5[ isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14.)
scatter      = axes[3].scatter(imz5[~isin], rmi5[~isin], marker='x', c=np.log10(phys5[~isin]), lw=0, s=5, cmap=cmap, vmin=9.5, vmax=14., alpha=0.7)


##  
axins3 = inset_axes(axes[3], width="70%", height="5%", loc='lower center')

cbar   = fig.colorbar(scatter, cax=axins3, orientation="horizontal", ticks=[9.5, 11.75, 14.0])

cbar.ax.set_yticklabels(['9.5', '11.75', '14.0'])

axins3.xaxis.set_ticks_position("top")

color_box(axes[1], 'u')
color_box(axes[2], 'g')
color_box(axes[3], 'r')

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

  ax.set_xlim(-1.5, 2.00)
  ax.set_ylim(-2.0, 3.55)

  ax.legend(frameon=False, loc=1)

plt.tight_layout()

pl.savefig('plots/colorcolor_obs.png')
