import matplotlib; matplotlib.use('PDF')
import numpy             as     np
import pylab             as     pl
import matplotlib.pyplot as     plt
import pickle

from   get_data          import get_pyloser


boxsize            = 100.
zz                 = 3.96

with open('bigdat/allfilters.pkl', 'rb') as handle:
    filters        = pickle.load(handle)

lssti              = filters['LSST i band 675.90 - 833.00']
lsstz              = filters['LSST z band 802.90 - 938.60']

fpath              = '/home/mjwilson/LBGSIMBA/pylosers/m100n1024/s50/pyloser_m100n1024_051.hdf5'
# retain           = ['i2300', 'i2800']                                                                                                                                                                                                   

retain             = []

wave, absf, ids    = get_pyloser(100., 3.963392, printit=False, magtype='abs', nodust=True, allfilters=True, fpath=fpath, retain=retain)
ww,   appf,  ii    = get_pyloser(100., 3.963392, printit=False, magtype='app', nodust=True, allfilters=True, fpath=fpath, retain=retain)

filters            = appf.columns
filters            = [x for x in filters if x[0] == 'M']

wave               = np.array([x[1:] for x in filters], dtype=np.float)
width              = 100.

absf               = absf[filters]
appf               = appf[filters]

fig, axes          = plt.subplots(2,1)

axes[0].axvspan(1500 - 225./2., 1500.+225./2, alpha=0.25, color='k')
axes[1].plot(lssti[:,0], 25. * lssti[:,1] / lssti[:,1].max())
axes[1].plot(lsstz[:,0], 25. * lsstz[:,1] / lsstz[:,1].max())

for i in range(len(absf)):
  sed  = []
  sed2 = []

  for j, f in enumerate(filters):
    sed.append(absf.at[i, f])
    sed2.append(appf.at[i, f])
    
  axes[0].errorbar(wave, np.array(sed), 0.0, xerr=50., lw=0.2, fmt='', elinewidth=0.4, c='k')
  axes[1].errorbar(wave * (1. + zz), np.array(sed2), 0.0, xerr=50., lw=0.2, fmt='', elinewidth=0.4, c='k')
  
  if i > 100:
    break

axes[0].set_xlim(500., 4500.)
axes[0].set_ylim(-28., -18.)

axes[1].set_xlim(4500., 9500.)
axes[1].set_ylim(21.,     30.)

axes[1].set_xlabel('Wavelength [$\AA$]')

pl.savefig('plots/narrows.pdf')
  
print('\n\nDone.\n\n')
