import matplotlib; matplotlib.use('PDF')

import pylab as pl
import numpy as np
import matplotlib.pyplot as plt

from   get_data import get_pyloser, get_pyloser_spectra, get_pyloser_AV, get_phys


##  3.003070, 3.963392, 5.0244
# zz               = 5.024400
zz                 = 3.963392
# zz               = 3.003070

waves, frame, ids  = get_pyloser(100., zz, nodust=False)
phys               = get_phys(100., zz)
wave,  spec        = get_pyloser_spectra(100., zz)

assert  len(phys['sfr']) == len(frame)

colors             = plt.rcParams['axes.prop_cycle'].by_key()['color']

lines              = [1215.24, 1264.738, 1309.276, 1533.431, 1660.809, 1666.150, 1908.734]

absorption         = [1260.4221, 1302.1685, 1304.3702, 1334.5323, 1393.76018, 1402.77291,\
                      1526.70698, 1548.204, 1550.781, 1608.45085, 1670.7886]

##
for line in lines:
  pl.axvline(x=(1. + zz)*line, linestyle='--', c='k', lw=0.5)

##
pl.axvline(x=(1. + zz)*1260.4221, linestyle='--', c='r', lw=0.5)
  
for line in absorption:
  pl.axvline(x=(1. + zz)*line, linestyle='--', c='g', lw=0.0)

for i, x in enumerate(spec):
  norm            = 10. ** (0.4 * (frame['LSST_y'][i] - 22.5))
  spec[i]        *= norm

name              = 'sfr'
prop              = phys[name]
cut               = np.percentile(prop, 90)
sample            = prop > cut

title             = 'Stacked for $\log_{10}$' + '({0:}) > {1:.2f} sample: {2:.2f}%'.format(name.upper(), np.log10(cut), 100. * np.count_nonzero(sample) / len(prop))

print(title)

stacked           = np.mean(spec[sample], axis=0)

pl.plot(wave, stacked, alpha=0.5, c='k')
pl.title(title)
pl.xlim(3.5e3, 1.e4)
  
pl.xlabel(r'Wavelength $[\AA]$')
pl.ylabel(r'$F_{\lambda}$ [ergs/$s$/cm$^2$/$\AA$]')

pl.savefig('plots/spectra.pdf')

pl.clf()

##  Extinction plot.
for zz in [2.024621, 3.003070, 3.963392, 5.024400]:
  AV              = get_pyloser_AV(100., zz) 
  RV              = 3.1
  EBV             = AV / RV

  pl.hist(EBV, bins=np.arange(-0.01, 1., 0.01), histtype='step', label=r'z=${:.2f}$'.format(zz), density=True)

pl.xlabel(r'$E(B-V \ | \ R_V=3.1)$')
pl.ylabel(r'$P(E(B-V))$')
pl.yscale('log')
pl.legend(frameon=False)
pl.savefig('plots/ebv.pdf')

print('\n\nDone.\n\n')
