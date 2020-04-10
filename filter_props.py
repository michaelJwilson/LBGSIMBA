import matplotlib; matplotlib.use('PDF')

import numpy as np
import pylab as pl
import pickle

from   collections import OrderedDict


## 
cols    = open('/home/mjwilson/LBGSIMBA/pylosers/fsps/data/allfilters.dat', 'r')
lines   = [x.replace('\n','') for x in list(cols.readlines())]

names   = [x for x in lines if x[0] == '#']
names   = [x.replace('#','').lstrip() for x in names]

# Remove dummy name at end.
print('Removing {}'.format(names[-1]))

names   = names[:-1]

print(names)

filters = OrderedDict()

for name in names:
    name = name.replace('#','').lstrip()    
    filters[name] = []
        
for line in lines:
    if line[0] == '#':
        name   = line.replace('#','').lstrip()

        print('Gathering {}'.format(name))
        
        if name not in names:
            print('Ending on {}.'.format(name))

            break
        
        toappend = filters[name]

    else:
        filters[name].append(filter(('').__ne__, line.split(' ')))

for key in filters.keys():
    filters[key] = np.array(filters[key], dtype=np.float64)

with open('bigdat/allfilters.pkl', 'wb') as handle:
    pickle.dump(filters, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
'''
for key in filters.keys():    
  _            = key.split(' ')

  if _[0] == 'LSST':
    wave       = filters[key][:,0]
    trans      = filters[key][:,1] / filters[key][:,1].max()

    dwave      = np.diff(wave)
    
    assert  np.all(np.diff(dwave) < 0.1)
    
    dwave      = dwave[0]

    bandwidth  = dwave * np.sum(trans)
    meanwave   = np.sum(wave * trans) / np.sum(trans)
    
    pl.plot(wave, trans, label=r'{} ${}$:  '.format(_[0], _[1]) + r' ${:.1f}$ +- ${:.1f}\AA$'.format(meanwave, bandwidth / 2.))

  if key == 'Idealized 1500A bandpass: rounded tophat centered at 1500A with 15% bandwidth, FWHM = 225A (MED)':
    wave       = filters[key][:,0]
    trans      = filters[key][:,1] / filters[key][:,1].max()

    dwave      = np.diff(wave)

    assert  np.all(np.diff(dwave) < 0.1)

    dwave      = dwave[0]

    bandwidth  = dwave * np.sum(trans)
    meanwave   = np.sum(wave * trans) / np.sum(trans)
    
    pl.plot(wave, trans, label=r'{}$\AA$: '.format(1500) + r'  ${:.1f}$ +- ${:.1f}\AA$'.format(meanwave, bandwidth / 2.))
    
    
pl.xlim(0000.0, 1.2e4)

pl.legend(frameon=True, loc=3, ncol=2)

pl.savefig('plots/filters.pdf')
'''
