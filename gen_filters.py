import matplotlib; matplotlib.use('PDF')

import numpy as np
import pylab as pl

from   collections import OrderedDict


## --- FILTER LIST ---
cols    = open('/home/mjwilson/LBGSIMBA/pylosers/fsps/data/allfilters/FILTER_LIST')
lines   = [x.replace('\n','') for x in list(cols.readlines())]

centers  = np.arange(500, 8000, 100)
ifilt    = 1 + np.arange(len(lines), len(lines) + len(centers), 1)
names    = ['M{}'.format(x) for x in centers]
comments = ['Idealised tophat centered on {}A with 100A width.'.format(x) for x in centers]

nnames   = []  

for i, _ in enumerate(ifilt):
  nnames += ['{} {} \t\t{}'.format(ifilt[i], names[i], comments[i])]
  lines  += ['{} {} \t\t{}'.format(ifilt[i], names[i], comments[i])]

print(','.join(names))
  
with open('/home/mjwilson/LBGSIMBA/pylosers/fsps/data/FILTER_LIST_v2', 'w') as f:
    for item in lines:
        f.write("{}\n".format(item))
  
## 
cols    = open('/home/mjwilson/LBGSIMBA/pylosers/fsps/data/allfilters/allfilters.dat', 'r')
lines   = [x.replace('\n','') for x in list(cols.readlines())]

keep    = lines[:-4]

for i, name in enumerate(names):
  keep  += ['# {}'.format(name)]

  cen    = centers[i]
  wave   = np.arange(cen - 65., cen + 65., 10.)

  trans  = np.ones_like(wave)

  trans[wave  > cen + 50.] = 0.0
  trans[wave  < cen - 50.] = 0.0

  for w, t in zip(wave, trans):
    keep += [' {:.1f}   {:.6f}'.format(w, t)]
    
keep    += lines[-4:]

#for l in keep:
#  print(l)

with open('/home/mjwilson/LBGSIMBA/pylosers/fsps/data/allfilters_v2.dat', 'w') as f:
    for item in keep:
        f.write("{}\n".format(item))

        # print("{}".format(item))
