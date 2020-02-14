import matplotlib;  matplotlib.use('PDF')

import itertools

import numpy  as     np
import pylab  as     pl

from   fithod import sat_model

 
ms     = np.logspace(9., 15., 50)

##  mcut, mone, alpha.
mcuts  = list(np.logspace(10.,  12.,  10.))
mones  = list(np.logspace(10.,  12.,  10.))
alphas = list(  np.arange(0.97, 1.24, 0.05))

runs   = zip(mcuts, mones, alphas)
runs   = list(itertools.product(*[mcuts, mones, alphas]))

print(runs)

for mcut, mone, alpha in runs:
    params = [mcut, mone, alpha]
    
    Ns     = sat_model(ms, params)

    label  = '{}  {}  {}'.format(params[0] / 1.e10, params[1] / 1.e10, params[2])
    label  = ''

    pl.loglog(ms, Ns, label=label)

#pl.legend(loc=2)

pl.savefig('satmodel_test.png')
