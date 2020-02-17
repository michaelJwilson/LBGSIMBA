import matplotlib;  matplotlib.use('PDF')

from   AbacusCosmos      import Halos

import matplotlib.pyplot as plt
import numpy             as np
import Corrfunc


fpath         = '/disk01/mjwilson/abacus/AbacusCosmos_720box_planck_products/AbacusCosmos_720box_planck_00-0_products/AbacusCosmos_720box_planck_00-0_FoF_halos/z1.500'

# The path to the catalog will be different on your system
cat           = Halos.make_catalog_from_dir(dirname=fpath, load_subsamples=False, load_pids=False, halo_type='FoF')

halos         = cat.halos

##  [-360., 360.] to [0., 720.]
halos['pos'] += 360.0

for field in sorted(halos.dtype.fields):
    print(field, ':', halos[field][:10])

##  
bins          = np.logspace(0., 2.23, 75)
rs            = (bins[:-1] + bins[1:]) / 2.

exit(0)

results       = Corrfunc.theory.xi(X=halos['pos'][:,0], Y=halos['pos'][:,1], Z=halos['pos'][:,2], boxsize=720., nthreads=4, binfile=bins)

plt.loglog(rs, rs * rs * results['xi'])
plt.legend()
plt.xlabel(r'$s$ [Mpc/h]')
plt.ylabel(r'$s^2 \cdot \xi(s)$')
plt.savefig('plots/abacus.pdf')
