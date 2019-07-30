import  matplotlib;  matplotlib.use('PDF')

import  numpy  as  np
import  pylab  as  pl
import  matplotlib.pyplot as plt


from    scipy.optimize  import  minimize 
from    scipy.special   import  erf


def cen_model(mhalo, params):
    mmin    = params[0]
    sigm    = params[1]

    result  = 1. + erf(np.log10(mhalo / mmin) / sigm)

    return  result / 2.
    
def sat_model(mhalo, params):
    mcut   = params[0]
    mone   = params[1]
    alpha  = params[2] 

    result = (mhalo - mcut) / mone
    result = result ##  ** alpha
    
    return  result

def chi2(params, args):
    model     = args['model']
    mhalo, Nx = args['mhalo'], args['Nx']

    err       = np.sqrt(Nx)
    result    = Nx - model(mhalo, params)

    return  np.sum((result / err) ** 2.)

    
if __name__ == '__main__':
    Mh, Nc, Ns = np.loadtxt('dat/hod.txt', unpack=True, dtype=[('Mh', np.float32), ('Nc', np.float32), ('Ns', np.float32)])

    ##  Centrals.
    args       = {'mhalo': Mh, 'Nx': Nc, 'model': cen_model}
    cenparams  = np.array([1.e12, 1.])
    
    result     = minimize(chi2, cenparams, args=args, options={'disp': True}, method='Nelder-Mead')

    print(result.x)
    print(result.success)
    print(result.message)
    
    pl.errorbar(Mh, Nc, yerr=np.sqrt(Nc), c='k', marker='^', linestyle='')
    pl.loglog(Mh, cen_model(Mh, result.x), c='k')

    np.savetxt('dat/hod-nc-params.txt', result.x, fmt='%.6le')
    
    ##  Satellites.
    args       = {'mhalo': Mh, 'Nx': Ns, 'model': sat_model}
    satparams  = np.array([1.e12, 5.e11, 1.2])

    result     = minimize(chi2, satparams, args=args, options={'disp': True}, method='Nelder-Mead')

    print(result.x)
    print(result.success)
    print(result.message)

    pl.errorbar(Mh, Ns, yerr=np.sqrt(Ns), c='darkcyan', marker='^', linestyle='')
    pl.loglog(Mh, sat_model(Mh, result.x), c='darkcyan')
    
    plt.tight_layout()
    
    pl.savefig('plots/fitted-hod.pdf')

    np.savetxt('dat/hod-ns-params.txt', result.x, fmt='%.6le')
