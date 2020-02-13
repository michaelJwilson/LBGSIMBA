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
    result = result**alpha
    
    return  result

def chi2(params, args):
    model     = args['model']
    mhalo, Nx = args['mhalo'], args['Nx']

    uerr      = args['uerr']
    
    err       = np.sqrt(Nx)
    result    = Nx - model(mhalo, params)

    if uerr:
        err   = np.ones_like(result)   
    
    return  np.sum((result / err) ** 2.)

    
def fit_hod(boxsize=100., getredshift=3.00307):
    fname      = 'dat/hod_{}.txt'.format(str(getredshift).replace('.', 'p'))

    Mh, Nc, Ns = np.loadtxt(fname, unpack=True, dtype=[('Mh', np.float32), ('Nc', np.float32), ('Ns', np.float32)])

    print(Mh)
    print(Nc)
    print(Ns)
    
    ##  Centrals.
    args       = {'mhalo': Mh, 'Nx': Nc, 'model': cen_model, 'uerr': 1}
    cenparams  = np.array([5.e11, 0.5])
    
    result     = minimize(chi2, cenparams, args=args, options={'disp': True, 'maxiter': 10000}, method='Nelder-Mead')

    print(result.x)
    print(result.success)
    print(result.message)
    
    pl.errorbar(Mh, Nc, yerr=np.sqrt(Nc), c='k', marker='^', linestyle='')
    pl.loglog(Mh, cen_model(Mh, result.x), c='k')

    np.savetxt('dat/hod-nc-params_{}.txt'.format(str(getredshift).replace('.', 'p')), result.x, fmt='%.6le')
    
    ##  Satellites.
    args       = {'mhalo': Mh, 'Nx': Ns, 'model': sat_model, 'uerr': 0}
    satparams  = np.array([5.2e11, 1.e12, 0.8])

    result     = minimize(chi2, satparams, args=args, options={'disp': True, 'maxiter': 10000}, method='Nelder-Mead')

    print(result.x)
    print(result.success)
    print(result.message)
    
    pl.errorbar(Mh, Ns, yerr=np.sqrt(Ns), c='darkcyan', marker='^', linestyle='')

    # pl.loglog(Mh, sat_model(Mh, result.x),  c='darkcyan')
    pl.loglog(Mh, sat_model(Mh, satparams), c='darkcyan')
    
    plt.tight_layout()
    
    pl.savefig('plots/fitted-hod.pdf')

    np.savetxt('dat/hod-ns-params_{}.txt'.format(str(getredshift).replace('.', 'p')), result.x, fmt='%.6le')


if __name__ == '__main__':
    print('\n\nWelcome to fit hod..\n\n')

    redshifts = [2.024621, 3.00307, 3.963392, 5.0244]

    for redshift in redshifts:
        fit_hod(100., getredshift=redshift)

        break
        
    print('\n\nDone.\n\n')
 
