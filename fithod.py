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

    result[mhalo < mcut] = 0.0
    
    return  result

def chi2(params, args):
    model     = args['model']
    mhalo, Nx = args['mhalo'], args['Nx']

    uerr      = args['uerr']
    err       = args['std']
    
    result    = Nx - model(mhalo, params)

    if uerr:
        err   = np.ones_like(result)   
    
    return  np.sum((result / err) ** 2.)

    
def fit_hod(boxsize=100., getredshift=3.00307, set_insample=0):
    print('\n\nLoading measured HOD.\n\n')

    fname                = 'dat/hod_{}_insample_{}.txt'.format(str(getredshift).replace('.', 'p'), set_insample)
    Mh, Nc, sNc, Ns, sNs = np.loadtxt(fname, unpack=True, dtype=[('Mh', np.float32), ('Nc', np.float32), ('sNc', np.float32), ('Ns', np.float32), ('sNs', np.float32)])
    
    isin                 = Ns > 0.0
    
    ##  Maximum of the two variance estimates. 
    sNs                  = np.maximum(sNs, np.sqrt(Ns))
    
    for i, _ in enumerate(Mh):
        print('{:e} \t {:e} \t {:e} \t {:e} \t {:e}'.format(_, Nc[i], sNc[i], Ns[i], sNs[i]))

    ##  Centrals.
    print('\n\nSolving for centrals.\n\n')

    args       = {'mhalo': Mh,\
                  'Nx': Nc,\
                  'std': sNc,\
                  'model': cen_model,\
                  'uerr': 1}

    cenparams  = np.array([5.e11, 0.5])
    
    result     = minimize(chi2, cenparams, args=args, options={'disp': True, 'maxiter': 10000}, method='Nelder-Mead')

    print('\n')
    print(result.x)
    print(result.success)
    print(result.message)

    pl.clf()
    pl.errorbar(Mh, Nc, yerr=sNc, c='k', marker='^', linestyle='')
    pl.loglog(Mh, cen_model(Mh, result.x), c='k')

    np.savetxt('dat/hod-nc-params_{}_insample_{}.txt'.format(str(getredshift).replace('.', 'p'), set_insample), result.x, fmt='%.6le')

    ##  Satellites.
    print('\n\nSolving for satellitess.')

    args       = {'mhalo': Mh[isin], 'Nx': Ns[isin], 'std': sNs[isin], 'model': sat_model, 'uerr': 0}
    satparams  = np.array([5.2e11, 1.e12, 0.96])
    
    result     = minimize(chi2, satparams, args=args, options={'disp': True, 'maxiter': 10000}, method='Nelder-Mead')

    print('\n')
    print(result.x)
    print(result.success)
    print(result.message)
    
    pl.errorbar(Mh, Ns, yerr=sNs, c='darkcyan', marker='^', linestyle='')

    pl.loglog(Mh, sat_model(Mh, result.x),  c='darkcyan')

    pl.ylim(1.e-2, 1.e2)
    
    plt.tight_layout()
    
    pl.savefig('plots/fitted-hod_{}_insample_{}.pdf'.format(str(getredshift).replace('.', 'p'), set_insample))

    np.savetxt('dat/hod-ns-params_{}_insample_{}.txt'.format(str(getredshift).replace('.', 'p'), set_insample), result.x, fmt='%.6le')


if __name__ == '__main__':
    print('\n\nWelcome to fit hod.')

    set_insample =  0
    redshifts    = [2.024621, 3.00307, 3.963392]

    for redshift in redshifts:
        print('\n\nSolving for redshift: {}'.format(redshift))
        
        fit_hod(100., getredshift=redshift, set_insample=set_insample)
        
    print('\n\nDone.\n\n')
 
