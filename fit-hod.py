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
    mcut   = params['mcut']
    mone   = params['mone']
    alpha  = params['alpha'] 

    result = (mhalo - mcut) / mone

    return  result ** alpha

def chi2(params, args):
    model     = args['model']
    mhalo, Nx = args['mhalo'], args['Nx']

    result    = Nx - model(mhalo, params)

    return  np.sum(result ** 2.)

    
if __name__ == '__main__':
    Mh, Nc, Ns = np.loadtxt('dat/hod.txt', unpack=True, dtype=[('Mh', np.float32), ('Nc', np.float32), ('Ns', np.float32)])

    args       = {'mhalo': Mh, 'Nx': Nc, 'model': cen_model}

    ##  np.array([1.e10, 1.0])   
    cenparams  = np.array([1.e12, 1.])

    for x, y in zip([1.e10, 1.e11, 1.e12], [1.0, 2.5, 5.0]):
        print(chi2(np.array([x,y]), args))

    
    result     = minimize(chi2, cenparams, args=args, options={'disp': True}, method='Nelder-Mead')

    print(result.x)
    print(result.success)
    print(result.message)
    
    pl.plot(Mh, Nc, 'k^')
    pl.plot(Mh, cen_model(Mh, result.x))

    plt.tight_layout()
    
    pl.savefig('plots/fitted-hod.pdf')
