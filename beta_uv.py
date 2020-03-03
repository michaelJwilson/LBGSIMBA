import  matplotlib;  matplotlib.use('PDF')

import  json
import  numpy              as     np
import  pylab              as     pl
import  pandas             as     pd
import  matplotlib.pyplot  as     plt

from    get_data           import get_pyloser, get_pyloser_fluxes
from    scipy.optimize     import minimize
from    sutils             import latexify
from    insample           import read_insample


colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']

def continuum_model(waves, f0, beta):
  return  f0 * waves**beta

def chi2(params, waves, flambdas, flambda_errs):
  f0     = params[0]
  beta   = params[1] 

  result = 0.0

  model  = continuum_model(waves, f0, beta)
  
  for i, x in enumerate(flambdas):
      result += ((x - model[i]) / flambda_errs[i])**2.

  return  result

def minimize_chi(chi2, waves, flambdas, flambda_errs):
    x0             = [1.e-10, 0.0]

    res            = minimize(chi2, x0, args=(waves, flambdas, flambda_errs), method='Nelder-Mead', tol=1e-6)

    chi2, f0, beta = res.fun, res.x[0], res.x[1]
    
    return  chi2, f0, beta

def det_bands(redshift, wave, bands, lim=1500.):
  ##  Eqn. (1) of https://arxiv.org/pdf/1109.0994.pdf                                                                                                                                                                          
  lowave  = lim * (1. + redshift)          ##  Angstroms
  
  bands   = np.array(list(wave.keys()))    ##  Ordered.                                                                                                                                                                          
  waves   = np.array(list(wave.values()))

  #print(bands)
  #print(waves)
  
  indx    = np.argsort(waves)

  waves   = waves[indx]
  bands   = bands[indx]

  #print(waves)
  #print(bands)
  
  retain  = [True if wave[band] >= lowave else False for band in bands]

  #print(retain)
  
  ##  Observed. mags. not affectded by IGM (rest-framce wave > 1500. A).                                                                                                                                                       
  bands   = bands[retain]
  waves   = waves[retain]

  return  waves, bands
  
def solve_chi2s(chi2, redshift, wave, frame):    
  print('\n\nSolving for betas.\n\n')

  ##  Eqn. (1) of https://arxiv.org/pdf/1109.0994.pdf
  lowave  = 1500. * (1. + redshift)      ##  Angstroms

  bands   = np.array(list(wave.keys()))  ##  Ordered. 
  waves   = np.array(list(wave.values()))

  indx    = np.argsort(waves)

  waves   = waves[indx]
  bands   = bands[indx]
    
  retain  = [True if wave[band] >= lowave else False for band in bands]

  ##  Observed. mags. not affectded by IGM (rest-framce wave > 1500. A). 
  bands   = bands[retain]
  waves   = waves[retain]

  print('Available Bands redder than rest-frame 1500A.: {}'.format(bands))

  # print(waves)
  
  # Prepare the return.
  X2      = np.zeros(len(frame))
  f0s     = np.zeros(len(frame))
  betas   = np.zeros(len(frame))
  
  for i in np.arange(0, len(frame), 1):
    row            = frame.iloc[[i]]

    fluxes         = ['FLAMBDA_{}'.format(x.split('_')[-1]) for x in bands]
    flux_errs      = ['FLAMBDAERR_{}'.format(x.split('_')[-1]) for x in bands]

    # print(row)
    
    fluxes         = row[fluxes].values[0]
    flux_errs      = row[flux_errs].values[0]

    # print(waves)
    # print(fluxes)
    # print(flux_errs)

    res, f0, beta  = minimize_chi(chi2, waves, fluxes, flux_errs)

    X2[i]          = res
    f0s[i]         = f0
    betas[i]       = beta
    
  return  X2, f0s, betas, ''.join(x.split('_')[-1] for x in bands)

def test_set():
  waves          = np.array([3570.,        4766.,        6215.,        7544.,        8707.])                                                                                                                                        
  
  flambdas       = np.array([7.575258e-11, 2.120440e-10, 2.448201e-10, 2.041831e-10, 2.643931e-10])                                                                                                                                 
  flambda_errs   = np.array([3.835001e-12, 3.214882e-12, 2.238559e-12, 1.599317e-12, 1.376071e-12]) 

  return  waves, flambdas, flambda_errs

def test():
  waves, flambdas, flambda_errs = test_set()

  chi2, f0, beta = minimize_chi(chi2, waves, flambdas, flambda_errs)                                                                                                                                                                
  
  print(chi2, f0, beta)  

  pl.errorbar(waves, flambdas, yerr=flambda_errs)
  pl.plot(waves, continuum_model(waves, f0, beta), c='r')

  pl.savefig('plots/beta.pdf')

def run(redshift, wave, frame):
  frame['X2'], frame['f0'], frame['beta'], bands = solve_chi2s(chi2, redshift, wave, frame)

  return  frame, bands

def plot_betas():
  latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

  betas     = []
  beta_err  = []

  sbetas    = []
  sbeta_err = []
  
  redshifts = [2.024621, 3.003070, 3.963392, 5.0244]
  
  for i, (redshift, name) in enumerate(zip(redshifts, ['two', 'three', 'four', 'five'])):
    frame   = pd.read_pickle('bigdat/betauv_{:.3f}.pkl'.format(redshift).replace('.', 'p'))
    bands   = json.load(open('dat/beta_bands_{:.3f}.txt'.format(redshift), 'r'))

    ##  Read sample selection.                                                                                                                                                                                                      
    lsst_sample =  read_insample(redshift)
    isin        =  lsst_sample['INSAMPLE'].values

    print('Percentage meeting selection at z={}: {}'.format(redshift, np.mean(isin)))
    
    betas.append(np.mean(frame['beta'].values))
    beta_err.append(np.std(frame['beta'].values))

    sbetas.append(np.mean(frame['beta'].values[isin]))
    sbeta_err.append(np.std(frame['beta'].values[isin]))
    
  ##
  betas     = np.array(betas)
  beta_err  = np.array(beta_err)

  sbetas    = np.array(sbetas)
  sbeta_err = np.array(sbeta_err)
  
  pl.errorbar(redshifts,  betas,  beta_err, markersize=3, c=colors[0], alpha=0.5)
  pl.errorbar(redshifts, sbetas, sbeta_err, markersize=3, c=colors[0])
  
  pl.xlabel(r'$z$')
  pl.ylabel(r'$\beta_{\rm{UV}}$')

  plt.tight_layout()

  pl.savefig('plots/beta.pdf')

def plot_betas2d():
  latexify(columns=1, equal=True, fontsize=8, ggplot=True, usetex=True)

  for i, redshift in enumerate([2.024621, 3.003070, 3.963392, 5.0244]):
    frame = pd.read_pickle('bigdat/betauv_{:.3f}.pkl'.format(redshift).replace('.', 'p'))
    bands = json.load(open('dat/beta_bands_{:.3f}.txt'.format(redshift), 'r'))

    pl.plot(np.log10(frame['f0']), frame['beta'], marker='x', markersize=3, c=colors[i], lw=0, label=r'$' + bands + '$')

    pl.legend(loc=1, frameon=False)

    pl.xlim(-11., -7.)
    pl.ylim(-1.0, 0.0)
    
    pl.xlabel(r'$\log_{10}|f_0|$')
    pl.ylabel(r'$\beta$')
    
    plt.tight_layout()                                                                                                                                                                                                           

    pl.savefig('plots/beta_{:.3f}.pdf'.format(redshift))

  
if __name__ == '__main__':
    print('\n\nWelcome to the beta continuum.\n\n')

    boxsize        = 100.

    nrows          = -1
    compute        = True

    if compute:
      for redshift in [2.024621, 3.003070, 3.963392, 5.0244]:
        print('\n\nSolving for redshift: {}'.format(redshift))

        wave,  frame = get_pyloser_fluxes(boxsize, redshift, nrows=nrows)
        frame, bands = run(redshift, wave, frame)

        print(frame)

        with open('dat/beta_bands_{:.3f}.txt'.format(redshift), 'w') as ofile:
          json.dump(bands, ofile)
                
        frame.to_pickle('bigdat/betauv_{:.3f}.pkl'.format(redshift).replace('.', 'p'))

    else:
      plot_betas()

      ##  plot_betas2d()
        
    print('\n\nDone.\n\n')
