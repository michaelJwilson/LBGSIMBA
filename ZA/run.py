import  pylab                as      pl
import  numpy                as      np

from    real_space_corr_fftw import  CorrelationFunction


# source activate ZA

k, Pk, PNL = np.loadtxt('dat/pklin_z502.txt', unpack=True)

# one_loop =False makes it Zeldovich
zelda      = CorrelationFunction(k, Pk, one_loop= False, threads=1, jn=10, shear=True)

# b1, b2, bs, counter-term.
r, xi      = zelda.combine_bias_terms(0.5, 0.5, 0.5, 0.0)

pl.loglog(r, xi)

pl.show()
