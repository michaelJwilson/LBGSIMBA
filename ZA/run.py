from real_space_corr_fftw import CorrelationFunction


##  source activate ZA

# one_loop =False makes it Zeldovich
zelda = CorrelationFunction(cosmo.k,cosmo.p, one_loop= False, threads=1, jn=10, shear=True)

r, xi = zelda.combine_bias_terms(b1, b2, bs, counterterm)
