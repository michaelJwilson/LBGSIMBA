#!/usr/bin/env python3
#
# Python code to compute the linear theory P(k) using CLASS.
#
import numpy  as      np
from   classy import Class


def make_pk(zlist):
    """
    Creates the P(k) files at redshifts in zlist.
    """
    OmM = 0.3
    OmB = 0.048
    hub = 0.68
    
    params = {
    'output': 'mPk',
    'P_k_max_h/Mpc': 150.,
    'z_pk': '0,10',
    'A_s': 2.1353046867600372e-09,
    'n_s': 0.97,
    'h': hub,
    'N_ur': 3.0328,
    'tau_reio': 0.06,
    'omega_b': OmB*hub**2,
    'omega_cdm': (OmM-OmB)*hub**2
    }  
    
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    
    # now put in a flag for HaloFit power spectra.
    params['non linear'] = 'halofit'
    nonlin = Class()
    nonlin.set(params)
    nonlin.compute()

    # First generate D(z) and f(z)
    fout = open("growth.txt","w")
    fout.write("# Growth factors vs. redshift.\n")
    fout.write("# OmM={:f}, OmB={:f}, hub={:f}\n".format(OmM,OmB,hub))
    fout.write("# ns={:f}, A_s={:e}, sig8={:f}\n".\
               format(params['n_s'],params['A_s'],cosmo.sigma8()))
    fout.write("# {:>3s} {:>8s} {:>8s}\n".format("z","D(z)","f(z)"))

    for zz in np.arange(0.0,9.51,0.25):
        fout.write("{:5.2f} {:8.4f} {:8.4f}\n".\
                  format(zz,cosmo.scale_independent_growth_factor(zz),\
                            cosmo.scale_independent_growth_factor_f(zz)))
    fout.close()

    # Now compute P(k) at each zz and write each to a file.
    kk = np.logspace(-4.0,2.0,300)

    for zz in zlist:
        iz = int(100*zz+0.001)
        pk = np.array([ cosmo.pk(k*params['h'],zz)*params['h']**3 for k in kk])
        hf = np.array([nonlin.pk(k*params['h'],zz)*params['h']**3 for k in kk])
        fout = open("pklin_z{:03d}.txt".format(iz),"w")
        fout.write("# Matter power spectra at z={:f}.\n".format(zz))
        fout.write("# D(z)={:f}, f(z)={:f}.\n".\
                   format(cosmo.scale_independent_growth_factor(zz),\
                          cosmo.scale_independent_growth_factor_f(zz)))
        fout.write("# {:>13s} {:>15s} {:>15s}\n".\
                   format("k[h/Mpc]","Plin(k)","HaloFit"))

        for i in range(kk.size):
            fout.write("{:15.5e} {:15.5e} {:15.5e}\n".format(kk[i],pk[i],hf[i]))
        fout.close()


if __name__=="__main__":
    zlist = np.arange(0.0, 5.0, 0.1)
    make_pk(zlist)
    
