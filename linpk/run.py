import matplotlib;  matplotlib.use('Agg')

import numpy             as      np
import matplotlib.pyplot as      plt

from   classy            import  Class


# create instance of the class "Class"
LambdaCDM = Class()

# Your value of omega_b=4.800000e-02 is out of the bounds [5.000000e-03 , 3.900000e-02]
# Planck (2015):  https://arxiv.org/pdf/1502.01589.pdf
# SIMBA cosmology:  https://arxiv.org/pdf/1901.10203.pdf
# 'YHe': 0.252

h          = 0.68
Ob         = 0.048
Om         = 0.3
Ocdm       = Om - Ob
s8         = 0.82
ns         = 0.97

simba      = {'omega_b': Ob * h * h, 'omega_cdm': Ocdm * h * h, 'h': h, 'n_s': ns, 'tau_reio': 0.0925, 'A_s': 2.215e-9}

# planck2015 = {'omega_b': 0.02222, 'omega_cdm': 0.1197, 'h': 0.68, 'A_s': 2.215e-9, 'n_s':0.9619, 'tau_reio': 0.0925}

params     = simba

# pass input parameters
LambdaCDM.set(params)

LambdaCDM.set({'output':'tCl,pCl,lCl,mPk','lensing':'yes','P_k_max_1/Mpc':3.0})

# run class
LambdaCDM.compute()

derived = LambdaCDM.get_current_derived_parameters(['sigma8'])


#
kk = np.logspace(-4,np.log10(3),1000) # k [h/Mpc].

Pk = []            # P(k) in (Mpc/h)**3
h  = LambdaCDM.h() # get reduced Hubble for conversions to 1/Mpc

for k in kk:
  Pk.append(LambdaCDM.pk(k*h,0.)*h**3) # function .pk(k,z)

Pk = np.array(Pk)
  
print(derived['sigma8'])
print(s8)

Pk = Pk * ( s8 / derived['sigma8'] )**2

del simba['A_s']

simba['sigma8'] = s8 
pstring         = ''.join(['%s:  %s;  ' % (key, value) for (key, value) in params.items()]).strip()

# Save.
np.savetxt('linpk.txt', np.c_[kk, Pk], fmt='%.6le', header=pstring)
    
# plot P(k)
plt.figure(2)

plt.xscale('log')
plt.yscale('log')

plt.xlim(kk[0],kk[-1])

plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')

plt.ylabel(r'$P(k) \,\,\,\, [\mathrm{Mpc}/h]^3$')

plt.plot(kk, Pk, 'k-')

plt.savefig('plots/linpk.pdf')

# optional: clear content of LambdaCDM (to reuse it for another model)
LambdaCDM.struct_cleanup()

# optional: reset parameters to default
LambdaCDM.empty()
