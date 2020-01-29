import numpy as np
import pylab as pl


def limber():
    # r[Mpc/h]     r^2.xiL       xir:1      xir:b1      xir:b2    xir:b1^2   xir:b1.b2    xir:b2^2
    
    for redshift in [2.024621, 3.00307, 3.963392, 5.0244]:
        # ZA                                                                                                                                                                                                                                                                              
        iz     = int(100 * redshift + 0.001)
        _      = np.loadtxt('/home/mjwilson/LBGSIMBA/dat/white/zeld_z{}.txt'.format(iz))

	b1, b2 = 1.0, 0.0

        cc     = np.array([0., 0., 1.0, b1, b2, b1**2, b1*b2, b2**2])

        rs     = _[:,0]
        _      = _[:,0:8]

        # Zeldovich correlation fn. 
        result = np.dot(_, cc)

    

if __name__ == '__main__':
    
    limber()
    
    print('\n\nDone.\n\n')
