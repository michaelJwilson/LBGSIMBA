from collections import OrderedDict


def  get_depths():
    depths           = OrderedDict()
    
    ##  LSST-Y10                                                                                                                                  
    depths['LSST_u'] = 25.30
    depths['LSST_g'] = 26.84
    depths['LSST_r'] = 27.04
    depths['LSST_i'] = 26.35
    depths['LSST_z'] = 25.22
    depths['LSST_y'] = 24.47

    ##  Euclid                                                                                                                                    
    depths['Y'] = 24.0
    depths['J'] = 24.0
    depths['H'] = 24.0
    depths['K'] = 24.0

    ##                                                                                                                                            
    depths['I'] = 25.5
    depths['V'] = 25.5
    depths['G'] = 25.5
    depths['B'] = 25.5
    depths['R'] = 25.5
    depths['U'] = 25.5

    depths['ACS_F435W']  = 27.0
    depths['ACS_F606W']  = 27.0
    depths['ACS_F775W']  = 27.0
    depths['ACS_F850LP'] = 27.0

    return  depths
