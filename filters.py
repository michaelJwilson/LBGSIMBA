from  collections import OrderedDict


filters = OrderedDict()


filters[0]  = 'V'
filters[1]  = 'SDSS-u'
filters[2]  = 'SDSS-g'
filters[3]  = 'SDSS-r'
filters[4]  = 'SDSS-i'
filters[5]  = 'SDSS-z'
filters[6]  = 'U'
filters[7]  = '2MASS-J'
filters[8]  = '2MASS-H'
filters[9]  = '2MASS-K'
filters[10] = 'F275W'
filters[11] = 'F475W'
filters[12] = 'F606W'
filters[13] = 'F814W'
filters[14] = 'F105W'
filters[15] = 'F140W'
filters[16] = 'F160W'
filters[17] = 'IRAC1'
filters[18] = 'IRAC2'
filters[19] = 'MUV'
filters[20] = 'JWST-F090W'
filters[21] = 'JWST-F115W'
filters[22] = 'JWST-F150W'
filters[23] = 'JWST-F200W'
filters[24] = 'NUV'
filters[25] = 'Y'
filters[26] = 'J'
filters[27] = 'H'
filters[28] = 'u'
filters[29] = 'g'
filters[30] = 'r'
filters[31] = 'i'
filters[32] = 'y'
filters[33] = 'z'
filters[34] = 'Subaru-B'
filters[35] = 'Subaru-V'
filters[36] = 'Subaru-R'
filters[37] = 'Subaru-i'
filters[38] = 'Subaru-z'
filters[39] = 'VISTA-Y'
filters[40] = 'WFCAM-J'
filters[41] = 'WFCAM-H'
filters[42] = 'WFCAM-K'


if __name__ == '__main__':
    print(filters.values())
