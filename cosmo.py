import astropy.units as u

from   astropy.cosmology  import FlatLambdaCDM


# FlatLambdaCDM(H0=68., Om0=0.3)
cosmo = FlatLambdaCDM(H0=68, Om0=0.3, Tcmb0=2.73, Neff=3.04, m_nu=[0.,0.,0.] * u.eV, Ob0=0.048)
