from astropy.cosmology import Planck18_arXiv_v2
import astropy.units as u
import numpy as np

s_per_solar_mass = 4.93e-6

z_ = np.expm1(np.linspace(np.log1p(1), np.log1p(100), 1024))
Rz_ = 4*np.pi*Planck18_arXiv_v2.differential_comoving_volume(z_).to(u.Gpc**3/u.sr).value/(1+z_)*(1+z_)**2.7/(1 + ((1+z_)/(1+1.9))**5.6)
R1 = np.trapz(Rz_, z_) # 1/yr
T1 = 3e7/R1 # s

def fdot_of_Mc(Mc, f):
    """Returns df/dt (1/s^2) for the given chirp mass and frequency."""
    return (Mc*s_per_solar_mass)**(5/3)*96/5*np.pi**(8/3)*f**(11/3)

def overlapping_fom(Mc, f, R):
    """..."""
    return fdot_of_Mc(Mc, f)*(T1/R)**2