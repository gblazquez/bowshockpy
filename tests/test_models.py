
import numpy as np

from astropy import units as u

from bowshockpy import models as bs

distpc = 300
L0 = (0.391 * distpc * u.au).to(u.km).value
zj = (4.58 * distpc * u.au).to(u.km).value 
vj = 111.5                                    
va = 0                                      
v0 = 22.9                                    
mass = 0.000231                               
rbf_obs = (0.75 * distpc * u.au).to(u.km).value
bsm = bs.NarrowJet(
    L0=L0, zj=zj, vj=vj, va=va,
    v0=v0, mass=mass, distpc=distpc, rbf_obs=rbf_obs
    )
bso = bs.ObsModel(
    bsm,
    i=20.0*np.pi/180,
    vsys=0,
    )

def test_narrowjet():
    # Test that NarrowJet produce the expected values
    assert np.isclose(bsm.mp0_solmassyr, 1.4013290507098549e-06), \
        "NarrowJet failed to produce the expected value for mp0"
    assert np.isclose(bsm.rhoa_gcm3, 35.043736038867355e-19), \
        "NarrowJet failed to produce the expected value for rhoa"
    assert np.isclose(bsm.mpamb_f_solmassyr, 3.176808287706813e-06), \
        "NarrowJet failed to produce the expected value for mpamb_f"

def test_numerical_vs_analytical():
    # Test that NarrowJet produces the same values when some quantities are
    # computed analytically and numerically
    rhoa_num = bsm.rhoa_fromintmass_sigma_simple(0, bsm.rbf, bsm.mass)
    assert np.isclose(rhoa_num, bsm.rhoa), \
        "Analytical Ambient density computation differs from the numerical"
    mass_halfradius_analytical = bsm.intmass_analytical(bsm.rbf/2)
    mass_halfradius_numerical = bsm.intmass_numerical(0, bsm.rbf/2)
    assert np.isclose(mass_halfradius_analytical, mass_halfradius_numerical), \
        "Analytical mass computation differs from the numerical"