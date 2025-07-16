
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
bsc = bs.BowshockCube(
    bso,
    nphis=100,
    nzs=100,
    nc=50,
    vch0=-10, 
    vchf=-120,
    xpmax=5,    
    nxs=50,
    nys=50, 
    refpix=[25, 10], 
    CIC=True,
    vt="2xchannel",
    tolfactor_vt=5,
    verbose=True,
    )

def test_cube_mass_consistency():
    bsc.makecube()
    massconsistent = bsc._check_mass_consistency(return_isconsistent=True)
    assert massconsistent, "Mass consistency check failed"

def test_cube_concatenate():
    ones = np.ones_like(bsc.cube)
    bsc.makecube(fromcube=ones)
    massconsistent = bsc._check_mass_consistency(return_isconsistent=True)
    assert massconsistent, "Mass consistency in the while concatenating cubes failed"