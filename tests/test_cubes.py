import copy

import numpy as np
from astropy import units as u

from bowshockpy.cube import BowshockCube
from bowshockpy.cubeproc import CubeProcessing
from bowshockpy.modelproj import ObsModel
from bowshockpy.models import BowshockModel

distpc = 300
L0 = (0.391 * distpc * u.au).to(u.km).value
zj = (4.58 * distpc * u.au).to(u.km).value
vj = 111.5
va = 0
v0 = 22.9
mass = 0.000231
rbf_obs = (0.75 * distpc * u.au).to(u.km).value
bsm = BowshockModel(
    L0=L0, zj=zj, vj=vj, va=va, v0=v0, mass=mass, distpc=distpc, rbf_obs=rbf_obs
)
bso = ObsModel(
    bsm,
    i_deg=20.0,
    vsys=0,
)
bsc1 = BowshockCube(
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
    cic=True,
    vt="2xchannel",
    tolfactor_vt=5,
    verbose=True,
)
bsc2 = copy.deepcopy(bsc1)
bsc1.makecube()
bsc3 = copy.deepcopy(bsc1)
bscp = CubeProcessing(
    [bsc1, bsc3],
    J=3,
    nu=345.79598990 * u.GHz,
    abund=8.5 * 10 ** (-5),
    meanmolmass=2.8,
    mu=0.112 * u.D,
    Tex=100 * u.K,
    Tbg=2.7 * u.K,
    coordcube="offset",
    bmin=0.1,
    bmaj=0.10,
    pabeam=-20.0,
    papv=bsc1.pa,
    parot=0,
    sigma_beforeconv=0.02,
    maxcube2noise=0,
)

bscp.calc_I()


def test_cube_mass_consistency():
    massconsistent = bsc1._check_mass_consistency()
    assert massconsistent, "Mass consistency check failed"


def test_makecube_fromcube():
    ones = np.ones_like(bsc1.cube)
    bsc2.makecube(fromcube=ones)
    massconsistent = bsc2._check_mass_consistency()
    assert (
        massconsistent
    ), "Mass consistency test failed while creating cube from an intial cube"


def test_combine_cubes():
    assert np.sum(bscp.cube) == np.sum(bsc1.cube) + np.sum(
        bsc3.cube
    ), "Mass consistency test failed while combining cubes"


def test_mass_index():
    # mass in control pixel
    mass_xyc = bscp.cube[36, 26, 25]
    assert np.isclose(
        mass_xyc, 2.8056152672732146e-07
    ), "CubeProcessing do not have the expected values of the masses"


def test_intensity_index():
    # intensity in control pixel
    intensity_xyc = bscp.cubes["I"][36, 26, 25]
    assert np.isclose(
        intensity_xyc, 0.02801070719875003
    ), f"CubeProcessing failed to obtain the expected values of the intensities"


def test_convolution():
    bscp.calc_I()
    bscp.convolve("I")
    assert np.isclose(
        np.sum(bscp.cubes["I"]), np.sum(bscp.cubes["I_c"])
    ), "Convolution failed: flux is not conserved"
