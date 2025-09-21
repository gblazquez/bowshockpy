import copy

import numpy as np
from astropy import units as u

from bowshockpy.cubemass import MassCube
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
    L0=L0,
    zj=zj,
    vj=vj,
    va=va,
    v0=v0,
    mass=mass,
    distpc=distpc,
    rbf_obs=rbf_obs,
)
bso = ObsModel(
    bsm,
    i_deg=20.0,
    vsys=0,
)
bsc1 = MassCube(
    bso,
    nphis=100,
    xpmax=5,
    vch0=-10,
    vchf=-120,
    chanwidth=None,
    nzs=100,
    nc=50,
    nxs=50,
    nys=50,
    refpix=[25, 10],
    cic=True,
    vt="2xchannel",
    tolfactor_vt=5,
    verbose=True,
)
bsc1.makecube()

bsc2 = MassCube(
    bso,
    nphis=100,
    xpmax=5,
    vch0=-10,
    vchf=None,
    chanwidth=bsc1.chanwidth,
    nzs=100,
    nc=50,
    nxs=50,
    nys=50,
    refpix=[25, 10],
    cic=True,
    vt=np.abs(bsc1.vt),
    tolfactor_vt=5,
    verbose=True,
)

bsc2.makecube()

bscp = CubeProcessing(
    [bsc1, bsc2],
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
    sigma_beforeconv=0.02,
    maxcube2noise=0,
)

bscp.calc_I()

mom0 = bscp.mom0(ck="I")
mom1 = bscp.mom1(ck="I")
mom2 = bscp.mom2(ck="I")
maxintens = bscp.maxintens(ck="I")
pv = bscp.pvalongz(ck="I", halfwidth=3)


def test_channel_consistency_vchf():
    assert np.isclose(bsc1.vchf, bsc1.velchans[-1])


def test_channel_consistency_chanwidth():
    assert np.isclose(bsc1.chanwidth, np.diff(bsc1.velchans)[0])


def test_cube_mass_consistency():
    massconsistent = bsc1._check_mass_consistency(tol=None)
    assert massconsistent, "Mass consistency check failed"


def test_makecube_fromcube():
    ones = np.ones_like(bsc1.cube)
    bsc_test = copy.deepcopy(bsc1)
    bsc_test.makecube(fromcube=ones)
    massconsistent = bsc_test._check_mass_consistency(tol=None)
    assert (
        massconsistent
    ), "Mass consistency test failed while creating cube from an intial cube"


def test_combine_cubes():
    masscombined = np.sum(bscp.cube)
    mass1 = np.sum(bsc1.cube)
    mass2 = np.sum(bsc2.cube)
    masstot = mass1 + mass2
    assert np.isclose(
        masscombined, masstot
    ),  f"{masscombined, masstot}"

    # f"Mass consistency test failed while combining cubes {masscombined, masstot}"


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
    ), "CubeProcessing failed to obtain the expected values of the intensities"


def test_convolution():
    bscp.convolve("I")
    assert np.isclose(
        np.sum(bscp.cubes["I"]), np.sum(bscp.cubes["I_c"])
    ), "Convolution failed: flux is not conserved"


def test_mom0():
    mom0_xy = mom0[26, 25]
    assert np.isclose(
        mom0_xy, 1.715026191836612
    ), "Failed to obtain the expected values of the moment 0"


def test_mom1():
    mom1_xy = mom1[26, 25]
    assert np.isclose(
        mom1_xy, -86.69768758752018
    ), "Failed to obtain the expected values of the moment 1"


def test_mom2():
    mom2_xy = mom2[26, 25]
    assert np.isclose(
        mom2_xy, 27.87445981444194
    ), "Fail to obtain the expected vaules of moment 2"


def test_maxintens():
    maxintens_xy = maxintens[26, 25]
    assert np.isclose(
        maxintens_xy, 0.06893409211230457
    ), "Fail to obtain the expected values of the max. intensity moment"


def test_pv():
    pv_xy = pv[13, 26]
    assert np.isclose(
        pv_xy, 0.03828409594116648
    ), "Fail to obtain the expected values of position-velocity diagram"


def test_cube_header():
    ck = "I"
    _ , ny, nx = np.shape(bscp.cubes[ck])
    hdrcube = bscp._create_header_cube(ck)
    hdrpv = bscp._create_header_pv(ck, bunit="")
    hdrmom = bscp._create_header_moment(ck, bunit="")
    ra_i_edge_py = (-0.5 - bscp.refpixs[ck][0]) * bscp.arcsecpix
    ra_f_edge_py = (nx - 0.5 - bscp.refpixs[ck][0]) * bscp.arcsecpix
    dec_i_edge_py = (-0.5 - bscp.refpixs[ck][1]) * bscp.arcsecpix
    dec_f_edge_py = (ny - 0.5 - bscp.refpixs[ck][1]) * bscp.arcsecpix
    vel_i_edge_py = bscp.velchans[0] - bscp.chanwidth / 2
    vel_f_edge_py = bscp.velchans[-1] + bscp.chanwidth / 2

    ra_i_edge_fits = (
        hdrcube["CDELT1"] * (0.5 - hdrcube["CRPIX1"]) + hdrcube["CRVAL1"]
    )
    ra_f_edge_fits = (
        hdrcube["CDELT1"] * (hdrcube["NAXIS1"] + 0.5 - hdrcube["CRPIX1"])
        + hdrcube["CRVAL1"]
    )

    dec_i_edge_fits = (
        hdrcube["CDELT2"] * (0.5 - hdrcube["CRPIX2"]) + hdrcube["CRVAL2"]
    )
    dec_f_edge_fits = (
        hdrcube["CDELT2"] * (hdrcube["NAXIS2"] + 0.5 - hdrcube["CRPIX2"])
        + hdrcube["CRVAL2"]
    )

    ra_i_edge_momfits = (
        hdrmom["CDELT1"] * (0.5 - hdrmom["CRPIX1"]) + hdrmom["CRVAL1"]
    )
    ra_f_edge_momfits = (
        hdrmom["CDELT1"] * (hdrmom["NAXIS1"] + 0.5 - hdrmom["CRPIX1"])
        + hdrmom["CRVAL1"]
    )

    dec_i_edge_momfits = (
        hdrmom["CDELT2"] * (0.5 - hdrmom["CRPIX2"]) + hdrmom["CRVAL2"]
    )
    dec_f_edge_momfits = (
        hdrmom["CDELT2"] * (hdrmom["NAXIS2"] + 0.5 - hdrmom["CRPIX2"])
        + hdrmom["CRVAL2"]
    )

    vel_i_edge_fits = (
        hdrcube["CDELT3"] * (0.5 - hdrcube["CRPIX3"]) + hdrcube["CRVAL3"]
    )

    vel_f_edge_fits = (
        hdrcube["CDELT3"] * (hdrcube["NAXIS3"] + 0.5 - hdrcube["CRPIX3"])
        + hdrcube["CRVAL3"]
    )

    vel_i_edge_pvfits = (
        hdrpv["CDELT2"] * (0.5 - hdrpv["CRPIX2"]) + hdrpv["CRVAL2"]
    )

    vel_f_edge_pvfits = (
        hdrpv["CDELT2"] * (hdrpv["NAXIS2"] + 0.5 - hdrpv["CRPIX2"])
        + hdrpv["CRVAL2"]
    )

    ra_i_ok = np.isclose(ra_i_edge_py, ra_i_edge_fits)
    ra_f_ok = np.isclose(ra_f_edge_py, ra_f_edge_fits)

    dec_i_ok = np.isclose(dec_i_edge_py, dec_i_edge_fits)
    dec_f_ok = np.isclose(dec_f_edge_py, dec_f_edge_fits)

    vel_i_ok = np.isclose(vel_i_edge_py, vel_i_edge_fits)
    vel_f_ok = np.isclose(vel_f_edge_py, vel_f_edge_fits)

    vel_i_ok_pv = np.isclose(vel_i_edge_py, vel_i_edge_pvfits)
    vel_f_ok_pv = np.isclose(vel_f_edge_py, vel_f_edge_pvfits)

    ra_i_momok = np.isclose(ra_i_edge_py, ra_i_edge_momfits)
    ra_f_momok = np.isclose(ra_f_edge_py, ra_f_edge_momfits)

    dec_i_momok = np.isclose(dec_i_edge_py, dec_i_edge_momfits)
    dec_f_momok = np.isclose(dec_f_edge_py, dec_f_edge_momfits)

    assert all(
        [
            ra_i_ok,
            ra_f_ok,
            dec_i_ok,
            dec_f_ok,
            vel_i_ok,
            vel_f_ok,
            vel_i_ok_pv,
            vel_f_ok_pv,
            ra_i_momok,
            ra_f_momok,
            dec_i_momok,
            dec_f_momok,
        ]
    ), f"{ra_i_ok, ra_f_ok, dec_i_ok, dec_f_ok}"
