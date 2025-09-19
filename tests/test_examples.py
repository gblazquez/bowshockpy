import bowshockpy.inputfiles.example1 as e1
import bowshockpy.inputfiles.example2 as e2
import bowshockpy.inputfiles.example3 as e3
import bowshockpy.inputfiles.example4 as e4
import bowshockpy.inputfiles.example5 as e5

list_required_params = [
    "modelname",
    "L0_1",
    "zj_1",
    "vj_1",
    "va_1",
    "v0_1",
    "rbf_obs_1",
    "mass_1",
    "i_1",
    "pa_1",
    "vsys",
    "distpc",
    "nphis",
    "nc",
    "vt",
    "vch0",
    "vchf",
    "chanwidth",
    "nxs",
    "nys",
    "nzs",
    "refpix",
    "xpmax",
    "papv",
    "bmaj",
    "bmin",
    "pabeam",
    "cic",
    "tolfactor_vt",
    "sigma_beforeconv",
    "maxcube2noise",
    "verbose",
    "J",
    "nu",
    "abund",
    "mu",
    "meanmolmass",
    "Tex",
    "Tbg",
    "ra_source_deg",
    "dec_source_deg",
    "coordcube",
]


def test_params_examples():
    for param in list_required_params:
        getattr(e1, param)
        getattr(e2, param)
        getattr(e3, param)
        getattr(e4, param)
        getattr(e5, param)
