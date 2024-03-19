import numpy as np

import astropy.units as u
import astropy.constants as c

import bowshock_params as p

from bowpy import bsmodels as bs
from bowpy import bsutils as bu

ps = { 
 "modelname": p.modelname,
 'rj': (p.rj * p.distpc * u.au).to(u.km).value,
 'L0': (p.L0 * p.distpc * u.au).to(u.km).value,   
 'zj': (p.zj * p.distpc * u.au).to(u.km).value,
 'vj': p.vj,
 'vw': p.vw,   
 'v0': p.v0,
 'rbf_obs': None, # TODO: (p.rbf_obs * p.distpc * u.au).to(u.km).value,
 'mass': p.mass,   
}

bsm = bs.NJ(ps)

psobs = { 
 'i': p.i * np.pi / 180,
 'vsys': p.vsys,
 'distpc': p.distpc,
 "nzs": p.nzs,
}

bsmobs = bs.ObsModel(ps, psobs)
if p.bs2Dplot:
    bs2Dplot = bs.Bowshock2DPlots(ps, psobs)
    bu.make_folder(f"models/{ps['modelname']}")
    bs2Dplot.fig_model.savefig(f"models/{ps['modelname']}/2D.pdf")

if len(p.outcubes) != 0:
    pscube = {
        "nphis": p.nphis,   
        "nc": p.nc,
        "vt": p.vt,
        "vch0": p.vch0,
        "vchf": p.vchf,   
        "nxs": p.nxs,
        "nys": p.nys,   
        "refpix": p.refpix, 
        "xpmax": p.xpmax,    
        "pa": p.pa,
        "ybeam": p.ybeam,
        "xbeam": p.xbeam,
        "pabeam": p.pabeam,
        "CIC": p.CIC,
        "tolfactor_vt": p.tolfactor_vt,
    }
    
    pscube["chanwidth"] = (pscube["vchf"] - pscube["vch0"]) / (pscube["nc"]-1)
    pscube["abschanwidth"] = np.abs(pscube["chanwidth"])
    pscube["vt"] = pscube["vt"] if type(pscube["vt"])!=str \
          else float(pscube["vt"].split("x")[0])*pscube["chanwidth"]
    pscube["arcsecpix"] = pscube["xpmax"] / float(pscube["nxs"])
    pscube["x_FWHM"] = pscube["xbeam"] / pscube["arcsecpix"]
    pscube["y_FWHM"] = pscube["ybeam"] / pscube["arcsecpix"]
    pscube["beamarea"] = np.pi * pscube["y_FWHM"] * pscube["x_FWHM"] / (4 * np.log(2))
    if pscube["refpix"] == None:
        if pscube["nxs"]%2 == 0:
            xref = pscube["nxs"] / 2
        else: 
            xref = (pscube["nxs"]-1) / 2
        if pscube["nys"]%2 == 0:
            yref = pscube["nys"] / 2
        else: 
            yref = (pscube["nys"]-1) / 2       
        pscube["refpix"] = [xref, yref]

    bsc = bs.BowshockCube(ps, psobs, pscube)
    bsc.makecube()

    mpars = {
        "muH2": p.muH2,
        "XCO": p.XCO,
        "meanmass": p.muH2 / (6.023*10**23) * u.g,
        "Tex": p.Tex*u.K,
        "Tbg": p.Tbg*u.K,
        "ra_source_deg": p.ra_source_deg,
        "dec_source_deg": p.dec_source_deg
    }

    bscs = bs.CubeProcessing(bsc, mpars)
    bscs.calc(["I_rc"])
    bscs.savecubes(["m", "tau", "I_rc", "NCO"])