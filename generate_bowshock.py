
import numpy as np

import astropy.units as u
import astropy.constants as c

import os

import importlib

bpf_str = input(
    """
    Enter the name of the file where the parameters are defined (default: bowshock_params): 
    """
)
bpf_str = bpf_str if bpf_str is not "" else "bowshock_params"
p = importlib.import_module(bpf_str.strip(".py"))

from bowpy import bsmodels as bs
from bowpy import bsutils as bu

pss = []
for i in range(p.nbowshocks):
    pss += [{
     "modelname": p.modelname,
     'rj':      (p.__getattribute__(f"rj_{i+1}") * p.distpc * u.au).to(u.km).value,
     'L0':      (p.__getattribute__(f"L0_{i+1}") * p.distpc * u.au).to(u.km).value,
     'zj':      (p.__getattribute__(f"zj_{i+1}") * p.distpc * u.au).to(u.km).value,
     'vj':       p.__getattribute__(f"vj_{i+1}"),
     'vw':       p.__getattribute__(f"vw_{i+1}"),
     'v0':       p.__getattribute__(f"v0_{i+1}"),
     'rbf_obs': (p.__getattribute__(f"rbf_obs_{i+1}") * p.distpc * u.au).to(u.km).value 
         if p.__getattribute__(f"rbf_obs_{i+1}") is not None 
         else p.__getattribute__(f"rbf_obs_{i+1}"),
     'mass':     p.__getattribute__(f"mass_{i+1}"),
    }]

psobs = {
 'i': p.i * np.pi / 180,
 'vsys': p.vsys,
 'distpc': p.distpc,
 "nzs": p.nzs,
}

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
        "bmaj": p.bmaj,
        "bmin": p.bmin,
        "pabeam": p.pabeam,
        "CIC": p.CIC,
        "tolfactor_vt": p.tolfactor_vt,
        "maxcube2noise": p.maxcube2noise,
        "verbose": p.verbose,
    }
    pscube["chanwidth"] = (pscube["vchf"] - pscube["vch0"]) / (pscube["nc"]-1)
    pscube["abschanwidth"] = np.abs(pscube["chanwidth"])
    pscube["vt"] = pscube["vt"] if type(pscube["vt"])!=str \
          else float(pscube["vt"].split("x")[0])*pscube["chanwidth"]
    pscube["arcsecpix"] = pscube["xpmax"] / float(pscube["nxs"])
    pscube["x_FWHM"] = pscube["bmin"] / pscube["arcsecpix"]
    pscube["y_FWHM"] = pscube["bmaj"] / pscube["arcsecpix"]
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
    mpars = {
        "muH2": p.muH2,
        "XCO": p.XCO,
        "meanmass": p.muH2 / (6.023*10**23) * u.g,
        "Tex": p.Tex*u.K,
        "Tbg": p.Tbg*u.K,
        "ra_source_deg": p.ra_source_deg,
        "dec_source_deg": p.dec_source_deg
    }

bscs = []
for i, ps in enumerate(pss):
    bsm = bs.NJ(ps)
    bsmobs = bs.ObsModel(ps, psobs)
    if p.bs2Dplot:
        bs2Dplot = bs.Bowshock2DPlots(ps, psobs)
        if i == 0:
            bu.make_folder(f"models/{ps['modelname']}")
        bs2Dplot.fig_model.savefig(f"models/{ps['modelname']}/2D_{i+1}.pdf")
    
    if len(p.outcubes) != 0:
        print(f"""
Generating bowshock {i+1}/{p.nbowshocks}
              """)
        if i == 0:
            bscs += [bs.BowshockCube(ps, psobs, pscube)]
            bscs[i].makecube()
        else:
            bscs += [bs.BowshockCube(ps, psobs, pscube)]
#            import pdb; pdb.set_trace()
            bscs[i].makecube(fromcube=bscs[i-1].cube)

bscp = bs.CubeProcessing(bscs[-1], mpars)
bscp.calc(p.outcubes)
bscp.savecubes(p.outcubes)

# Save the file with all the parameters used to generate the bowshocks
os.system("cp {bpf_str.strip('.py')}.py models/{p.modelname}")