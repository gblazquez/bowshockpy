import numpy as np

import astropy.units as u

import os

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

from bowshockpy import models as bs
from bowshockpy import utils as ut
from bowshockpy.version import __version__


def generate_bowshock(p):
    print(
    f"""

--------------------------------------------
BowshockPy v{__version__}

https://bowshockpy.readthedocs.io/en/latest/
--------------------------------------------

Parameters read from {p.filename}
    """
    )
    pss = []
    psobss = []
    for i in range(p.nbowshocks):
        pss += [{
         "modelname": p.modelname,
         'L0':      (p.__getattribute__(f"L0_{i+1}") * p.distpc * u.au).to(u.km).value,
         'zj':      (p.__getattribute__(f"zj_{i+1}") * p.distpc * u.au).to(u.km).value,
         'vj':       p.__getattribute__(f"vj_{i+1}"),
         'va':       p.__getattribute__(f"va_{i+1}"),
         'v0':       p.__getattribute__(f"v0_{i+1}"),
         'rbf_obs': (p.__getattribute__(f"rbf_obs_{i+1}") * p.distpc * u.au).to(u.km).value
             if p.__getattribute__(f"rbf_obs_{i+1}") is not None
             else p.__getattribute__(f"rbf_obs_{i+1}"),
         'mass':     p.__getattribute__(f"mass_{i+1}"),
         'pa':       p.__getattribute__(f"pa_{i+1}") * np.pi / 180,
        }]

        psobss += [{
         'i': p.__getattribute__(f"i_{i+1}") * np.pi / 180,
         'vsys': p.vsys,
         'distpc': p.distpc,
         "nzs": p.nzs,
        }]

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
            "parot": p.parot,
            "papv": p.papv,
            "bmaj": p.bmaj,
            "bmin": p.bmin,
            "pabeam": p.pabeam,
            "CIC": p.CIC,
            "tolfactor_vt": p.tolfactor_vt,
            "sigma_beforeconv": p.sigma_beforeconv,
            "maxcube2noise": p.maxcube2noise,
            "verbose": p.verbose,
        }
        mpars = {
            "muH2": p.muH2,
            "J": p.J,
            "XCO": p.XCO,
            "meanmass": p.muH2 / (6.023*10**23) * u.g,
            "Tex": p.Tex*u.K,
            "Tbg": p.Tbg*u.K,
            "ra_source_deg": p.ra_source_deg,
            "dec_source_deg": p.dec_source_deg,
            "coordcube": p.coordcube
        }

    bscs = []
    for i, (ps,psobs) in enumerate(zip(pss,psobss)):
        bsm = bs.NarrowJet(
            L0=ps["L0"],
            zj=ps["zj"],
            vj=ps["vj"],
            va=ps["va"],
            v0=ps["v0"],
            mass=ps["mass"],
            distpc=psobs["distpc"],
            rbf_obs=ps["rbf_obs"]
            )
        bsmobs = bs.ObsModel(
            model=bsm,
            i=psobs["i"],
            vsys=psobs["vsys"],
            )
        if p.bs2Dplot:
            plt_model = bsm.modelplot(
                modelname=ps["modelname"],
            )
            plt_model.plot()
            if i == 0:
                ut.make_folder(f"models/{ps['modelname']}")
            plt_model.savefig(
                f"models/{ps['modelname']}/bowshock_model_{i+1}.pdf",
                )
            plt_obsmodel = bsmobs.modelplot(
                )
            plt_obsmodel.plot()
            plt_obsmodel.savefig(
                f"models/{ps['modelname']}/bowshock_projected_{i+1}.jpg",
                dpi=300,
            )
            # bs2Dplot = bs.Bowshock2DPlots(bsmobs, ps['modelname'])
            # if i == 0:
            #     ut.make_folder(f"models/{ps['modelname']}")
            # bs2Dplot.savefig(
            #     f"models/{ps['modelname']}/2D_{i+1}.pdf",
            #     )
        if len(p.outcubes) != 0:
            print(f"""

Generating bowshock {i+1}/{p.nbowshocks}
                  """)
            if i == 0:
                bscs += [
                    bs.BowshockCube(
                        obsmodel=bsmobs,
                        nphis=pscube["nphis"],
                        vch0=pscube["vch0"],
                        vchf=pscube["vchf"],
                        xpmax=pscube["xpmax"],
                        nc=pscube["nc"],
                        nxs=pscube["nxs"],
                        nys=pscube["nys"],
                        refpix=pscube["refpix"],
                        CIC=pscube["CIC"],
                        vt=pscube["vt"],
                        tolfactor_vt=pscube["tolfactor_vt"],)
                    ]
                bscs[i].makecube()
                print(f"""
Channel width: {bscs[i].abschanwidth:.3} km/s
Pixel size: {bscs[i].arcsecpix:.4} arcsec/pix
     """)
            else:
                bscs += [
                    bs.BowshockCube(
                        obsmode=bsmobs,
                        nphis=pscube["nphis"],
                        nc=pscube["nc"],
                        vch0=pscube["vch0"],
                        vchf=pscube["vchf"],
                        xpmax=pscube["xpmax"],
                        nxs=pscube["nxs"],
                        nys=pscube["nys"],
                        refpix=pscube["refpix"],
                        CIC=pscube["CIC"],
                        vt=pscube["vt"],
                        tolfactor_vt=pscube["tolfactor_vt"],)
                    ]
                bscs[i].makecube(fromcube=bscs[i-1].cube)
    bscs[-1].plot_channels(
        vmax=np.percentile(bscs[-1], 99.9),
        savefig=f"models/{ps['modelname']}/bowshock_cube_{i+1}.pdf",
    )

    print(
f"""
The masses has been computed!

The cubes are going to be processed in order to get the desired outputs
specified in {p.filename}. The outputs will be saved in fits format. The
filename of each cube indicate its quantity and the operations applied to the
cube ("<quantity>_<operations>.fits"). Some abbreviations will be used:

Abbreviations for quantities are:        Abbreviations for the operations are:
    m: mass [SolarMass]                      s: add_source
    I: Intensity [Jy/beam]                   r: rotate
    Ithin: Intensity [Jy/beam]               n: add_noise
    NCO: CO column density [cm-2]            c: convolve
    tau: Opacity
"""
    )
    bscp = bs.CubeProcessing(
        bscs[-1],
        modelname=ps["modelname"],
        J=mpars["J"],
        XCO=mpars["XCO"],
        meanmass=mpars["meanmass"],
        Tex=mpars["Tex"],
        Tbg=mpars["Tbg"],
        coordcube=mpars["coordcube"],
        ra_source_deg=mpars["ra_source_deg"],
        dec_source_deg=mpars["dec_source_deg"],
        bmin=pscube["bmin"],
        bmaj=pscube["bmaj"],
        pabeam=pscube["pabeam"],
        papv=pscube["papv"],
        parot=pscube["parot"],
        sigma_beforeconv=pscube["sigma_beforeconv"],
        maxcube2noise=pscube["maxcube2noise"],
            )
    bscp.calc(p.outcubes)
    bscp.savecubes(p.outcubes)
    bscp.momentsandpv_and_params_all(
         bscs,
         savefits=p.savefits,
         saveplot=p.saveplot,
         mom1clipping=p.mom1clipping,
         mom2clipping=p.mom2clipping,
         mom0values=p.mom0values,
         mom1values=p.mom1values,
         mom2values=p.mom2values,
         mom8values=p.mom8values,
         pvvalues=p.pvvalues,
         )

    # Save the file with all the parameters used to generate the bowshocks
    os.system(f"cp {p.filename.rstrip('.py')}.py models/{p.modelname}")


def main():
    # bpf_str = input(
    # """
    # Enter the name of the file where the parameters are defined (default: bowshock_params):
    # """
    # )
    import argparse
    import runpy

    description = """
Bowshockpy is a Python package that generates synthetic spectral cubes,
position-velocity diagrams, and moment images for a simple analytical jet-driven
bowshock model, using the prescription for protostellar jets presented in
Ostriker et al. (2001) and Tabone et al. (2018). Please, see the documentation
at:

https://bowshockpy.readthedocs.io/en/latest/

    """

    parser = argparse.ArgumentParser(
        description=description
    )
    parser.add_argument(
        "-r", "--read",
        dest="parameters_file",
        type=str,
        help="Reads a configuration file to generate the bowshock model",
        default="None"
        )
    parser.add_argument(
        "-p", "--print-example",
        dest="inputfile_example",
        type=str,
        help="""
        Prints an example of input file. Write the number of the corresponding
        example that is closer to your needs. There are 3 examples: write 1 to
        print an example of input file of a redshifted bowshock model, write 2
        for a model including two redshifted bowshocks, write 3 for a
        blueshifted bowshock. See https://bowshockpy.readthedocs.io/en/latest/
        for a detailed documentation of the examples.  """,
        default="None"
        )

    args = parser.parse_args()
    filename = args.parameters_file
    nexample = args.inputfile_example
    if filename != "None":
        parameters = runpy.run_path(filename)
        p = ut.VarsInParamFile(parameters)
        generate_bowshock(p)
    if nexample != "None":
        ut.print_example(nexample)

if __name__ == "__main__":
    main()
