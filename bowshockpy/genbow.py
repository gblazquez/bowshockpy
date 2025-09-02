"""This module contains all the workflow of BowshockPy when it is run from the
terminal"""

import os

import astropy.units as u

import argparse
import runpy

from bowshockpy import utils as ut
from bowshockpy.cubemass import BowshockCube
from bowshockpy.cubeproc import CubeProcessing
from bowshockpy.modelproj import ObsModel
from bowshockpy.models import BowshockModel
from bowshockpy.version import __version__


def generate_bowshock(p):
    """
    Generates the bowshock model from an input file

    Parameters:
    -----------
    p : dict
        Dictionary with all the parameters from the input file
    """
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
        pss += [
            {
                "modelname": p.modelname,
                "L0": (p.__getattribute__(f"L0_{i+1}") * p.distpc * u.au)
                .to(u.km)
                .value,
                "zj": (p.__getattribute__(f"zj_{i+1}") * p.distpc * u.au)
                .to(u.km)
                .value,
                "vj": p.__getattribute__(f"vj_{i+1}"),
                "va": p.__getattribute__(f"va_{i+1}"),
                "v0": p.__getattribute__(f"v0_{i+1}"),
                "rbf_obs": (
                    (p.__getattribute__(f"rbf_obs_{i+1}") * p.distpc * u.au)
                    .to(u.km)
                    .value
                    if p.__getattribute__(f"rbf_obs_{i+1}") is not None
                    else p.__getattribute__(f"rbf_obs_{i+1}")
                ),
                "mass": p.__getattribute__(f"mass_{i+1}"),
            }
        ]

        psobss += [
            {
                "i_deg": p.__getattribute__(f"i_{i+1}"),
                "pa_deg": p.__getattribute__(f"pa_{i+1}"),
                "vsys": p.vsys,
                "distpc": p.distpc,
            }
        ]

    make_output_cubes = len(p.outcubes) != 0

    if make_output_cubes:
        pscube = {
            "nphis": p.nphis,
            "nc": p.nc,
            "vt": p.vt,
            "vch0": p.vch0,
            "vchf": p.vchf,
            "chanwidth": p.chanwidth,
            "nxs": p.nxs,
            "nys": p.nys,
            "nzs": p.nzs,
            "refpix": p.refpix,
            "xpmax": p.xpmax,
            "parot": p.parot,
            "papv": p.papv,
            "bmaj": p.bmaj,
            "bmin": p.bmin,
            "pabeam": p.pabeam,
            "cic": p.cic,
            "tolfactor_vt": p.tolfactor_vt,
            "sigma_beforeconv": p.sigma_beforeconv,
            "maxcube2noise": p.maxcube2noise,
            "verbose": p.verbose,
        }
        mpars = {
            "J": p.J,
            "nu": p.nu * u.GHz,
            "abund": p.abund,
            "mu": p.mu * u.Debye,
            "meanmolmass": p.meanmolmass,
            "Tex": p.Tex * u.K,
            "Tbg": p.Tbg * u.K,
            "tau_custom_function": (
                p.tau_custom_function if "tau_custom_function" in p.__dict__ else None
            ),
            "Inu_custom_function": (
                p.Inu_custom_function if "Inu_custom_function" in p.__dict__ else None
            ),
            "ra_source_deg": p.ra_source_deg,
            "dec_source_deg": p.dec_source_deg,
            "coordcube": p.coordcube,
        }
    else:
        pscube = {}
        mpars = {}

    bscs = []
    for i, (ps, psobs) in enumerate(zip(pss, psobss)):
        bsm = BowshockModel(
            L0=ps["L0"],
            zj=ps["zj"],
            vj=ps["vj"],
            va=ps["va"],
            v0=ps["v0"],
            mass=ps["mass"],
            distpc=psobs["distpc"],
            rbf_obs=ps["rbf_obs"],
        )
        bsmobs = ObsModel(
            model=bsm,
            i_deg=psobs["i_deg"],
            pa_deg=psobs["pa_deg"],
            vsys=psobs["vsys"],
        )
        if i == 0:
            ut.make_folder(f"models/{ps['modelname']}")
        plt_model = bsm.get_modelplot(
            modelname=ps["modelname"] + f" bowshock_{i+1}",
        )
        plt_model.plot()
        plt_model.savefig(
            f"models/{ps['modelname']}/bowshock_model_{i+1}.pdf",
        )
        plt_obsmodel = bsmobs.get_obsmodelplot(
            modelname=ps["modelname"] + f" bowshock_{i+1}"
        )
        plt_obsmodel.plot()
        plt_obsmodel.savefig(
            f"models/{ps['modelname']}/bowshock_projected_{i+1}.jpg",
            dpi=300,
        )
        if make_output_cubes:
            print(
                f"""

Generating bowshock {i+1}/{p.nbowshocks}
                  """
            )
            bscs += [
                BowshockCube(
                    obsmodel=bsmobs,
                    nphis=pscube["nphis"],
                    xpmax=pscube["xpmax"],
                    vch0=pscube["vch0"],
                    vchf=pscube["vchf"],
                    chanwidth=pscube["chanwidth"],
                    nzs=pscube["nzs"],
                    nc=pscube["nc"],
                    nxs=pscube["nxs"],
                    nys=pscube["nys"],
                    refpix=pscube["refpix"],
                    cic=pscube["cic"],
                    vt=pscube["vt"],
                    tolfactor_vt=pscube["tolfactor_vt"],
                    verbose=pscube["verbose"],
                )
            ]
            bscs[i].makecube()
            print(
                f"""
Central velocity of the first channel: {bscs[i].vch0:.2f} km/s
Central velocity of the last channel: {bscs[i].vchf:.2f} km/s
Channel width: {bscs[i].chanwidth:.3f} km/s
Pixel size: {bscs[i].arcsecpix:.5f} arcsec/pix
     """
            )

    print(
        f"""
The masses have been computed!

The cubes are going to be processed in order to get the desired outputs
specified in {p.filename}. The outputs will be saved in fits format. The
filename of each cube indicate its quantity and the operations applied to the
cube ("<quantity>_<operations>.fits"). Some abbreviations will be used in the
name of the fits files:

Abbreviations for quantities are:             Abbreviations for the operations are:
   m: mass [SolarMass]                             s: add_source
   I: Intensity [Jy/beam]                          r: rotate
   Ntot: Total column density [cm-2]               n: add_noise
   Nmol: Column density emitting molecule [cm-2]   c: convolve
   tau: Opacity
"""
    )
    bscp = CubeProcessing(
        bscs,
        modelname=p.modelname,
        J=mpars["J"],
        nu=mpars["nu"],
        abund=mpars["abund"],
        meanmolmass=mpars["meanmolmass"],
        mu=mpars["mu"],
        Tex=mpars["Tex"],
        Tbg=mpars["Tbg"],
        tau_custom_function=mpars["tau_custom_function"],
        Inu_custom_function=mpars["Inu_custom_function"],
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
    for ck in bscp.listmompvs:
        bscp.plot_channels(
            ck,
            savefig=f"models/{p.modelname}/bowshock_cube_{ck}.pdf",
            add_beam=True,
        )

    bscp.momentsandpv_and_params_all(
        savefits=p.savefits,
        saveplot=p.saveplot,
        mom1clipping=p.mom1clipping,
        mom2clipping=p.mom2clipping,
        mom0values=p.mom0values,
        mom1values=p.mom1values,
        mom2values=p.mom2values,
        mom8values=p.mom8values,
        pvvalues=p.pvvalues,
        add_beam=True,
    )

    # Save the file with all the parameters used to generate the bowshocks
    os.system(f"cp {p.filename.rstrip('.py')}.py models/{p.modelname}")


def main():
    """
    This is called when BowshockPy is run from the terminal
    """

    description = """
BowshockPy is a Python package that generates synthetic spectral cubes,
position-velocity diagrams, and moment images for a simple analytical jet-driven
bowshock model, using the prescription for protostellar jets presented in
Ostriker et al. (2001) and Tabone et al. (2018). Please, see the documentation
at:

https://bowshockpy.readthedocs.io/en/latest/

    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-r",
        "--read",
        dest="parameters_file",
        type=str,
        help="Reads a configuration file to generate the bowshock model",
        default="None",
    )
    parser.add_argument(
        "-p",
        "--print",
        dest="inputfile_example",
        type=str,
        help="""
        Prints an example of input file. Write the number of the corresponding
        example that is closer to your needs. There are 3 examples: write 1 to
        print an example of input file of a redshifted bowshock model, write 2
        for a model including two redshifted bowshocks, write 3 for a
        blueshifted bowshock. See https://bowshockpy.readthedocs.io/en/latest/
        for a detailed documentation of the examples.  """,
        default="None",
    )

    examples_available = [
        "example1.py",
        "example2.py",
        "example3.py",
        "example4.py",
        "example5.py",
    ]
    args = parser.parse_args()
    filename = args.parameters_file
    _example = args.inputfile_example
    if filename != "None":
        parameters = runpy.run_path(filename)
        p = ut.VarsInParamFile(parameters)
        generate_bowshock(p)
    if _example != "None":
        example = _example if _example.endswith(".py") else f"{_example}.py"
        if example in examples_available:
            ut.print_example(example)
            print(f"{example} has been created")
        else:
            print(f"{example} file is not available and could not be created")


if __name__ == "__main__":
    main()
