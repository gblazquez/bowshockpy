"""This module contains the tools that allows to compute the column densities,
opacities, and intensities of the bowshock model spectral cubes"""

import warnings
from datetime import datetime

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.gridspec import GridSpec
from scipy.ndimage import rotate

import bowshockpy.plots as pl
import bowshockpy.radtrans as rt
import bowshockpy.rotlinearmol as rlm
import bowshockpy.utils as ut
from bowshockpy import moments
from bowshockpy.cubemass import MassCube
from bowshockpy.version import __version__

warnings.formatwarning = ut.formatwarning
warnings.filterwarnings("error", category=ut.UserError)
warnings.filterwarnings("always", category=UserWarning)


class CubeProcessing(MassCube):
    """
    Process a MassCube instance

    Parameters
    ----------
    bscube :  class instance
        Instance of MassCube
    modelname : str, optional
        Name of the folder in /models where the outputs will be saved
    J : int, optional
        Upper level of the rotational transition (e.g. 3 for transition
        "3-2")
    nu : float | astropy.units.Quantity, optional
        Frequency of the transition. If float, the units of nu should be GHz
    abund : float, optional
        Abundance relative to the molecular hydrogen
    meanmolmass : float, optional
        Mean mass per H molecule.
    mu : float | astropy.unit.Quantity, optional
        Permanent dipole moment of the molecule. If float, the units of mu
        should be Debye.
    Tex : float | astropy.unit.Quantity, optional
        Excitation temperature. If float, the units of Tex should be Kelvin.
    Tbg : float | astropy.unit.Quantity, optional
        Excitation temperature. If float, the units of Tex should be Kelvin.
    tau_custom_function : callable, optional
        By default, BowshockPy will compute the opacities from a rotational
        transition of a linear molecule (neglecting vibrational excited states
        and centrifugal distortion of the molecule). You can provide a custom
        function to compute the opacities from the column densities per
        velocity bin. This function should accept only the column densities per
        velocity bin, and return the opacities.
    Inu_custom_function : callable, optional
        Custom function to compute the intensities from the opacities
    coordcube : str, optional
        Set to "sky" if you would like to set the cube headers in sky
        coordinates, or "offset" if you prefer them in offsets relative to the
        origin (the source).
    ra_source_deg : float, optional
        Source right ascension [deg]
    dec_source_deg : float, optional
        Source declination [deg]
    bmin : float | None, optional
        Beam minor axis [arcsec]
    bmaj : float | None, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    papv : float
        Position angle used to calculate the PV [degrees]
    parot : float
        Angle to rotate the image [degrees]
    sigma_beforeconv : float
        Standard deviation of the noise of the map, before convolution. Set to
        None if maxcube2noise is used.
    maxcube2noise : float
        Standard deviation of the noise of the map, before convolution, relative
        to the maximum pixel in the cube. The actual noise will be computed
        after convolving. This parameter would not be used if sigma_beforeconve
        is not None.
    verbose : bolean, optional
        Set True to verbose messages about the computation
    kwargs : optional
        Keyword arguments into `~bowshockpy.plot.plot_channel`

    Attributes
    ----------
    x_FWHM : float | None
        Full width half maximum of the Gaussian beam for the x direction
        [pixel]
    y_FWHM : float | None
        Full width half maximum of the Gaussian beam for the y direction
        [pixel]
    beamarea : float | None
        Area of the beam [pixel^2]
    cubes : dict
        Dictionary of the processed cubes. Keys are abbreviations of the
        quantity of the cube and the operations performed to it
    refpixs : dict
        Dictionary of the reference pixel of the cubes.  Keys are abbreviations
        of the quantity of the cube and the operations performed to it
    hdrs : dict
        Dictionary of the headers `astropy.io.fits.header.Header` of each cube.
        The headers are generated when savecube method is used.
    areapix_cm : float
        Area of a pixel in cm.
    beamarea_sr : `astropy.units.Quantity` | None
        Area of the beam in stereoradians. If no beam is provided, beamarea_sr
        will be None, and all the intensities will be expressed in Jy/arcsec^2
    listmompvs : list
        List of cubes to which the moments and the position velocity diagrams
        are going to performed when the method self.momentsandpv_all and
        self.momentsandpv_and_params_all are called
    """

    default_kwargs = {}
    btypes = {
        "m": "mass",
        "I": "Intensity",
        "Ntot": "Total column density",
        "Nmol": "Emitting molecule column density",
        "tau": "Opacity",
    }
    btypes_colorbar = {
        "m": r"mass / (pixel $\times$ channel)",
        "I": "Intensity",
        "Ntot": "Total column density / channel",
        "Nmol": "Emitting molecule col. dens. / channel",
        "tau": "Opacity",
    }
    bunits_default = {
        "m": "solMass",
        "I": "Jy arcsec-2",
        "Ntot": "cm-2",
        "Nmol": "cm-2",
        "tau": "-",
    }
    dos = {
        "s": "add_source",
        "r": "rotate",
        "n": "add_noise",
        "c": "convolve",
    }
    momtol_clipping = 10 ** (-4)
    attribs_to_get_from_cubes = [
        "arcsecpix",
        "nxs",
        "nys",
        "nc",
        "vch0",
        "velchans",
        "vchf",
        "xpmax",
        "distpc",
        "refpix",
        "chanwidth",
        "abschanwidth",
        "vsys",
    ]

    def __init__(
        self,
        modelcubes,
        modelname="none",
        J=3,
        nu=345.79598990 * u.GHz,
        abund=8.5 * 10 ** (-5),
        meanmolmass=2.8,
        mu=0.112 * u.D,
        Tex=100 * u.K,
        Tbg=2.7 * u.K,
        tau_custom_function=None,
        Inu_custom_function=None,
        coordcube="offset",
        ra_source_deg=None,
        dec_source_deg=None,
        bmin=None,
        bmaj=None,
        pabeam=0,
        papv=0.0,
        parot=0.0,
        sigma_beforeconv=None,
        maxcube2noise=None,
        verbose=True,
        **kwargs,
    ):

        if isinstance(modelcubes, list):
            self.nmodels = len(modelcubes)
            self.combine_cubes(modelcubes)
        else:
            self.nmodels = 1
            modelcubes = [modelcubes]

        for att in self.attribs_to_get_from_cubes:
            setattr(self, att, modelcubes[0].__getattribute__(att))

        self.ies = np.array([mc.i for mc in modelcubes])
        self.L0s = np.array([mc.L0_arcsec for mc in modelcubes])
        self.zjs = np.array([mc.zj_arcsec for mc in modelcubes])
        self.vjs = np.array([mc.vj for mc in modelcubes])
        self.vas = np.array([mc.va for mc in modelcubes])
        self.v0s = np.array([mc.v0 for mc in modelcubes])
        self.rbfs = np.array([mc.rbf_arcsec for mc in modelcubes])
        self.tjs = np.array([mc.tj_yr for mc in modelcubes])
        self.masss = np.array([mc.mass for mc in modelcubes])
        self.rhoas = np.array([mc.rhoa_gcm3 for mc in modelcubes])
        self.m0s = np.array([mc.mp0_solmassyr for mc in modelcubes])
        self.mwfs = np.array([mc.mpamb_f_solmassyr for mc in modelcubes])

        self.modelname = modelname
        self.J = J
        self.nu = ut.make_astropy_units(nu, u.GHz)
        self.rottrans = f"{int(J)}-{int(J-1)}"
        self.abund = abund
        self.meanmolmass = meanmolmass
        self.mu = ut.make_astropy_units(mu, u.Debye)
        self.Tex = ut.make_astropy_units(Tex, u.K)
        self.Tbg = ut.make_astropy_units(Tbg, u.K)
        self.tau_custom_function = tau_custom_function
        self.Inu_custom_function = Inu_custom_function
        self.coordcube = coordcube
        self.ra_source_deg = ra_source_deg
        self.dec_source_deg = dec_source_deg
        self.bmin = bmin
        self.bmaj = bmaj
        self.pabeam = pabeam
        self.papv = papv
        self.parot = parot
        self.sigma_beforeconv = sigma_beforeconv
        self.maxcube2noise = maxcube2noise
        self.verbose = verbose

        for kwarg in self.default_kwargs:
            kwarg_attr = (
                kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            )
            setattr(self, kwarg, kwarg_attr)

        self.cubes = {}
        self.cubes["m"] = self.cube
        self.sigma_noises = {}
        self.sigma_noises["m"] = 0
        self.noisychans = {}
        self.noisychans["m"] = np.zeros_like(self.cube[0])
        self.refpixs = {}
        self.refpixs["m"] = self.refpix
        self.hdrs = {}
        self.listmompvs = []

        self.areapix_cm = None
        self.x_FWHM = None
        self.y_FWHM = None
        self.beamarea = None
        self.beamarea_sr = None
        self.bunits = self.bunits_default

        self._calc_beamarea_sr()
        self._calc_areapix_cm()

    @staticmethod
    def _newck(ck, s):
        return f"{ck}_{s}" if "_" not in ck else ck + s

    @staticmethod
    def _q(ck):
        return ck.split("_")[0] if "_" in ck else ck

    def _getunitlabel(self, ck):
        if self.bunits[self._q(ck)] == "-":
            unitlabel = f"{self.btypes_colorbar[self._q(ck)]}"
        else:
            unitlabel = f"{self.btypes_colorbar[self._q(ck)]} [{self.bunits[self._q(ck)]}]"
        return unitlabel

    def _calc_beamarea_sr(self):
        if self.bmin is not None and self.bmaj is not None:
            self.x_FWHM = self.bmin / self.arcsecpix
            self.y_FWHM = self.bmaj / self.arcsecpix
            self.beamarea = np.pi * self.y_FWHM * self.x_FWHM / (4 * np.log(2))
            self.beamarea_sr = ut.mb_sa_gaussian_f(
                self.bmaj * u.arcsec, self.bmin * u.arcsec
        )
            self.bunits["I"] = "Jy/beam"
        else:
            self.x_FWHM = None
            self.y_FWHM = None
            self.beamarea = None
            self.beamarea_sr = None
            self.bunits["I"] = "Jy arcsec-2"

    def _calc_areapix_cm(self):
        self.areapix_cm = ((self.arcsecpix * self.distpc * u.au) ** 2).to(u.cm**2)

    def _check_combine_possibility(self, modelcubes):
        for att in self.attribs_to_get_from_cubes:
            if not ut.allequal([mc.__getattribute__(att) for mc in modelcubes]):
                raise ValueError(f"Trying to combine cubes with different {att}")

    def combine_cubes(self, modelcubes):
        """
        Combines (sums) a list of cubes

        Parameters
        ----------
        modelcubes : list
            List of cubes to combine
        """
        self._check_combine_possibility(modelcubes)
        self.cube = np.sum([modelcube.cube for modelcube in modelcubes], axis=0)

    def tau_f(self, dNmoldv):
        if self.tau_custom_function is not None:
            return self.tau_custom_function(dNmoldv=dNmoldv)
        return rlm.tau_linearmol(
            dNmoldv=dNmoldv,
            J=self.J,
            nu=self.nu,
            Tex=self.Tex,
            mu=self.mu,
        )

    def I_f(self, tau):
        if self.Inu_custom_function is not None:
            return self.Inu_custom_function(tau=tau)
        return rt.Inu_func(
            tau=tau,
            nu=self.nu,
            Tex=self.Tex,
            Tbg=self.Tbg,
        )

    def calc_Ntot(
        self,
    ):
        """
        Computes the total (molecular hydrogen + heavier components) column
        densities of the model cube
        """
        if self.verbose:
            print("\nComputing column densities...")

        self.cubes["Ntot"] = (
            rt.column_density_tot(
                m=self.cubes["m"] * u.solMass,
                area=self.areapix_cm,
                meanmolmass=self.meanmolmass,
            )
            .to(u.cm ** (-2))
            .value
        )

        self.refpixs["Ntot"] = self.refpixs["m"]
        self.noisychans["Ntot"] = self.noisychans["m"]
        self.sigma_noises["Ntot"] = self.sigma_noises["m"]
        if self.verbose:
            print("column densities has been calculated (Ntot cube)\n")

    def calc_Nmol(self):
        """
        Computes the emitting molecule column densities of the model cube
        """
        if "Ntot" not in self.cubes:
            self.calc_Ntot()
        if self.verbose:
            print("\nComputing column densities of the emitting molecule...")

        self.cubes["Nmol"] = (
            rt.column_density_mol(
                Ntot=self.cubes["Ntot"] * u.cm ** (-2),
                abund=self.abund,
            )
            .to(u.cm ** (-2))
            .value
        )

        self.refpixs["Nmol"] = self.refpixs["m"]
        self.noisychans["Nmol"] = self.noisychans["m"]
        self.sigma_noises["Nmol"] = self.sigma_noises["m"]
        if self.verbose:
            print(
                "The column densities of the emitting molecule have been calculated (Nmol cube)\n"
            )

    def calc_tau(self):
        """
        Computes the opacities of the model cube
        """
        if "Nmol" not in self.cubes:
            self.calc_Nmol()
        if self.verbose:
            print("\nComputing opacities...")

        self.cubes["tau"] = (
            self.tau_f(
                dNmoldv=self.cubes["Nmol"]
                * u.cm ** (-2)
                / (self.abschanwidth * u.km / u.s)
            )
            .to("")
            .value
        )

        self.refpixs["tau"] = self.refpixs["m"]
        self.noisychans["tau"] = self.noisychans["m"]
        self.sigma_noises["tau"] = self.sigma_noises["m"]
        if self.verbose:
            print("Opacities has been calculated (tau cube)\n")

    def calc_I(self):
        """
        Calculates the intensity [Jy/beam] of the model cube.
        """
        if "tau" not in self.cubes:
            self.calc_tau()
        if self.verbose:
            print("\nComputing intensities...")

        if self.beamarea_sr is not None:
            self.cubes["I"] = (
                (self.I_f(tau=self.cubes["tau"]) * self.beamarea_sr).to(u.Jy).value
            )
        else:
            self.cubes["I"] = (
                (self.I_f(tau=self.cubes["tau"])).to(u.Jy / u.arcsec**2).value
            )

        self.refpixs["I"] = self.refpixs["m"]
        self.noisychans["I"] = self.noisychans["m"]
        self.sigma_noises["I"] = self.sigma_noises["m"]
        if self.verbose:
            print("Intensities has been calculated (I cube)\n")

    def add_source(self, ck="m", value=None):
        """
        Adds a source to the cube in the reference pixel

        Parameters
        -----------
        ck : str, optional
            Key of the cube to add the source
        value : float, optional
            Pixel value of the source. If None, the maximum of the cube will be
            considered
        """
        nck = self._newck(ck, "s")
        if self.verbose:
            print(f"\nAdding source to {ck}...")
        self.cubes[nck] = np.copy(self.cubes[ck])
        value = value if value is not None else np.max(self.cubes[ck])
        if self.refpixs[ck][1] >= 0 and self.refpixs[ck][0] >= 0:
            self.cubes[nck][:, self.refpix[1], self.refpix[0]] = value
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = self.noisychans[ck]
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"""
{nck} cube has been created by adding a point source to {ck}, in pix
[{self.refpixs[nck][0]:.2f}, {self.refpixs[nck][1]:.2f}] pix\n
"""
            )

    def rotate(self, ck="m", forpv=False):
        """
        Rotates the cube an angle self.parot

        Parameters
        -----------
        ck : str, optional
            Key of the cube to rotate
        forpv : bool, optional
            If True, the image is rotated to calculate the PV along the
            bowshock axis
        """

        nck = self._newck(ck, "r") if not forpv else self._newck(ck, "R")
        if self.verbose:
            if forpv:
                print(f"\nRotatng {ck} in order to compute the PV diagram...")
            else:
                print(f"\nRotating {ck}...")
        # before allowing rotation of the model and not the cube
        # angle = -self.pa-90 if not forpv else self.pa+90
        # after allowing the model to be rotated
        angle = -self.parot if not forpv else self.papv + 90
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            self.cubes[nck][chan] = rotate(
                self.cubes[ck][chan], angle=angle, reshape=False, order=1
            )
        ang = angle * np.pi / 180
        centerx = (self.nxs - 1) / 2
        centery = (self.nys - 1) / 2
        rp_center_x = self.refpixs[ck][0] - centerx
        rp_center_y = self.refpixs[ck][1] - centery
        self.refpixs[nck] = [
            +rp_center_x * np.cos(ang) + rp_center_y * np.sin(ang) + centerx,
            -rp_center_x * np.sin(ang) + rp_center_y * np.cos(ang) + centery,
        ]
        self.noisychans[nck] = rotate(
            self.noisychans[ck],
            angle=angle,
            reshape=False,
            order=1,
        )
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"""
{nck} cube has been created by rotating {ck} cube an angle {angle} deg to
compute the PV diagram
""")

    def add_noise(self, ck="m"):
        """
        Adds Gaussian noise to the cube.

        Parameters
        -----------
        ck : str, optional
            Key of the cube to rotate
        """
        nck = self._newck(ck, "n")
        if self.verbose:
            print(f"\nAdding noise to {ck}...")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            # sigma_noise = self.target_noise * 2 * np.sqrt(np.pi) \
            #          * np.sqrt(self.x_FWHM*self.y_FWHM) / 2.35
            sigma_noise = (
                self.sigma_beforeconv
                if self.sigma_beforeconv is not None
                else np.max(self.cubes[ck]) * self.maxcube2noise
            )
            noise_matrix = np.random.normal(
                0, sigma_noise, size=np.shape(self.cubes[ck][chan])
            )
            self.cubes[nck][chan] = self.cubes[ck][chan] + noise_matrix
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = noise_matrix
        self.sigma_noises[nck] = sigma_noise
        if self.verbose:
            print(f"""
{nck} cube has been created by adding Gaussian noise to {ck} cube
"""
                  )

    def convolve(self, ck="m"):
        """
        Convolves the cube with the defined Gaussian kernel (self.bmaj,
        self.bmin, self.pabeam)

        Parameters
        -----------
        ck : str, optional
            Key of the cube to convolve
        """

        if self.beamarea is None:
            warnings.warn(
                message="""
The major axis and minor axis of the beam to be used for the convolution (bmaj
and bmin parameters, respectively) were not provided . Please, instantiate the
class CubeProcessing again providing values for both beam axis.
""",
                category=ut.UserError,
            )
        else:
            nck = self._newck(ck, "c")
            if self.verbose:
                print(f"\nConvolving {ck}... ")
            self.cubes[nck] = np.zeros_like(self.cubes[ck])

            if self.verbose:
                ts = []
                ut.progressbar_bowshock(
                    0, self.nc, length=50, timelapsed=0, intervaltime=0
                )
            for chan in range(np.shape(self.cubes[ck])[0]):
                if self.verbose:
                    t0 = datetime.now()
                self.cubes[nck][chan] = ut.gaussconvolve(
                    self.cubes[ck][chan],
                    x_FWHM=self.x_FWHM,
                    y_FWHM=self.y_FWHM,
                    pa=self.pabeam,
                    return_kernel=False,
                )
                if self.verbose:
                    tf = datetime.now()
                    intervaltime = (tf - t0).total_seconds()
                    ts.append(intervaltime)
                    ut.progressbar_bowshock(
                        chan + 1, self.nc, np.sum(ts), intervaltime, length=50
                    )
            self.refpixs[nck] = self.refpixs[ck]
            self.noisychans[nck] = ut.gaussconvolve(
                self.noisychans[ck],
                x_FWHM=self.x_FWHM,
                y_FWHM=self.y_FWHM,
                pa=self.pabeam,
                return_kernel=False,
            )
            self.sigma_noises[nck] = np.std(self.noisychans[nck])
            if self.verbose:
                print(
                    f"""
{nck} cube has been created by convolving {ck} cube with a Gaussian kernel of
size [{self.x_FWHM:.2f}, {self.y_FWHM:.2f}] pix and PA of {self.pabeam:.2f}deg
"""
                )
                if "n" in nck:
                    print(
                        f"""
The rms of the convolved image is {self.sigma_noises[nck]:.5} {self.bunits[self._q(nck)]}
"""
                    )

    def _useroutputcube2dostr(self, userdic):
        # TODO: Move to utils or a function above the class?
        dictrad = {
            "mass": "m",
            "intensity": "I",
            "emitting_molecule_column_density": "Nmol",
            "total_column_density": "Ntot",
            "opacity": "tau",
            "add_source": "s",
            "rotate": "r",
            "add_noise": "n",
            "convolve": "c",
        }
        dostrs = []
        for userkey in userdic:
            q = dictrad[userkey]
            ops = userdic[userkey]
            calcmompv = "moments_and_pv" in ops
            if calcmompv:
                if len(ops) > 1:
                    ss = "".join(
                        [
                            dictrad[s_user]
                            for s_user in userdic[userkey]
                            if s_user != "moments_and_pv"
                        ]
                    )
                    dostr = [f"{q}_{ss}"]
                elif len(ops) == 1:
                    dostr = [f"{q}"]
                else:
                    dostr = []
                dostrs += dostr
                self.listmompvs += dostr
            else:
                if len(ops) != 0:
                    ss = "".join([dictrad[s_user] for s_user in userdic[userkey]])
                    dostrs += [f"{q}_{ss}"]
                else:
                    dostrs += [f"{q}"]
        return dostrs

    def calc(self, userdic):
        """
        Computes the quantities and the operations to the cubes.

        Parameters
        -----------
        userdic : dict
            Dictionary indicating the desired output spectral cubes and the
            operations performed over them. The keys of the dictionary are
            strings indicating the quantities of the desired cubes. These are
            the available quantities of the spectral cubes:

            - "mass": Total mass of molecular hydrogen in solar mass.
            - "mol_column_density": Column density of the emitting molecule in
              cm-2.
            - "intensity": Intensity in Jy/beam.
            - "tau": Opacities.

            The values of the dictionary are lists of strings indicating the
            operations to be performed over the cube. These are the available
            operations:

            - "add_source": Add a source at the reference pixel, just for
              spatial reference purposes.
            - "rotate": Rotate the whole spectral cube by an angle given by
              parot parameter.
            - "add_noise": Add Gaussian noise, defined by maxcube2noise
              parameter.
            - "convolve": Convolve with a Gaussian defined by the parameters
              bmaj, bmin, and pabeam.
            - "moments_and_pv": Computes the moments 0, 1, and 2, the maximum
              intensity and the PV diagram.

            The operations will be performed folowing the order of the strings
            in the list (from left to right). The list can be left empty if no
            operations are desired.

        Example:
        --------
        >>> cp = CubeProcessing(...)
        >>> outcubes = {
        >>>    "intensity": ["add_noise", "convolve", "moments_and_pv"],
        >>>    "opacity": [],
        >>>    "mol_column_density": ["convolve"],
        >>>    "mass": [],
        >>> }
        >>> cp.calc(outcubes)

        will save 4 spectral cubes in fits format. The first one are the
        intensities with Gaussian noise added, it will be convolved, and the
        moments and PV diagrams will be computed; the second cube will be the
        opacity; the third will be the mol_column_density, which will be
        convolved; and the forth cube will be the masses. The first spectral
        cube will be named I_nc.fits, the second tau.fits, the third
        Nmol_c.fits, and the fourth m.fits.
        """
        dostrs = self._useroutputcube2dostr(userdic)
        for ds in dostrs:
            _split = ds.split("_")
            q = _split[0]
            if q not in self.cubes:
                self.__getattribute__(f"calc_{q}")()
            if len(_split) > 1:
                ss = _split[1]
                for i, s in enumerate(ss):
                    ck = q if i == 0 else f"{q}_{ss[:i]}"
                    if self._newck(ck, s) not in self.cubes:
                        self.__getattribute__(self.dos[s])(ck=ck)

    def savecube(self, ck, fitsname=None):
        """
        Saves the cube in fits format

        Parameters
        -----------
        ck : str
            Key of the cube to convolve
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}.fits. If the path does not
            exist, it will be created.
        """
        hdr = ut.get_default_hdr(naxis=3)
        if self.coordcube == "offset":
            ctype1 = "OFFSET"
            ctype2 = "OFFSET"
            cunit1 = "arcsec"
            cunit2 = "arcsec"
            crval1 = 0.0
            crval2 = 0.0
            cdelt1 = self.arcsecpix
            cdelt2 = self.arcsecpix
        else:
            ctype1 = "RA---SIN"
            ctype2 = "DEC--SIN"
            cunit1 = "deg"
            cunit2 = "deg"
            crval1 = self.ra_source_deg
            crval2 = self.dec_source_deg
            cdelt1 = -self.arcsecpix / 3600
            cdelt2 = self.arcsecpix / 3600
        hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
        hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
        hdr["NAXIS3"] = np.shape(self.cubes[ck])[0]
        if self.beamarea is not None:
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
        hdr["BTYPE"] = self.btypes[self._q(ck)]
        hdr["OBJECT"] = f"{self.modelname}"
        hdr["BUNIT"] = self.bunits[self._q(ck)]
        hdr["CTYPE1"] = ctype1
        hdr["CRVAL1"] = crval1
        hdr["CDELT1"] = cdelt1
        hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
        hdr["CUNIT1"] = cunit1
        hdr["CTYPE2"] = ctype2
        hdr["CRVAL2"] = crval2
        hdr["CDELT2"] = cdelt2
        hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
        hdr["CUNIT2"] = cunit2
        hdr["CTYPE3"] = "VRAD"
        hdr["CRVAL3"] = self.velchans[0]
        hdr["CDELT3"] = self.velchans[1] - self.velchans[0]
        hdr["CRPIX3"] = 1.0
        hdr["CUNIT3"] = "km/s"
        hdr["RESTFRQ"] = self.nu.to(u.Hz).value

        self.hdrs[ck] = hdr
        hdu = fits.PrimaryHDU(self.cubes[ck])
        hdul = fits.HDUList([hdu])
        hdu.header = self.hdrs[ck]
        if fitsname is None:
            savefolder = f"models/{self.modelname}/fits"
            ut.make_folder(foldername=savefolder)
            fitsname = f"models/{self.modelname}/fits/{ck}.fits"

        hdul.writeto(fitsname, overwrite=True)
        if self.verbose:
            print(f"{fitsname} saved")

    def savecubes(self, userdic):
        """
        Saves the cubes specified by userdic

        Parameters
        -----------
        userdic : dict
            Dictionary indicating the desired output spectral cubes and the
            operations performed over them.

        Example:
        --------
        >>> cp = CubeProcessing(...)
        >>> outcubes = {
        >>>    "intensity": ["add_noise", "convolve", "moments_and_pv"],
        >>>    "opacity": [],
        >>>    "mol_column_density": ["convolve"],
        >>>    "mass": [],
        >>> }
        >>> cp.savecubes(outputcubes)

        will save 4 spectral cubes in fits format. The first one are the
        intensities with Gaussian noise added, it will be convolved, and the
        moments and PV diagrams will be computed; the second cube will be the
        opacity; the third will be the mol_column_density, which will be
        convolved; and the forth cube will be the masses. The first spectral
        cube will be named I_nc.fits, the second tau.fits, the third
        Nmol_c.fits, and the fourth m.fits.
        """
        cks = self._useroutputcube2dostr(userdic)
        for ck in cks:
            self.savecube(ck)

    def plot_channel(
        self,
        ck,
        chan,
        vmax=None,
        vmin=None,
        cmap="inferno",
        savefig=None,
        add_beam=False,
        return_fig_axs=False,
    ):
        """
        Plots a channel map of a cube

        Parameters
        ----------
        ck : str
            Key of the cube to plot (see keys of self.cubes dictionary)
        chan : int
            Channel map to plot
        vmax : float, optional
            Maximum value of the colormap. If None (default), the maximum value
            of the channel is chosen.
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value
            of the channel is chosen.
        cmap : str, optional
            Label of the colormap, by default "inferno".
        savefig : str, optional String of the full path to save the figure. If
            None, no figure is saved. By default, None.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        return_fig_axs : bool, optional
            If True, returns a tuple of the ax of the channel map and the
            colorbar.  If False, does not return anything.

        Returns
        --------
        (fig, ax, cbax) : tuple of matplotlib.axes.Axes Axes of the channel map
            and the colorbar, only returns if return_fig_axs=True.
        """
        add_beam = add_beam if "c" in ck else False
        fig, axs, cbax = pl.plot_channel(
            cube=self.cubes[ck],
            chan=chan,
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            vmax=vmax,
            vmin=vmin,
            cmap=cmap,
            units=self._getunitlabel(ck),
            refpix=self.refpix,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            return_fig_axs=True,
        )
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")
        if return_fig_axs:
            return fig, axs, cbax

    def plot_channels(
        self, ck, savefig=None, add_beam=False, return_fig_axs=False, **kwargs
    ):
        """
        Plots several channel map of a spectral cube.

        Parameters
        ----------
        ck : str
            Key of the cube to plot (see keys of self.cubes dictionary)
        ncol : int, optional
            Number of columns in the figure, by default 4
        nrow : int, optional
            Number of rows of the figure, by default 4
        figsize : tuple, optional
            Size of the figure. If None, an optimal size will be computed. By
            default None.
        wspace : float, optional
            Width space between the channel plots, by default 0.05
        hspace : float, optional
            Height space between the cannel plots, by default 0.0
        vmax : float, optional
            Maximum value of the colormap. If None (default), the maximum value
            of the channel is chosen.
        vcenter : _type_, optional
            _description_, by default None
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value
            of the channel is chosen.
        cmap : str, optional
            Label of the colormap, by default "inferno".
        units : str, optional
            Units of the values of the cube, by default "Mass [Msun]"
        xmajor_locator : float, optional
            Major locator in x-axis, by default 1
        xminor_locator : float, optional
            Minor locator in x-axis, by default 0.2
        ymajor_locator : float, optional
            Major locator in y-axis, by default 1
        yminor_locator : float, optional
            Minor locator in y-axis, by default 0.2
        refpix : list, optional
            Pixel of reference, by default [0,0]
        savefig : str, optional String of the full path to save the figure. If
            None, no figure is saved. By default, None.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        return_fig_axs : bool, optional
            If True, returns the figure, axes of the channel map, and the axes the
            colorbar.  If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs=True.
        """
        add_beam = add_beam if "c" in ck else False
        fig, axs, cbax = pl.plot_channels(
            cube=self.cubes[ck],
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            units=self._getunitlabel(ck),
            refpix=self.refpix,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            return_fig_axs=True,
            **kwargs,
        )
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")
        if return_fig_axs:
            return fig, axs, cbax

    def pvalongz(self, ck, halfwidth=0, savefits=False, fitsname=None):
        """
        Performs the position velocity diagram along the self.papv direction

        Parameters
        -----------
        ck : str
            Key of the cube to perform the PV-diagram.
        halfwidth : int, optional
            Number of pixels around xpv that will be taking into account to
            compute the PV-diagram.
        savefits : boolean
            If True, save the PV-diagram in fits format.
        fitsname : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        pvimage : numpy.ndarray
            Position velocity diagram
        """
        pvimage = moments.pv(
            self.cubes[ck], int(self.refpixs[ck][1]), halfwidth=halfwidth, axis=1
        )
        if savefits:
            hdrpv = ut.get_default_hdr(naxis=2, beam=False, pv=True)
            hdrpv["NAXIS"] = 2
            hdrpv["NAXIS1"] = np.shape(self.cubes[ck])[1]
            hdrpv["NAXIS2"] = np.shape(self.cubes[ck])[0]
            hdrpv["BTYPE"] = self.btypes[self._q(ck)]
            hdrpv["OBJECT"] = f"{self.modelname}"
            hdrpv["BUNIT"] = self.bunits[self._q(ck)]
            hdrpv["CTYPE1"] = "OFFSET"
            hdrpv["CRVAL1"] = 0.0
            hdrpv["CDELT1"] = self.arcsecpix
            hdrpv["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdrpv["CUNIT1"] = "arcsec"
            hdrpv["CTYPE2"] = "VRAD"
            hdrpv["CRVAL2"] = self.velchans[0]
            hdrpv["CDELT2"] = self.chanwidth
            hdrpv["CRPIX2"] = 1.0
            hdrpv["CUNIT2"] = "km/s"
            hdrpv["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(pvimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdrpv
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_pv.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_pv.fits saved")
        return pvimage

    def sumint(self, ck, chan_range=None, savefits=False, fitsname=None):
        """
        Computes the image of the summation of pixels of the cube along the
        velocity axis

        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        savefits : boolean
            If True, save the PV-diagram in fits format.
        fitsname : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_sumint.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        sumint : numpy.ndarray
            Image of the summation of the pixels of the cube along the velocty
            axis

        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        sumintimage = moments.sumint(self.cubes[ck], chan_range=chan_range)
        if savefits:
            if self.coordcube == "offset":
                ctype1 = "OFFSET"
                ctype2 = "OFFSET"
                cunit1 = "arcsec"
                cunit2 = "arcsec"
                crval1 = 0.0
                crval2 = 0.0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = "RA---SIN"
                ctype2 = "DEC--SIN"
                cunit1 = "deg"
                cunit2 = "deg"
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600

            hdr = ut.get_default_hdr(naxis=2)
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            if self.beamarea is not None:
                hdr["BMAJ"] = self.bmaj / 3600
                hdr["BMIN"] = self.bmin / 3600
                hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f"{self.modelname}"
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
            hdr["CUNIT2"] = cunit2
            # hdr["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(sumintimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_sumint.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_sumint.fits saved")
        return sumintimage

    def mom0(self, ck, chan_range=None, savefits=False, fitsname=None):
        """
        Computes the 0th order moment along the velocity axis

        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        savefits : boolean
            If True, save the PV-diagram in fits format.
        fitsname : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom0.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom0 : numpy.ndarray
            Moment 0 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0] : chan_range[-1]]
        mom0 = moments.mom0(
            self.cubes[ck],
            chan_vels=chan_vels,
            chan_range=chan_range,
        )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = "OFFSET"
                ctype2 = "OFFSET"
                cunit1 = "arcsec"
                cunit2 = "arcsec"
                crval1 = 0.0
                crval2 = 0.0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = "RA---SIN"
                ctype2 = "DEC--SIN"
                cunit1 = "deg"
                cunit2 = "deg"
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = ut.get_default_hdr(naxis=2)
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            if self.beamarea is not None:
                hdr["BMAJ"] = self.bmaj / 3600
                hdr["BMIN"] = self.bmin / 3600
                hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f"{self.modelname}"
            hdr["BUNIT"] = "Jy/beam.km/s"
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
            hdr["CUNIT2"] = cunit2
            # hdr["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(mom0)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_mom0.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_mom0.fits saved")
        return mom0

    def mom1(self, ck, chan_range=None, clipping=0, savefits=False, fitsname=None):
        """
        Computes the 1th order moment along the velocity axis

        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter
            will be masked with 0 values.
        savefits : boolean, optional
            If True, save the PV-diagram in fits format.
        fitsname : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom1.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom1 : numpy.ndarray
            Moment 1 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0] : chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = (
            clipping if clipping != 0 else self.momtol_clipping * np.max(self.cubes[ck])
        )
        cube_clipped[cube_clipped < clipping] = 0
        mom1 = np.nan_to_num(
            moments.mom1(cube_clipped, chan_vels=chan_vels, chan_range=chan_range)
        )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = "OFFSET"
                ctype2 = "OFFSET"
                cunit1 = "arcsec"
                cunit2 = "arcsec"
                crval1 = 0.0
                crval2 = 0.0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = "RA---SIN"
                ctype2 = "DEC--SIN"
                cunit1 = "deg"
                cunit2 = "deg"
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = ut.get_default_hdr(naxis=2)
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            if self.beamarea is not None:
                hdr["BMAJ"] = self.bmaj / 3600
                hdr["BMIN"] = self.bmin / 3600
                hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f"{self.modelname}"
            hdr["BUNIT"] = "km/s"
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
            hdr["CUNIT2"] = cunit2
            # hdr["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(mom1)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_mom1.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_mom1.fits saved")
        return mom1

    def mom2(self, ck, chan_range=None, clipping=0, savefits=False, fitsname=None):
        """
        Computes the 2th order moment along the velocity axis

        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter
            will be masked with 0 values.
        savefits : boolean, optional
            If True, save the PV-diagram in fits format.
        fitsname : str, optional
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom2.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom2 : numpy.ndarray
            Moment 2 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0] : chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = (
            clipping if clipping != 0 else self.momtol_clipping * np.max(self.cubes[ck])
        )
        cube_clipped[cube_clipped < clipping] = 0
        mom2 = np.nan_to_num(
            moments.mom2(cube_clipped, chan_vels=chan_vels, chan_range=chan_range)
        )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = "OFFSET"
                ctype2 = "OFFSET"
                cunit1 = "arcsec"
                cunit2 = "arcsec"
                crval1 = 0.0
                crval2 = 0.0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = "RA---SIN"
                ctype2 = "DEC--SIN"
                cunit1 = "deg"
                cunit2 = "deg"
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = ut.get_default_hdr(naxis=2)
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            if self.beamarea is not None:
                hdr["BMAJ"] = self.bmaj / 3600
                hdr["BMIN"] = self.bmin / 3600
                hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f"{self.modelname}"
            hdr["BUNIT"] = "km/s"
            hdr["RADESYS"] = "ICRS"
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
            hdr["CUNIT2"] = cunit2
            # hdr["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(mom2)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_mom1.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_mom2.fits saved")
        return mom2

    def maxintens(self, ck, chan_range=None, clipping=0, savefits=False, fitsname=None):
        """
        Computes the maximum value of the cube along the velocity axis

        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the moment.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter
            will be masked with 0 values.
        savefits : boolean, optional
            If True, save the moment in fits format.
        fitsname : str, optional
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_maxintens.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        maxintens : numpy.ndarray
            Maximum value of the pixels of the cubes along the velocity axis

        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        # chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = (
            clipping if clipping != 0 else self.momtol_clipping * np.max(self.cubes[ck])
        )
        cube_clipped[cube_clipped < clipping] = 0
        maxintens = np.nan_to_num(
            moments.maxintens(
                cube_clipped,
                chan_range=chan_range,
            )
        )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = "OFFSET"
                ctype2 = "OFFSET"
                cunit1 = "arcsec"
                cunit2 = "arcsec"
                crval1 = 0.0
                crval2 = 0.0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = "RA---SIN"
                ctype2 = "DEC--SIN"
                cunit1 = "deg"
                cunit2 = "deg"
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = ut.get_default_hdr(naxis=2)
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            if self.beamarea is not None:
                hdr["BMAJ"] = self.bmaj / 3600
                hdr["BMIN"] = self.bmin / 3600
                hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f"{self.modelname}"
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.0
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1.0
            hdr["CUNIT2"] = cunit2
            # hdr["RESTFRQ"] = self.nu.to(u.Hz).value

            hdu = fits.PrimaryHDU(maxintens)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if fitsname is None:
                ut.make_folder(foldername=f"models/{self.modelname}/fits")
                fitsname = f"models/{self.modelname}/fits/{ck}_mom1.fits"
            hdul.writeto(fitsname, overwrite=True)
            if self.verbose:
                print(f"models/{self.modelname}/fits/{ck}_maxintens.fits saved")
        return maxintens

    def plotpv(
        self,
        ck,
        halfwidth,
        ax=None,
        cbax=None,
        savefits=False,
        fitsname=None,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the position velocity diagram.

        Parameters
        -----------
        ck : str
            Key of the cube to which the PV diagram will be computed.
        halfwidth : int, optional
            Number of pixels around xpv that will be taking into account to
            compute the PV-diagram.
        ax : matplotlib.axes.Axes, optional The matplotlib.axes.Axes` instance
            in which the position velodity diagram is drawn. If None, it will
            create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        savefits : bool
            If True, the position velocity diagram will be saved in fits format
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        pvimage : numpy.ndarray
            Position velocity diagram
        """
        ckpv = ck + "R"
        if ckpv not in self.cubes:
            self.rotate(ck, forpv=True)

        pvimage = self.pvalongz(
            ckpv,
            halfwidth=halfwidth,
            savefits=savefits,
            fitsname=fitsname,
        )
        rangex = (
            np.array(
                [-0.5 - self.refpixs[ckpv][0], self.nxs - 0.5 - self.refpixs[ckpv][0]]
            )
            * self.arcsecpix
        )
        fig, axs, cbax = pl.plotpv(
            pvimage,
            rangex=rangex,
            chan_vels=self.velchans,
            ax=ax,
            cbax=cbax,
            cbarlabel=self._getunitlabel(ckpv),
            return_fig_axs=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, pvimage
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plotsumint(
        self,
        ck,
        chan_range=None,
        ax=None,
        cbax=None,
        add_beam=False,
        savefits=False,
        fitsname=False,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the sum of the pixels of the cubes along the velocity axis.

        Parameters
        -----------
        ck : str
            Key of the cube to which the moment will be computed.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        ax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes` instance in which the position velodity
            diagram is drawn. If None, it will create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        savefits : bool
            If True, the moment will be saved in fits format
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_sumint.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        sumint : numpy.ndarray
            Sum of all the pixels along the velocity axis.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        sumint = self.sumint(
            ck,
            chan_range=chan_range,
            savefits=savefits,
            fitsname=fitsname,
        )
        extent = (
            np.array(
                [
                    -(-0.5 - self.refpixs[ck][0]),
                    -(self.nxs - 0.5 - self.refpixs[ck][0]),
                    (-0.5 - self.refpixs[ck][1]),
                    (self.nys - 0.5 - self.refpixs[ck][1]),
                ]
            )
            * self.arcsecpix
        )
        fig, axs, cbax = pl.plotsumint(
            sumint,
            extent=extent,
            interpolation=None,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Integrated " + self._getunitlabel(ck).rstrip("]") + " km/s]",
            return_fig_axs=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, sumint
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plotmom0(
        self,
        ck,
        chan_range=None,
        ax=None,
        cbax=None,
        add_beam=False,
        savefits=False,
        fitsname=None,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the moment 0 (integrated intensity).

        Parameters
        -----------
        ck : str
            Key of the cube to which the moment will be computed.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        ax : matplotlib.axes.Axes, optional
           The matplotlib.axes.Axes` instance in which the position velodity
           diagram is drawn. If None, it will create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom
            corner.
        savefits : bool
            If True, the moment will be saved in fits format
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom0.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        mom0 : numpy.ndarray
            Integrated intensity
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        mom0 = self.mom0(
            ck,
            chan_range=chan_range,
            savefits=savefits,
            fitsname=fitsname,
        )
        extent = (
            np.array(
                [
                    -(-0.5 - self.refpixs[ck][0]),
                    -(self.nxs - 0.5 - self.refpixs[ck][0]),
                    (-0.5 - self.refpixs[ck][1]),
                    (self.nys - 0.5 - self.refpixs[ck][1]),
                ]
            )
            * self.arcsecpix
        )
        fig, axs, cbax = pl.plotmom0(
            mom0,
            extent=extent,
            interpolation=None,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Integrated " + self._getunitlabel(ck).rstrip("]") + " km/s]",
            return_fig_axs=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, mom0
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plotmom1(
        self,
        ck,
        chan_range=None,
        mom1clipping=0,
        ax=None,
        cbax=None,
        add_beam=False,
        savefits=False,
        fitsname=None,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the moment 1 (Intensity weighted mean velocity field).

        Parameters
        -----------
        ck : str
            Key of the cube to which the moment will be computed.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        mom1clipping : float
            Clipping to in order to compute the moment 1. Pixels with values
            smaller than the one given by clipping parameter will be masked
            with 0 values.
        ax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes` instance in which the position velodity
            diagram is drawn. If None, it will create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        savefits : bool
            If True, the moment will be saved in fits format.
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom1.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        mom1 : numpy.ndarray
            Intensity weighted velocity field.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        clipping = (
            float(mom1clipping.split("x")[0]) * self.sigma_noises[ck]
            if mom1clipping != 0
            else 0
        )
        add_beam = add_beam if "c" in ck else False
        mom1 = self.mom1(
            ck,
            chan_range=chan_range,
            clipping=clipping,
            savefits=savefits,
            fitsname=fitsname,
        )
        extent = (
            np.array(
                [
                    -(-0.5 - self.refpixs[ck][0]),
                    -(self.nxs - 0.5 - self.refpixs[ck][0]),
                    (-0.5 - self.refpixs[ck][1]),
                    (self.nys - 0.5 - self.refpixs[ck][1]),
                ]
            )
            * self.arcsecpix
        )
        fig, axs, cbax, velcmap = pl.plotmom1(
            mom1,
            extent=extent,
            interpolation=None,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Mean velocity [km/s]",
            return_fig_axs_velcmap=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, mom1
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plotmom2(
        self,
        ck,
        chan_range=None,
        mom2clipping=0,
        ax=None,
        cbax=None,
        add_beam=False,
        savefits=False,
        fitsname=None,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the moment 2 (velocity dispersion).

        Parameters
        -----------
        ck : str
            Key of the cube to which the moment will be computed.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        mom2clipping : float
            Clipping to in order to compute the moment 2. Pixels with values
            smaller than the one given by clipping parameter will be masked
            with 0 values.
        ax : matplotlib.axes.Axes, optional The matplotlib.axes.Axes` instance
            in which the position velodity diagram is drawn. If None, it will
            create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        savefits : bool
            If True, the moment will be saved in fits format
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_mom2.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        mom2 : numpy.ndarray
            Velocity dispersion
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        clipping = (
            float(mom2clipping.split("x")[0]) * self.sigma_noises[ck]
            if mom2clipping != 0
            else 0
        )
        add_beam = add_beam if "c" in ck else False
        mom2 = self.mom2(
            ck,
            chan_range=chan_range,
            clipping=clipping,
            savefits=savefits,
            fitsname=fitsname,
        )
        extent = (
            np.array(
                [
                    -(-0.5 - self.refpixs[ck][0]),
                    -(self.nxs - 0.5 - self.refpixs[ck][0]),
                    (-0.5 - self.refpixs[ck][1]),
                    (self.nys - 0.5 - self.refpixs[ck][1]),
                ]
            )
            * self.arcsecpix
        )
        fig, axs, cbax, velcmap = pl.plotmom2(
            mom2,
            extent=extent,
            interpolation=None,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Velocity dispersion [km/s]",
            return_fig_axs_velcmap=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, mom2
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plotmaxintens(
        self,
        ck,
        chan_range=None,
        ax=None,
        cbax=None,
        add_beam=False,
        savefits=False,
        fitsname=None,
        savefig=None,
        return_fig_axs_im=False,
        **kwargs,
    ):
        """
        Plots the moment 8 (peak intensity).

        Parameters
        -----------
        ck : str
            Key of the cube to which the moment diagram will be computed.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the moment. If None, the whole cube will be considered.
        ax : matplotlib.axes.Axes, optional The matplotlib.axes.Axes` instance
            in which the position velodity diagram is drawn. If None, it will
            create one. By default, None
        cbax : matplotlib.axes.Axes, optional
            The matplotlib.axes.Axes instance in which the color bar is drawn.
            If None, it will create one. By default, None
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        savefits : bool
            If True, the moments and the position velocity diagram will be
            saved in fits format
        fitsname : str
            Relative path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_maxintens.fits. If the path does not
            exist, it will be created.
        savefig : str, optional
            String of the full path to save the figure. If None, no figure is
            saved. By default, None.
        return_fig_axs_im : bool, optional
            If True, returns the figure, axes of the channel map, the axes the
            colorbar, and the image. If False, does not return anything.

        Returns:
        --------
        fig : matplotlib.figure.Figure
            Figure instance, only return_fig_axs_im=True
        ax : matplotlib.axes.Axes
            Axes of the channel map, only return_fig_axs_im=True
        cbax : tuple of matplotlib.axes.Axes
            Axes of the channel map and the colorbar, only returns if
            return_fig_axs_im=True.
        maxintens : numpy.ndarray
            Peak intensity
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        maxintens = self.maxintens(
            ck,
            chan_range=chan_range,
            savefits=savefits,
            fitsname=fitsname,
        )
        extent = (
            np.array(
                [
                    -(-0.5 - self.refpixs[ck][0]),
                    -(self.nxs - 0.5 - self.refpixs[ck][0]),
                    (-0.5 - self.refpixs[ck][1]),
                    (self.nys - 0.5 - self.refpixs[ck][1]),
                ]
            )
            * self.arcsecpix
        )
        fig, axs, cbax = pl.plotmaxintens(
            maxintens,
            extent=extent,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Peak " + self._getunitlabel(ck),
            return_fig_axs=True,
            **kwargs,
        )
        if return_fig_axs_im:
            return fig, axs, cbax, maxintens
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")


    def momentsandpv_and_params(
        self,
        ck,
        savefits=False,
        saveplot=False,
        mom1clipping=0,
        mom2clipping=0,
        verbose=True,
        chan_range=None,
        halfwidth_pv=0,
        add_beam=False,
        mom0values={v: None for v in ["vmax", "vcenter", "vmin"]},
        mom1values={v: None for v in ["vmax", "vcenter", "vmin"]},
        mom2values={v: None for v in ["vmax", "vcenter", "vmin"]},
        maxintensvalues={v: None for v in ["vmax", "vcenter", "vmin"]},
        pvvalues={v: None for v in ["vmax", "vcenter", "vmin"]},
    ):
        """
        Computes the moments and position velocity diagram including also the
        main parameters of the model listed in the first ax

        Parameters
        -----------
        ck : str
            Key of the cube to which the moments and the position velocity
            diagram will be computed.
        savefits : bool
            If True, the moments and the position velocity diagram will be
            saved in fits format
        saveplot : bool
            If True, a plot of the moments and position velocity diagrams will
            be saved
        mom1clipping : float
            Clipping to in order to compute the moment 1. Pixels with values
            smaller than the one given by clipping parameter will be masked
            with 0 values.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom
            corner.
        mom2clipping : float
            Clipping to in order to compute the moment 2. Pixels with values
            smaller than the one given by clipping parameter will be masked
            with 0 values.
        maxintensclipping : float
            Clipping to in order to compute the maximum value of the pixels
            along the velocity axis. Pixels with values smaller than the one
            given by clipping parameter will be masked with 0 values.
        mom0values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 0. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value
            of the moment image will be considered, respectively. Example:
            mom0values = {"vmax": None, "vcenter": None, "vmin": 0,}.
        mom1values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 1. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value
            of the moment image will be considered, respectively. Example:
            mom1values = {"vmax": None, "vcenter": None, "vmin": 0,}.
        mom2values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 2. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value
            of the moment image will be considered, respectively. Example:
            mom1values = {"vmax": None, "vcenter": None, "vmin": 0,}.
        maxintensvalues : dict or None
           Dictionary with the maximum, central, and minimum value to show in
           the plot of the maximum value along the velocity axis. If the
           dictionary value is None for vmax, vcenter, or vmin, then the
           maximum, central, or the minimum value of the moment image will be
           considered, respectively. Example: maxintensvalues = {"vmax": None,
           "vcenter": None, "vmin": None,}.
        """

        if verbose:
            print(
                """
\nComputing moments and the PV-diagram along the jet axis
"""
            )

        fig = plt.figure(figsize=(14, 10))
        gs = GridSpec(
            2,
            3,
            width_ratios=[1] * 3,  # + [0.85]*2,
            height_ratios=[1] * 2,
            hspace=0.3,
            wspace=0.25,
        )
        axs = {}
        cbaxs = {}
        gss = {}

        ik = "text"
        axs[ik] = plt.subplot(gs[0, 0])
        axs[ik].set_axis_off()

        ik = "mom0"
        gss[ik] = gs[0, 1].subgridspec(
            2,
            1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        axs[ik] = plt.subplot(gss[ik][1, 0])
        cbaxs[ik] = plt.subplot(gss[ik][0, 0])

        ik = "maxintens"
        gss[ik] = gs[0, 2].subgridspec(
            2,
            1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        axs[ik] = plt.subplot(gss[ik][1, 0])
        cbaxs[ik] = plt.subplot(gss[ik][0, 0])

        ik = "pv"
        gss[ik] = gs[1, 0].subgridspec(
            2,
            1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        axs[ik] = plt.subplot(gss[ik][1, 0])
        cbaxs[ik] = plt.subplot(gss[ik][0, 0])

        ik = "mom1"
        gss[ik] = gs[1, 1].subgridspec(
            2,
            1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        axs[ik] = plt.subplot(gss[ik][1, 0])
        cbaxs[ik] = plt.subplot(gss[ik][0, 0])

        ik = "mom2"
        gss[ik] = gs[1, 2].subgridspec(
            2,
            1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        axs[ik] = plt.subplot(gss[ik][1, 0])
        cbaxs[ik] = plt.subplot(gss[ik][0, 0])

        ak = "text"

        showtext = rf"""
        {self.modelname}
        Number of bowshocks: {self.nmodels}
        Tex = {self.Tex.value} K
        $i = {{{ut.list2str(self.ies*180/np.pi)}}}^\circ$
        $v_\mathrm{{sys}} = {self.vsys}$ km/s
        $v_\mathrm{{j}} = {{{ut.list2str(self.vjs)}}}$ km/s
        $v_0 = {{{ut.list2str(self.v0s)}}}$ km/s
        $v_a = {{{ut.list2str(self.vas)}}}$ km/s
        $L_0 = {{{ut.list2str(self.L0s)}}}$ arcsec
        $z_\mathrm{{j}} = {{{ut.list2str(self.zjs)}}}$ arcsec
        $r_\mathrm{{b,f}} = {{{ut.list2str(self.rbfs)}}}$ arcsec
        $t_\mathrm{{j}} = {{{ut.list2str(self.tjs)}}}$ yr
        mass $= {{{ut.list2str(self.masss*10**4)}}}\times 10^{{-4}}$ M$_\odot$
        $\rho_a = {{{ut.list2str(self.rhoas*10**20)}}}\times 10^{{-20}}$ g cm$^{{-3}}$
        $\dot{{m}}_0 = {{{ut.list2str(self.m0s*10**6)}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
        $\dot{{m}}_{{a,f}} = {{{ut.list2str(self.mwfs*10**6)}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
        """
        for n, line in enumerate(showtext.split("\n")):
            axs["text"].text(
                0,
                0.99 - 0.06 * n,
                line,
                fontsize=12 - self.nmodels,
                transform=axs["text"].transAxes,
            )

        add_beam = add_beam if "c" in ck else False
        ak = "mom0"
        self.plotmom0(
            ck=ck,
            chan_range=chan_range,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_fig_axs_im=False,
            add_beam=add_beam,
            **mom0values,
        )

        ak = "mom1"
        self.plotmom1(
            ck=ck,
            chan_range=chan_range,
            mom1clipping=mom1clipping,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_fig_axs_im=False,
            add_beam=add_beam,
            **mom1values,
        )

        ak = "mom2"
        self.plotmom2(
            ck=ck,
            chan_range=chan_range,
            mom2clipping=mom2clipping,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_fig_axs_im=False,
            add_beam=add_beam,
            **mom2values,
        )

        ak = "pv"
        self.plotpv(
            ck=ck,
            halfwidth=halfwidth_pv,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_fig_axs_im=False,
            **pvvalues,
        )

        ak = "maxintens"
        self.plotmaxintens(
            ck=ck,
            chan_range=chan_range,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_fig_axs_im=False,
            add_beam=add_beam,
            **maxintensvalues,
        )

        if saveplot:
            fig.savefig(
                f"models/{self.modelname}/momentsandpv_and_params_{ck}.pdf",
                bbox_inches="tight",
            )

    def momentsandpv_and_params_all(self, **kwargs):
        """
        Computes all the moments and pv to the cubes listed in self.listmompvs,
        including a list of values of the main parameters of the model in the
        first ax
        """
        for ck in self.listmompvs:
            self.momentsandpv_and_params(ck, **kwargs)
