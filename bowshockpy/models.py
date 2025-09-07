"""This module contains the implementation of the momentum conserving bowshock
model"""

import astropy.units as u
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

import bowshockpy.plots as pl
from bowshockpy.version import __version__


class BaseModel:
    """
    Parent class of Models.
    """
    def stoyr(self, value):
        """
        Converts seconds to years

        Parameters
        ----------
        value : float
            Seconds to convert to years

        Returns
        -------
        float
            Years
        """
        return value * u.s.to(u.yr)

    def solMasskm2togcm2(self, value):
        """
        Converts solar masses/km**2 to g/cm**2

        Parameters
        ----------
        value : float
            Solar masses per km^2 to convert to g/cm^2

        Returns
        -------
        float
            g/cm2
        """
        return value * (u.solMass / u.km**2).to(u.g / u.cm**2)

    def solMasskm3togcm3(self, value):
        """
        Converts solar masses/km**3 to g/cm**3

        Parameters
        ----------
        value : float
            Solar masses per km^3 to convert to g/cm^3

        Returns
        -------
        float
            g/cm3
        """
        return value * (u.solMass / u.km**3).to(u.g / u.cm**3)

    def km2arcsec(self, value):
        """
        Converts km to arcsec

        Parameters
        ----------
        value : float
            km to convert to arcsec

        Returns
        -------
        float
            g/cm3
        """
        return value * u.km.to(u.au) / self.distpc


class BowshockModel(BaseModel):
    """
    Bowshock model for a negligible internal working surface radius.

    Parameters
    -----------
    L0 : float
        Characteristic scale of the bowshock [km]
    zj : float
        Distance between the source and the internal working surface [km]
    vj : float
        Velocity of the internal working surface [km/s]
    va : float
        Velocity of the ambient gas surrounding the jet [km/s]
    v0 : float
        Velocity at which the material is ejected sideways from the internal
        working surface [km/s]
    mass : float
        Total mass of the bowshock shell [Solar masses]
    distpc : float
        Distance between the source and the observer [pc]
    rbf_obs: float, optional
        Final radius of the bowshock [km]. If None, the theoretical final
        radius is calculated.

    Attributes:
    -----------
    rbf : float
        Final radius of the bowshock [km]
    zbf : float
        z coordinate at the final radius rbf [km]
    tj : float
        Dynamical time defined as zj / vj [s].
    tj_yr : float
        Dynamical time in yrs.
    rhoa : float
        Ambient density [Msun km^3]
    rhoa_gcm3 : float
        Ambient density [g cm^-3]
    mp0 : float
        Mass rate at which the jet material is ejected sideways from the
        internal working surface [Msun / s]
    mp0_solmassyr : float
        Mass rate at which the jet material is ejected sideways from the
        internal working surface [Msun / yr]
    mpamb_f : float
        Mass rate of ambient material is incorporated into the bowshock
        shell [Msun / s]
    mpamb_f_solmassyr : float
        Mass rate of ambient material is incorporated into the
        bowshock shell [Msun / yr]
    zj_arcsec : float
        Distance from the source to the internal working surface [arcsec]
    L0_arcsec : float
        Bowshock characteristic scale [arcsec]
    rbf_arcsec : float
        Final radius of the bowshock [arcsec]
    zbf_arcsec : float
        z-coordinate at the final radius [arcsec]

    References:
    -----------
    [1] Tabone, B., Raga, A., Cabrit, S. & Pineau des Forêts, G. "Interaction
    between a pulsating jet and a surrounding disk wind. A hydrodynamical
    perspective." Astron. Astrophys. 614, A119 (2018).

    [2] Ostriker, E. C., Lee, C.-F., Stone, J. M. & Mundy, L. G. A Ballistic
    Bow Shock Model for Jet-driven Protostellar Outflow Shells. Astrophys. J.
    557, 443–450 (2001).

    [3] Blazquez-Calero, G., et al. (in prep)

    """

    default_kwargs = {
        "rbf_niters": 1000,
    }

    def __init__(self, L0, zj, vj, va, v0, mass, distpc, rbf_obs=None, **kwargs):
        # for param in ps:
        #     setattr(self, param, ps[param])
        self.L0 = L0
        self.zj = zj
        self.vj = vj
        self.va = va
        self.v0 = v0
        self.mass = mass
        self.rbf_obs = rbf_obs
        self.distpc = distpc
        # TODO remove kwargs
        for kwarg in self.default_kwargs:
            kwarg_attr = (
                kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            )
            setattr(self, kwarg, kwarg_attr)
        if self.rbf_obs is None:
            self.rbf = self.rbf_calc()
        else:
            self.rbf = self.rbf_obs
        self.zbf = self.zb_r(self.rbf)

        self.tj = self.zj / self.vj
        self.tj_yr = self.stoyr(self.tj)
        self.rhoa = self.rhoa_fromintmass_analytical(self.rbf, self.mass)
        self.rhoa_gcm3 = self.solMasskm3togcm3(self.rhoa)
        self.mp0 = self.mp0_calc(self.rhoa)
        self.mp0_solmassyr = self.mp0_calc(self.rhoa) / self.stoyr(1)
        self.mpamb_f = self.mpamb_f_calc(self.rhoa)
        self.mpamb_f_solmassyr = self.mpamb_f / self.stoyr(1)
        self.zj_arcsec = self.km2arcsec(self.zj)
        self.L0_arcsec = self.km2arcsec(self.L0)
        self.rbf_arcsec = self.km2arcsec(self.rbf)
        self.zbf_arcsec = self.km2arcsec(self.zbf)

    def gamma(self):
        """
        Computes the gamma parameter, defined as gamma = (vj - va) / v0

        Returns
        -------
        float
            gamma parameter
        """
        return (self.vj - self.va) / self.v0

    def rb(self, zb):
        """
        Bowshock radius for a given z coordinate of the bowshock

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Bowshock radius at zb [km]
        """
        return (self.L0**2 * (self.zj - zb)) ** (1 / 3)

    def vr(self, zb):
        """
        Computes the transversal component of the velocity (along the r-axis,
        perpendicular to the symmetry axis of the bowshock)

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Transversal component of the velocity [km/s]
        """
        return self.v0 * (1 + 3 * self.rb(zb) ** 2 / self.gamma() / self.L0**2) ** (-1)

    def vz(self, zb):
        """
        Computes the longitudinal component of the velocity (along z-axis, the
        symmetry axis of the bowshock)

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Longitudinal component of the velocity [km/s]
        """
        return self.va + (self.vj - self.va) * (
            1 + 3 * self.rb(zb) ** 2 / self.gamma() / self.L0**2
        ) ** (-1)

    def vtot(self, zb):
        """
        Computes the total speed of the velocity

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Speed of the bowshock at zb [km/s]
        """
        return np.sqrt(self.vr(zb) ** 2 + self.vz(zb) ** 2)

    def alpha(self, zb):
        """
        Computes the angle between the bowshock axis and the velocity

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Angle between the bowshock axis and the velocity [radians]
        """
        return np.arctan(self.vr(zb) / self.vz(zb))

    def alpha2(self, zb):
        """
        Computes the angle between the bowshock axis and the local tangent to
        the shell surface at the z-coordinate of the bowshock zb

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Angle between the bowshock axis and the local tangent to the shell
            surface at the z-coordinate of the bowshock zb [radians]
        """
        return np.arcsin((1 + 9 * self.rb(zb) ** 4 / self.L0**4) ** (-0.5))

    def alpha2_rb(self, rb):
        """
        Computes the angle between the bowshock axis and the local tangent to
        the shell surface at the z-coordinate of the bowshock zb

        Parameters
        ----------
        rb : float
            r coordinate of the bowshock [km]

        Returns
        -------
        float
            Angle between the bowshock axis and the local tangent to the shell
            surface at the z-coordinate of the bowshock zb [radians]
        """
        return np.arcsin((1 + 9 * rb**4 / self.L0**4) ** (-0.5))

    def theta(self, zb):
        """
        Computes the polar angle of the position vector

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            Polar angle of the position vector in radians
        """
        return np.arctan(self.rb(zb) / zb)

    def rbf_0(self, rr):
        """
        This is the target function to minimize by rbf_calc() to find the
        theoretical final radius

        Parameters
        ----------
        rr : float
            r coordinate of the bowshock [km]

        Returns
        -------
        float
            Value to minimize (it will be zero when rr is the final radius of
            the bowshock)
        """
        return (
            1 / self.gamma() * (rr / self.L0) ** 3
            + rr / self.L0
            - self.v0 * self.zj / self.L0 / self.vj
        )

    def rbf_calc(self, ns=None, use_minimize=True):
        """
        Computes numerically the bowshock final radius

        Parameters
        ----------
        use_minimize : bool
            If True, scipy.optimize.minimize_scalar is used. If False, rbf is
            computed by brute force
        ns : int
            number of iteractions to compute rbf by brute force if
            use_minimize=False

        Returns:
        --------
        float
            Solution for the final radius of the bowshock [km]
        """
        if use_minimize:
            bounds = (0, self.rb(0))
            return minimize_scalar(
                lambda x: np.abs(self.rbf_0(x)), method="bounded", bounds=bounds
            ).x
        else:
            ns = self.rbf_niters if ns is None else ns
            rrs = np.linspace(0, self.rb(0), ns)
            trials = np.array([np.abs(self.rbf_0(rr)) for rr in rrs])
            return rrs[np.argmin(trials)]

    def zb_r(self, rr):
        """
        Bowshock z coordinate for a given radius of the bowshock

        Parameters
        ----------
        rr : float
            radius of the bowshock [km]

        Returns
        -------
        float
            z coordinate of the bowshock [km]
        """

        return self.zj - rr**3 / self.L0**2

    # def surfdens(self, rr):
    #     return rhoa/2/rb * (ps["L0"]**2*gamma(ps)/3 + rb**2)**2 \
    #        / (rb**2*np.cos(alpha2(rb,ps)) + (ps["L0"]**2/3)*np.sin(alpha2(rb,ps)))

    def surfdens(self, zb):
        """
        Computes the surface density of the bowshock

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]

        Returns
        -------
        float
            surface density at zb [solmass/km^2]
        """
        cosa = np.cos(self.alpha2(zb))
        tana = np.tan(self.alpha2(zb))
        sd = 0.5 * self.rhoa * cosa * (self.gamma() * tana + 1) ** 2 * self.rb(zb)
        return sd

    def dr_func(self, zb, dz):
        """
        Differential of r given a differential of z

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        dz : float
            increment in z [km]

        Returns
        -------
        float
            differential of r [km]
        """
        return 1 / 3 * (self.L0 / self.rb(zb)) ** 2 * dz

    def dz_func(self, zb, dr):
        """
        Differential of r given a differential of z

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        dr : float
            increment in radius [km]

        Returns
        -------
        float
            differential of r [km]
        """
        return 3 * (self.rb(zb) / self.L0) ** 2 * dr

    def dsurf_func(self, zb, dz, dphi):
        """
        Differential of surface given a differential in z and phi (azimuthal
        angle)

        Parameters
        -----------
        zb : float
            z coordinate of the bowshock [km]
        dz : float
            increment in radius [km]
        dphi : float
            increment in azimuthal angle [radians]

        Returns
        -------
        float
            differential of surface [Msun/km^2]
        """
        # sina = np.sin(self.alpha2(zb))
        sina = (1 + 9 * self.rb(zb) ** 4 / self.L0**4) ** (-0.5)
        dr = self.dr_func(zb, dz)
        return self.rb(zb) * dr * dphi / sina

    def dmass_func(self, zb, dz, dphi):
        """
        Differential of mass given a differential in z and phi (azimuthal angle)

        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        dz : float
            increment in radius [km]
        dphi : float
            increment in azimuthal angle [radians]

        Returns
        -------
        float
            differential of mass [Msun]
        """
        return self.surfdens(zb) * self.dsurf_func(zb, dz, dphi)

    def intmass_analytical(self, rbf):
        """
        Computes the total mass of the bowshock shell

        Parameters
        ----------
        rbf : float
            final radius of the bowshock [km]

        Returns
        -------
        float
            Total mass of the shell [Msun]
        """
        uu = rbf / self.L0 * (3 / self.gamma()) ** 0.5
        analit_int = uu**5 / 5 + 2 * uu**3 / 3 + uu
        massint = (
            self.rhoa
            * (self.L0 / np.sqrt(3)) ** 3
            * np.pi
            * self.gamma() ** (5 / 2)
            * analit_int
        )
        return massint

    def intmass_numerical(self, r0, rbf, return_residual=False):
        """
        Computes numerically the total mass of the bowshock shell in a range of
        radius from r0 to rbf

        Parameters
        ----------
        r0 : float
            initial radius of the bowshock [km]
        rbf : float
            final radius of the bowshock [km]

        Returns
        -------
        float
            Total mass of the shell [Msun]
        """

        def integrand(rb):
            tana = np.tan(self.alpha2_rb(rb))
            return (self.gamma() * tana + 1) ** 2 / tana * rb**2

        integ = quad(integrand, r0, rbf)
        massint = self.rhoa * np.pi * integ[0]
        if return_residual:
            return massint, integ[1]
        else:
            return massint

    def rhoa_fromintmass_analytical(self, rb, massint):
        """
        Computes the ambient density given the integrated mass of the bowshock
        at rb

        Parameters
        ----------
        rb : float
            radius of the bowshock [km]
        massint : float
            integrated mass of the bowshock shell up to rb [Msun]

        Returns
        -------
        float
            Density of the ambient [Msun/km^3]
        """

        uu = rb / self.L0 * (3 / self.gamma()) ** 0.5
        analit_int = uu**5 / 5 + 2 * uu**3 / 3 + uu
        rhoa = massint * (
            (self.L0 / np.sqrt(3)) ** 3 * np.pi * self.gamma() ** (5 / 2) * analit_int
        ) ** (-1)
        return rhoa

    def rhoa_fromintmass_sigma_simple(self, R0, Rb, massint, return_residual=False):
        """
        Computes numerically the ambient density taken into account the
        integrated mass in a range of radii from R0 to Rb

        Parameters
        ----------
        R0 : float
            initial radius of the bowshock [km]
        Rb : float
            radius of the bowshock [km]
        massint : float
            integrated mass of the bowshock shell from R0 to Rb [Msun]

        Returns
        -------
        float
            Density of the ambient [Msun/km^3]
        """

        def integrand(rb):
            tana = np.tan(self.alpha2_rb(rb))
            return (self.gamma() * tana + 1) ** 2 / tana * rb**2

        integ = quad(integrand, R0, Rb)
        rhoa = massint / np.pi / integ[0]
        if return_residual:
            return rhoa, integ[1]
        else:
            return rhoa

    def mp0_calc(self, rhoa):
        """
        Computes the mass rate at which the jet material is ejected sideways
        from the internal working surface

        Parameters
        ----------
        rhoa : float
            Ambient density

        Returns
        -------
        float
            Mass rate at which the jet material is ejected sideways from the
            internal working surface [Msun/s]
        """
        mp0 = rhoa * np.pi * self.L0**2 * (self.vj - self.va) ** 2 / 3 / self.v0
        return mp0

    def mpamb_f_calc(self, rhoa):
        """
        Computes the mass rate of ambient material incorporated into the
        bowshock shell

        Parameters
        ----------
        rhoa : float
            Ambient density

        Returns
        -------
        float
            Mass rate of the ambient incorporated into the bowshock [Msun/s]
        """
        mpamb_f = np.pi * rhoa * (self.vj - self.va) * self.rbf**2
        return mpamb_f

    def get_modelplot(
        self,
        modelname="none",
        nzs=200,
        figsize=(16, 3),
        narrows=10,
        v_arrow_ref=100,
        linespacing=0.08,
        textbox_widthratio=0.7,
        **kwargs,
    ):
        """
        Plot a figure including the main parameters of the bowshock model, its
        morphology and kinematics, and the distribution of the surface density

        Parameters
        ----------
        modelname : str, optional
            Name of the model to include in the plot
        nzs : int, optional
            Number of points used to compute the model solutions
        figsize : tuple, optional
            Tuple passed to `matplotib.pyplot.figure` to define the
            dimensions of the figure
        narrows : int, optional
            Number of arrows to show in order to indicate the velocity at
            each symmetrical half of the model.
        v_arrow_ref : float, optional
            Velocity in km/s to use as reference in the reference arrow
        linespacing : float, optional
            Spacing between the text lines
        textbox_widthratio : float, optional
            Width ratio of the text ax to pass to GridSpec
        kwargs : optional
            Keyword arguments into `~bowshockpy.plot.BowshockModelPlot`

        Returns
        --------
        modelplot : `~bowshockpy.plot.BowshockModelPlot` class instance
            An instance of a class BowshockModelPlot, which contains
            information on the figure and the model data
        """
        modelplot = pl.BowshockModelPlot(
            self,
            modelname=modelname,
            nzs=nzs,
            figsize=figsize,
            narrows=narrows,
            v_arrow_ref=v_arrow_ref,
            linespacing=linespacing,
            textbox_widthratio=textbox_widthratio,
            **kwargs,
        )
        return modelplot
