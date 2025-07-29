import numpy as np

from scipy.optimize import minimize_scalar
from scipy.integrate import quad
from scipy.ndimage import rotate

import astropy.units as u
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from datetime import datetime

import sys

import bowshockpy.utils as ut
import bowshockpy.radtrans as rt
import bowshockpy.moments as moments
import bowshockpy.plots as pl
from bowshockpy.version import __version__


class BowshockModel():
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
        Final radius of the bowshock [km]. If None, the theoretical final radius
        is calculated.
    
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
     
    [2] Ostriker, E. C., Lee, C.-F., Stone, J. M. & Mundy, L. G. A Ballistic Bow
    Shock Model for Jet-driven Protostellar Outflow Shells. Astrophys. J. 557,
    443–450 (2001).
     
    [3] Blazquez-Calero, G., et al. (in prep)

    """

    default_kwargs = {
        "rbf_niters": 1000,
    }

    def __init__(
            self, L0, zj, vj, va, v0, mass, 
            distpc, rbf_obs=None, **kwargs):
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
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs \
            else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)
        if self.rbf_obs is None:
            self.rbf = self.rbf_calc()
        else:
            self.rbf = self.rbf_obs
        self.zbf = self.zb_r(self.rbf)

        self.tj = self.zj / self.vj
        self.tj_yr = self.stoyr(self.tj)
        self.rhoa = self.rhoa_fromintmass_analytical(
            self.rbf,
            self.mass
        )
        self.rhoa_gcm3 = self.solMasskm3togcm3(self.rhoa)
        self.mp0 = self.mp0_calc(self.rhoa)
        self.mp0_solmassyr = self.mp0_calc(self.rhoa) / self.stoyr(1)
        self.mpamb_f = self.mpamb_f_calc(self.rhoa)
        self.mpamb_f_solmassyr = self.mpamb_f / self.stoyr(1)
        self.zj_arcsec = self.km2arcsec(self.zj)
        self.L0_arcsec = self.km2arcsec(self.L0)
        self.rbf_arcsec = self.km2arcsec(self.rbf)
        self.zbf_arcsec = self.km2arcsec(self.zbf)
 
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
        return value * (u.solMass/u.km**2).to(u.g/u.cm**2)

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
        return value * (u.solMass/u.km**3).to(u.g/u.cm**3)

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

    def gamma(self):
        """
        Computes the gamma parameter, defined as gamma = (vj - va) / v0

        Returns
        -------
        float
            gamma parameter
        """
        return (self.vj-self.va) / self.v0

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
        return (self.L0**2*(self.zj-zb))**(1/3)

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
        return self.v0*(1 + 3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

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
        return self.va + (self.vj-self.va)*(1+3*self.rb(zb)**2/self.gamma()/self.L0**2)**(-1)

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
        return np.sqrt(self.vr(zb)**2 + self.vz(zb)**2)

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
        return np.arcsin((1+9*self.rb(zb)**4/self.L0**4)**(-0.5))

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
        return np.arcsin((1+9*rb**4/self.L0**4)**(-0.5))

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
        return 1 / self.gamma() * (rr/self.L0)**3 + rr/self.L0 - self.v0*self.zj/self.L0/self.vj

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
            return minimize_scalar(lambda x: np.abs(self.rbf_0(x)),
                                   method="bounded", bounds=bounds).x
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
        sd = 0.5 * self.rhoa * cosa * (self.gamma()*tana+1)**2 * self.rb(zb)
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
        return 1/3 * (self.L0 / self.rb(zb))**2 * dz

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
        return 3 * (self.rb(zb) / self.L0)**2 * dr

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
        sina = (1+9*self.rb(zb)**4/self.L0**4)**(-0.5)
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
        uu = rbf / self.L0 * (3/self.gamma())**0.5
        analit_int = uu**5 / 5 + 2*uu**3/3 + uu
        massint = self.rhoa * (self.L0/np.sqrt(3))**3 * np.pi * self.gamma()**(5/2) * analit_int
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
            return (self.gamma()*tana+1)**2 / tana * rb**2
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
 
        uu = rb / self.L0 * (3/self.gamma())**0.5
        analit_int = uu**5 / 5 + 2*uu**3/3 + uu
        rhoa = massint * ((self.L0/np.sqrt(3))**3 * np.pi * self.gamma()**(5/2) * analit_int)**(-1)
        return rhoa

    def rhoa_fromintmass_sigma_simple(
            self, R0, Rb, massint, return_residual=False
            ):
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
            return (self.gamma()*tana+1)**2 / tana * rb**2

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
        mp0 = rhoa * np.pi * self.L0**2 * (self.vj-self.va)**2 / 3 / self.v0
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

    def get_modelplot(self, **kwargs):
        """
        Plot a figure including the main parameters of the bowshock model, its
        morphology and kinematics, and the distribution of the surface density
        
        Parameters
        -----------
        kwargs : optional
            Keyword arguments into `~bowshockpy.plot.BowshockModelPlot`

        Returns
        --------
        modelplot : `~bowshockpy.plot.BowshockModelPlot` class instance
            An instance of a class BowshockModelPlot, which contains information
            on the figure and the model data
        """
        modelplot = pl.BowshockModelPlot(self, **kwargs)
        return modelplot


class ObsModel(BowshockModel):
    """
    Computes the projected morphology and kinematics of a BowshockModel model

    Parameters
    -----------
    model : class instance
        instance of BowshockModel model to get the attributes
    i_deg : float
        Inclination angle between the bowshock axis and the line-of-sight
        [degrees] 
    pa_deg : float, optional
        Position angle, default 0 [degrees]
    vsys : float, optional
        Systemic velocity of the source, default 0 [km/s]
    """
    def __init__(self, model, i_deg, pa_deg=0, vsys=0, **kwargs):
        self.__dict__ = model.__dict__
        self.i_deg = i_deg
        self.i = i_deg * np.pi/180
        self.pa_deg = pa_deg
        self.pa = pa_deg * np.pi/180
        self.vsys = vsys
        # for param in model.__dict__:
        #     setattr(self, param, getattr(model, param))
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs \
            else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

    def vzp(self, zb, phi):
        """
        Calculates the line-of-sight velocity for a point of the bowshock shell
        with (zb, phi)
        
        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        phi : float
            azimuthal angle [radians]

        Returns
        -------
        float
            Line-of-sight velocity [km/s]
        """
        a = self.alpha(zb)
        return self.vtot(zb) * (
            np.cos(a)*np.cos(self.i) - np.sin(a)*np.cos(phi)*np.sin(self.i)
            )

    def xp(self, zb, phi):
        """
        Calculates the xp coordinate for a point of the bowshock shell
        with (zb, phi)
        
        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        phi : float
            azimuthal angle [radians]

        Returns
        -------
        float
            xp coordinate in the plane-of-sky [km]
        """
        return self.rb(zb)*np.cos(phi)*np.cos(self.i) + zb*np.sin(self.i)

    def yp(self, zb, phi):
        """
        Calculates the yp coordinate for a point of the bowshock shell
        with (zb, phi)
        
        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        phi : float
            azimuthal angle [radians]

        Returns
        -------
        float
            yp coordinate in the plane-of-sky [km]
        """
        return self.rb(zb) * np.sin(phi)

    def zp(self, zb, phi):
        """
        Calculates the xp coordinate for a point of the bowshock shell
        with (zb, phi)
        
        Parameters
        ----------
        zb : float
            z coordinate of the bowshock [km]
        phi : float
            azimuthal angle [radians]

        Returns
        -------
        float
            zp coordinate, along the line-of-sight direction [km]
        """
        return -self.rb(zb)*np.cos(phi)*np.sin(self.i) + zb*np.cos(self.i)


    def get_modelplot(self, **kwargs):
        """
        Plot a figure including the main parameters of the bowshock model, its
        morphology and kinematics, and the distribution of the surface density
        
        Parameters
        -----------
        kwargs : optional
            Keyword arguments into `~bowshockpy.plot.BowshockModelPlot`

        Returns
        --------
        modelplot : `~bowshockpy.plot.BowshockModelPlot` class instance
            An instance of a class BowshockModelPlot, which contains information
            on the figure and the model data
        """
        modelplot = pl.BowshockObsModelPlot(self, **kwargs)
        return modelplot



class BowshockCube(ObsModel):
    """
    Computes the spectral cube of the bowshock model

    Parameters
    -----------
    obsmodel : class instance
        Instance of ObsModel
    nphis : int
        Number of azimuthal angles phi to calculate the bowshock solution
    vch0 : float
        Central velocity of the first channel map [km/s]
    vchf : float
        Central velocity of the last channel map [km/s]
    xpmax : float
        Physical size of the channel maps along the x axis [arcsec]
    nzs : int, optional 
        Number of points used to compute the model solutions
    nc : int, optional
        Number of spectral channel maps
    nxs : int, optional
        Number of pixels in the right ascension axis.
    nys : int, optional
        Number of pixels in the declination axis. 
    refpix : list | None, optional
        Pixel coordinates (zero-based) of the source, i.e., the origin from
        which the distances are measured. The first index is the R.A. axis, the
        second is the  Dec. axis [[int, int] or None] 
    CIC : bolean, optional
        Set to True to perform Cloud in Cell interpolation [1].
    vt : str | float, optional
        Thermal+turbulent line-of-sight velocity dispersion [km/s] If
        thermal+turbulent line-of-sight velocity dispersion is smaller than the
        instrumental spectral resolution, vt should be the spectral resolution.
        It can be also set to a integer times the channel width (e.g.,
        "2xchannel")
    tolfactor_vt : float, optional
        Neighbour channel maps around a given channel map with vch will stop
        being populated when their difference in velocity with respect to vch is
        higher than this factor times vt. The lower the factor, the quicker will
        be the code, but the total mass will be underestimated. If vt is not
        None, compare the total mass of the output cube with the mass parameter
        that the user has defined 
    verbose : bolean, optional
        Set True to verbose messages about the computation
    kwargs : optional
        Keyword arguments into `~bowshockpy.plot.plot_channel`

    Attributes:
    -----------
    nrs : int
        Number of model points to which the solution has been computed.
    rs : numpy.ndarray
        Array of the radii of the model.
    dr : float
        Increment of radii between the points, which is constant.
    zs : numpy.ndarray
        Array of the z-coordinates of the model.
    dzs : numpy.ndarray
        Increment of z-coordinates between the points.
    phis : numpy.ndarray
        Array of the azimuthal angles of the model.
    dphi : float
        Increment in azimuthal angle of the points of the model.
    vs : numpy.ndarray
        Array with the velocities of the points of the model.
    velchans : numpy.ndarray
        Array with the line-of-sight velocities of the channels of the spectral
        cube.
    cube : numpy.ndarray
        Spectral cube of the masses of the bowshock model.

    References:
    -----------
    [1] Fehske, H., Schneider, R., & Weiße, A. (2008), Computational
    Many-Particle Physics, Vol.  739 (Springer), doi: 10.1007/978-3-540-74686-7.

    """

    def __init__(
            self, obsmodel, nphis, vch0, vchf, xpmax, nzs=200, nc=50,
            nxs=200, nys=200, refpix=[0,0], CIC=True, vt="2xchannel",
            tolfactor_vt=None, verbose=True, **kwargs):
        self.__dict__ = obsmodel.__dict__
        self.nphis = nphis
        self.vch0 = vch0
        self.vchf = vchf
        self.xpmax = xpmax
        self.nzs = nzs
        self.nc = nc
        self.nxs = nxs
        self.nys = nys
        self.refpix = refpix
        self.CIC = CIC
        self.vt = vt
        self.tolfactor_vt = tolfactor_vt
        self.verbose = verbose
        self._calc_params_init()
        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs \
                else self.default_kwargs[kwarg]
            setattr(self, kwarg, kwarg_attr)

        self.nrs = None
        self.rs = np.array([])
        self.dr = None
        self.zs = np.array([])
        self.dzs = np.array([])

        self.phis = np.array([])
        self.dphi = None

        self.vs = np.array([])
        self.velchans = np.array([])

        self._fromcube_mass = 0
        self.cube = None
        self.cube_sampling = None

    def _DIMENSION_ERROR(self, fromcube):
            sys.exit(f"""
ERROR: The provided cube into which the model is to be build has dimensions
{np.shape(fromcube)} but the dimensions of the desired model cube is {(self.nc,
self.nys, self.nxs)}. Please, provide a cube with the right dimensions or do not
provide any cube.
""")
  
    def _OUTSIDEGRID_WARNING(self,):
        print("""
WARNING: Part of the model lie outside the grid of the spectral cube! The model
will be truncated or not appearing at all in your spectral cube. This is due to
at least one of three reasons: 
    - The image is too small. Try to make the image larger by increasing the
    number of pixels (parameters nxs and nys), or increase the physical size of
    the image (parameter xpmax).
    - The model is far away from the image center. Try to change the reference
    pixel where the physical center (the source) is found (parameter refpix).
    - The model is outside your velocity coverage. Try to change the range of
    velocity channels of the spectral cube (parameters vch0 and vchf, consider
    negative floats if the model is blueshifted).\n
""")

    def _MASS_CONSISTENCY_WARNING(self, massloss):
        print(rf"""
WARNING: The integrated mass of the cube is {massloss:.1e} % less than the input
total mass of the bowshock. This can be due to several factors:
    - Part of the model lie outside the grid of the spectral cube. If this is
    not intended, try to solve it by making the maps larger, changing the
    reference pixel to center the model in the maps, or increasing the velocity
    coverage of the spectral cube.  
    - The difference between the integrated mass of the cube and the input total
    mass of the bowshock model is due to numerical errors. If you think that the
    difference is too big, you can reduce it by increasing the number of points
    of the model (inceasng nzs or/and nphis parameters).
    - The masses corresponding to a channel maps are spread along the cube in
    the velocity axis following a Gaussian distribution, being sigma equal to vt
    parameter. This distribution is truncated at vt*tolfactor_vt in order to
    make the computation substatially faster, but it can result in an
    underestimation of the integrated mass of the spectral cube. Try to make
    tolfactor_vt larger.
""")

    def _SAMPLING_XY_WARNING(self,):
        print("""
WARNING: It is possible that the model is not well sampled in the plane of sky
given the cube dimensions and the number of model points. You can ensure a
better sampling by increasing the number of model points (nzs
parameter) or decreasing the pixel size (nxs and nys parameters).
""")

    def _SAMPLING_V_WARNING(self,):
        print("""
WARNING: It is possible that the model is not well sampled in the velocity direction 
given the cube dimensions and the number of model points. You can ensure a
better sampling by increasing the number of model points in the azimuthal direction (nphis parameter) or decreasing the pixel size (nxs and nys parameters).
""")

    def _SAMPLING_PHI_WARNING(self,):
        print("""
WARNING: It is possible that the model is not well sampled in the plane of sky
given the cube dimensions and the number of azimuthal points of the model
(nphis). You can ensure a better sampling by increasing the number of model
points in the azimuthal direction (nphis parameter) or decreasing the pixel
size (nxs and nys parameters).
""")

    def _calc_params_init(self):
        self.chanwidth = (self.vchf - self.vch0) / (self.nc-1)
        self.abschanwidth = np.abs(self.chanwidth)
        self.vt = self.vt if type(self.vt)!=str \
              else float(self.vt.split("x")[0])*self.chanwidth
        self.arcsecpix = self.xpmax / float(self.nxs)
        if self.refpix == None:
            if self.nxs%2 == 0:
                xref = int(self.nxs / 2)
            else:
                xref = int((self.nxs-1) / 2)
            if self.nys%2 == 0:
                yref = int(self.nys / 2)
            else:
                yref = int((self.nys-1) / 2)
            self.refpix = [xref, yref]

    def _cond_populatechan(self, diffv):
        if self.tolfactor_vt is not None:
            return diffv < np.abs(self.vt)*self.tolfactor_vt
        else:
            return True

    def _wvzp(self, diffv, dmass):
        """
        Weight the masses across the velocity axis using a Gaussian
        distribution
        """
        normfactor = np.abs(self.chanwidth) / (np.sqrt(np.pi)*np.abs(self.vt))
        em = dmass * np.exp(-(diffv/self.vt)**2) * normfactor
        return em

    # def _wxpyp(self, chan, vchan, xpix, ypix, dxpix, dypix, vzp, dmass):
    def _doCIC(self, chan, diffv, xpix, ypix, dxpix, dypix, dmass):
        """Cloud In Cell method"""
        em = self._wvzp(diffv, dmass)
        self.cube[chan, ypix, xpix] += em * (1-dxpix) * (1-dypix)
        self.cube[chan, ypix, xpix+1] += em * dxpix * (1-dypix)
        self.cube[chan, ypix+1, xpix] += em * (1-dxpix) * dypix
        self.cube[chan, ypix+1, xpix+1] += em * dxpix * dypix
    
    def _doNGP(self, chan, diffv, xpix, ypix, dxpix, dypix, dmass):
        """Nearest Grid Point method"""
        em = self._wvzp(diffv, dmass)
        self.cube[chan, ypix, xpix] += em 

    def _sampling(self, chan, xpix, ypix):
        self.cube_sampling[chan, ypix, xpix] += 1
    
    def _check_mass_consistency(self, return_isconsistent=False):
        print("Checking total mass consistency...")
        intmass_cube = np.sum(self.cube)
        intmass_model = self.mass+self._fromcube_mass
        mass_consistent = np.isclose(
            intmass_cube, intmass_model)
        massloss = (intmass_model-intmass_cube) / self.mass * 100
        if mass_consistent:
            print(rf"""
Mass consistency test passed: The input total mass of the bowshock model
coincides with the total mass of the cube. 
""")
# (only a small fraction of mass, {massloss:.1e} %, is lost due to numerical errors
        else:
            self._MASS_CONSISTENCY_WARNING(massloss)
        if return_isconsistent:
            return mass_consistent

    def _check_sampling(self):
        maxdz = np.max(self.km2arcsec(np.abs(self.dzs)))
        maxdr = np.max(self.km2arcsec(np.abs(self.dr)))
        maxdvs = np.max(np.abs(np.diff(self.vs)))
        maxdphi = np.max(np.abs(self.dphi))
        maxds = maxdphi * self.km2arcsec(np.max(self.rs))

        zsamp = maxdz > self.arcsecpix or maxdz > self.arcsecpix
        rsamp = maxdr > self.arcsecpix or maxdr > self.arcsecpix
        vsamp = maxdvs > self.abschanwidth
        phisamp = maxds > self.arcsecpix or maxds > self.arcsecpix

        if zsamp or rsamp:
            self._SAMPLING_XY_WARNING()
        if vsamp:
            self._SAMPLING_V_WARNING()
        if phisamp:
            self._SAMPLING_PHI_WARNING()

    def makecube(self, fromcube=None):
        """
        Makes the spectral cube of the model

        Parameters
        -----------
        fromcube : numpy.ndarray, optional
            Cube that will be populated with the model data. If None, and empty
            cube will be considered. 
        """
        if self.verbose:
            ts = []
            print("\nComputing masses in the spectral cube...")

        self.nrs = self.nzs
        self.rs = np.linspace(self.rbf, 0, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.zb_r(self.rs)
        self.dzs = self.dz_func(self.zb_r(self.rs), self.dr)

        self.phis = np.linspace(0, 2*np.pi, self.nphis+1)[:-1]
        self.dphi = self.phis[1] - self.phis[0]

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.velchans = np.linspace(self.vch0, self.vchf, self.nc)
        minvelchans = np.min(self.velchans)
        maxvelchans = np.max(self.velchans)

        if fromcube is None:
            self.cube = np.zeros((self.nc, self.nys, self.nxs))
        elif (fromcube is not None) and np.shape(fromcube)==((self.nc, self.nys, self.nxs)):
            self.cube = np.copy(fromcube)
            self._fromcube_mass = np.sum(fromcube)
        else:
            self._DIMENSION_ERROR(fromcube)

        self.cube_sampling = np.zeros((self.nc, self.nys, self.nxs))

        ci = np.cos(self.i)
        si = np.sin(self.i)
        cpa = np.cos(self.pa)
        spa = np.sin(self.pa)

        outsidegrid_warning = True
        ut.progressbar_bowshock(0, self.nzs, length=50, timelapsed=0, intervaltime=0)
        particle_in_cell = self._doCIC if self.CIC else self._doNGP
        for iz, z in enumerate(self.zs):
            if self.verbose:
                t0 = datetime.now()

            if iz == 0:
                # Treat outer boundary
                intmass = self.intmass_analytical(self.rbf)
                intmass_halfdr = self.intmass_analytical(self.rbf-self.dr/2)
                dmass =  (intmass - intmass_halfdr) / self.nphis
            elif iz == len(self.zs)-1:
                # Treat head boundary
                dmass = self.intmass_analytical(self.dr/2) / self.nphis
            else:
                # Treat the rest of the bowshock
                dmass = self.dmass_func(z, self.dzs[iz], self.dphi)

            for phi in self.phis:
            #for phi in self.phis+self.dphi*np.random.rand():
                _xp = self.rs[iz] * np.sin(phi)
                _yp = self.rs[iz] * np.cos(phi) * ci + z * si
                xp = _xp * cpa - _yp * spa
                yp = _xp * spa + _yp * cpa
                vzp = -self.vzp(z, phi)
                vlsr = vzp + self.vsys

                xpixcoord = self.km2arcsec(xp) / self.arcsecpix + self.refpix[0]
                ypixcoord = self.km2arcsec(yp) / self.arcsecpix + self.refpix[1]
                xpix = int(xpixcoord)
                ypix = int(ypixcoord)
                # Conditions model point inside cube
                condition_inside_map = \
                    (xpix+1<self.nxs) and (ypix+1<self.nys) \
                    and (xpix>0) and (ypix>0)
                condition_inside_velcoverage = \
                    vlsr <= maxvelchans and vlsr >= minvelchans
                if condition_inside_map and condition_inside_velcoverage:
                    dxpix = xpixcoord - xpix
                    dypix = ypixcoord - ypix
                    for chan, vchan in enumerate(self.velchans):
                        diffv = np.abs(vlsr-vchan)
                        if self._cond_populatechan(diffv):
                            particle_in_cell(
                                chan, diffv, xpix, ypix, dxpix, dypix, dmass)
                            if diffv < self.abschanwidth/2:
                                self._sampling(chan, xpix, ypix)
                else:
                    if outsidegrid_warning:
                        self._OUTSIDEGRID_WARNING()
                        outsidegrid_warning = False
            if self.verbose:
                tf = datetime.now()
                intervaltime = (tf-t0).total_seconds()
                ts.append(intervaltime)
                ut.progressbar_bowshock(
                    iz+1, self.nzs, np.sum(ts), intervaltime, length=50)
        self._check_mass_consistency()
        self._check_sampling()

    def makecube_variablephi(self, fromcube=None):
        """
        Makes the spectral cube of the model. This make the cube with a variable
        number of phi angle per z point. This is in principle quicker if one is
        not interested in reaching a ~0.5% of accuracy in masses.

        Parameters
        -----------
        fromcube : numpy.ndarray, optional
            Cube that will be populated with the model data. If None, and empty
            cube will be considered. 
       """
        if self.verbose:
            ts = []
            print("\nComputing masses in the spectral cube...")

        self.nrs = self.nzs
        self.rs = np.linspace(self.rbf, 0, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.zb_r(self.rs)
        self.dzs = self.dz_func(self.zb_r(self.rs), self.dr)

        self.phis = np.linspace(0, 2*np.pi, self.nphis+1)[:-1]
        self.dphi = self.phis[1] - self.phis[0]

        nphis0 = self.nphis
        phis0 = np.linspace(0, 2*np.pi, nphis0+1)[:-1]
        dphi0 = phis0[1] - phis0[0]
        ds = self.rbf * dphi0

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.velchans = np.linspace(self.vch0, self.vchf, self.nc)
        minvelchans = np.min(self.velchans)
        maxvelchans = np.max(self.velchans)

        if fromcube is None:
            self.cube = np.zeros((self.nc, self.nys, self.nxs))
        elif (fromcube is not None) and np.shape(fromcube)==((self.nc, self.nys, self.nxs)):
            self.cube = np.copy(fromcube)
            self._fromcube_mass = np.sum(fromcube)
        else:
            self._DIMENSION_ERROR(fromcube)

        self.cube_sampling = np.zeros((self.nc, self.nys, self.nxs))

        ci = np.cos(self.i)
        si = np.sin(self.i)
        cpa = np.cos(self.pa)
        spa = np.sin(self.pa)

        outsidegrid_warning = True
        ut.progressbar_bowshock(0, self.nzs, length=50, timelapsed=0, intervaltime=0)
        particle_in_cell = self._doCIC if self.CIC else self._doNGP
        for iz, z in enumerate(self.zs):
            if self.verbose:
                t0 = datetime.now()

            if iz == 0:
                # Treat outer boundary
                phis = phis0
                dphi = dphi0
                nphis = nphis0
                intmass = self.intmass_analytical(self.rbf)
                intmass_halfdr = self.intmass_analytical(self.rbf-self.dr/2)
                dmass =  (intmass - intmass_halfdr) / nphis
            elif iz == len(self.zs)-1:
                phis = phis0
                dphi = dphi0
                nphis = nphis0
                # Treat head boundary
                dmass = self.intmass_analytical(self.dr/2) / nphis
            else:
                dphi = ds / self.rs[iz]
                phis = np.arange(0, 2*np.pi, dphi)
                nphis = len(phis)
                # Treat the rest of the bowshock
                dmass = self.dmass_func(z, self.dzs[iz], dphi)

            for phi in self.phis:
            # for phi in phis+dphi*np.random.rand():
                _xp = self.rs[iz] * np.sin(phi)
                _yp = self.rs[iz] * np.cos(phi) * ci + z * si
                xp = _xp * cpa - _yp * spa
                yp = _xp * spa + _yp * cpa
                vzp = -self.vzp(z, phi)
                vlsr = vzp + self.vsys

                xpixcoord = self.km2arcsec(xp) / self.arcsecpix + self.refpix[0]
                ypixcoord = self.km2arcsec(yp) / self.arcsecpix + self.refpix[1]
                xpix = int(xpixcoord)
                ypix = int(ypixcoord)
                # Conditions model point inside cube
                condition_inside_map = \
                    (xpix+1<self.nxs) and (ypix+1<self.nys) \
                    and (xpix>0) and (ypix>0)
                condition_inside_velcoverage = \
                    vlsr <= maxvelchans and vlsr >= minvelchans
                if condition_inside_map and condition_inside_velcoverage:
                    dxpix = xpixcoord - xpix
                    dypix = ypixcoord - ypix
                    for chan, vchan in enumerate(self.velchans):
                        diffv = np.abs(vlsr-vchan)
                        if self._cond_populatechan(diffv):
                            particle_in_cell(
                                chan, diffv, xpix, ypix,
                                dxpix, dypix, dmass)
                            if diffv < self.abschanwidth/2:
                                self._sampling(chan, xpix, ypix)
                else:
                    if outsidegrid_warning:
                        self._OUTSIDEGRID_WARNING()
                        outsidegrid_warning = False
            if self.verbose:
                tf = datetime.now()
                intervaltime = (tf-t0).total_seconds()
                ts.append(intervaltime)
                ut.progressbar_bowshock(
                    iz+1, self.nzs, np.sum(ts), intervaltime, length=50)
        self._check_mass_consistency()
        self._check_sampling()

    def plot_channel(self, chan, vmax=None, vmin=None,
        cmap="inferno", savefig=None):
        """
        Plots a channel map of a cube

        Parameters
        ----------
        chan : int
            Channel map to plot
        vmax : float, optional
            Maximum value of the colormap. If None (default), the maximum value of
            the channel is chosen.
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value of
            the channel is chosen.
        cmap : str, optional
            Label of the colormap, by default "inferno".
        savefig : str, optional String of the full path to save the figure. If
            None, no figure is saved. By default, None.
        """
        fig, axs, cbax = pl.plot_channel(
            cube=self.cube,
            chan=chan,
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            vmax=vmax,
            vmin=vmin,
            cmap=cmap,
            units="Mass [Msun]",
            refpix=self.refpix,
            return_fig_axs=True
        )
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plot_channels(self, savefig=None, **kwargs):
        """
        Plots several channel map of a spectral cube.
    
        Parameters
        ----------
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
            Maximum value of the colormap. If None (default), the maximum value of
            the channel is chosen.
        vcenter : _type_, optional
            _description_, by default None
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value of
            the channel is chosen.
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
        """
        fig, ax, cbax = pl.plot_channels(
            cube=self.cube,
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            units="Mass [Msun]",
            refpix=self.refpix,
            return_fig_axs=True,
            **kwargs,
        )
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")


class CubeProcessing(BowshockCube):
    """
    Process a BowshockCube instance

    Parameters
    ----------
    bscube :  class instance
        Instance of BowshockCube
    modelname : str, optional
        Name of the folder in /models where the outputs will be saved
    J : int, optional
        Upper level of the CO rotational transition (e.g. 3 for transition
        "3-2")
    XCO : str, optional
        CO abundance relative to the molecular hydrogen
    meanmolmass : astropy.unit.Quantity, optional
        Mean mass per H molecule
    Tex : astropy.unit.Quantity, optional
        Excitation temperature
    Tbg : astropy.unit.Quantity, optional
        Excitation temperature
    coordcube : str, optional
        Set to "sky" if you would like to set the cube headers in sky
        coordinates, or "offset" if you prefer them in offsets relative to the
        origin (the source).
    ra_source_deg : float, optional
        Source right ascension [deg]
    dec_source_deg : float, optional
        Source declination [deg]
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
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
    x_FWHM : float or None
        Full width half maximum of the gaussian beam for the x direction [pixel]
    y_FWHM : float or None
        Full width half maximum of the gaussian beam for the y direction [pixel]
    beamarea : float or None
        Area of the beam [pixel^2]
    cubes : dict
        Dictionary of the processed cubes. Keys are abbreviations of the
        quantity of the cube and the operations performed to it 
    refpixs : dict
        Dictionary of the reference pixel of the cubes.  Keys are abbreviations
        of the quantity of the cube and the operations performed to it 
    hdrs : dict
        Dictionary of the headers `astropy.io.fits.header.Header` of each cube. The headers are generated when savecube method is used.
    areapix_cm : float 
        Area of a pixel in cm.
    beamarea_sr : `astropy.units.Quantity`
        Area of the beam in stereoradians.
    listmompvs : list
        List of cubes to which the moments and the position velocity diagrams are going to performed when the method self.momentsandpv_all and self.momentsandpv_and_params_all are called
    """
    default_kwargs = { }
    btypes = {
        "m": "mass",
        "I": "Intensity",
        "Ithin": "Intensity",
        "Ntot": "Total column density",
        "NCO": "CO column density",
        "tau": "Opacity"
    }
    bunits = {
        "m": "solMass",
        "I": "Jy/beam",
        "Ithin": "Jy/beam",
        "Ntot": "cm-2",
        "NCO": "cm-2",
        "tau": "-"
    }
    dos = {
        "s": "add_source",
        "r": "rotate",
        "n": "add_noise",
        "c": "convolve",
    }
    momtol_clipping = 10**(-4)
    attribs_to_get_from_cubes = [
        "arcsecpix", "nxs", "nys", "nc", "vch0", "velchans", "vchf", "xpmax",
        "distpc", "refpix", "chanwidth", "abschanwidth", "vsys"
    ]

    def __init__(
            self, modelcubes, modelname="none", J=3, XCO=8.5*10**(-5), meanmolmass=2.8, Tex=100*u.K, Tbg=2.7*u.K, coordcube="offset",
            ra_source_deg=None, dec_source_deg=None, bmin=None, bmaj=None,
            pabeam=None, papv=None, parot=None, sigma_beforeconv=None,
            maxcube2noise=None, verbose=True, **kwargs):

        if type(modelcubes) != list:
            self.nmodels = 1
            modelcubes = [modelcubes]
        if type(modelcubes) == list:
            self.nmodels = len(modelcubes)
            self.concatenate_cubes(modelcubes)

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
        self.rottrans = f"{int(J)}-{int(J-1)}"
        self.XCO = XCO
        self.meanmolmass = meanmolmass
        self.Tex = Tex
        self.Tbg = Tbg
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
        if bmin is not None and bmaj is not None:
            self.x_FWHM = self.bmin / self.arcsecpix
            self.y_FWHM = self.bmaj / self.arcsecpix
            self.beamarea = np.pi * self.y_FWHM * self.x_FWHM / (4 * np.log(2))
        else:
            self.x_FWHM = None
            self.y_FWHM = None
            self.beamarea = None

        for kwarg in self.default_kwargs:
            kwarg_attr = kwargs[kwarg] if kwarg in kwargs \
                else self.default_kwargs[kwarg]
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
        self.beamarea_sr = None
        self._calc_beamarea_sr()
        self._calc_areapix_cm()

    @staticmethod
    def _newck(ck, s):
        return f"{ck}_{s}" if "_" not in ck else ck+s

    @staticmethod
    def _q(ck):
        return ck.split("_")[0] if "_" in ck else ck

    def _getunitlabel(self, ck):
        if self.bunits[self._q(ck)] == "-":
            unitlabel = f"{self.btypes[self._q(ck)]}"
        else:
            unitlabel = f"{self.btypes[self._q(ck)]} [{self.bunits[self._q(ck)]}]"
        return unitlabel

    def _calc_beamarea_sr(self):
        self.beamarea_sr = ut.mb_sa_gaussian_f(
            self.bmaj*u.arcsec,
            self.bmin*u.arcsec
        )

    def _calc_areapix_cm(self):
        self.areapix_cm = ((self.arcsecpix * self.distpc * u.au)**2).to(u.cm**2)
    
    def _check_concat_possibility(self, modelcubes):
        for att in self.attribs_to_get_from_cubes:
            if not ut.allequal(
                [mc.__getattribute__(att) for mc in modelcubes]
                ):
                raise ValueError(
                    f"Trying to concatenate cubes with different {att}"
                    )

    def concatenate_cubes(self, modelcubes):
        self._check_concat_possibility(modelcubes)
        self.cube = np.sum(
            [modelcube.cube for modelcube in modelcubes], axis=0)

    def calc_Ntot(self,):
        """
        Computes the total (molecular hydrogen + heavier components) column
        densities of the model cube
        """
        if self.verbose:
            print(f"\nComputing column densities...")
        self.cubes["Ntot"] = rt.column_density_tot(
            m=self.cubes["m"] * u.solMass,
            meanmolmass=self.meanmolmass,
            area=self.areapix_cm,
        ).to(u.cm**(-2)).value
        self.refpixs["Ntot"] = self.refpixs["m"]
        self.noisychans["Ntot"] = self.noisychans["m"]
        self.sigma_noises["Ntot"] = self.sigma_noises["m"]
        if self.verbose:
            print(f"column densities has been calculated\n")

    def calc_NCO(self,):
        """
        Computes the CO column densities of the model cube
        """
        if self.verbose:
            print(f"\nComputing CO column densities...")
        self.cubes["NCO"] = rt.column_density_CO(
            m=self.cubes["m"] * u.solMass,
            meanmolmass=self.meanmolmass,
            area=self.areapix_cm,
            XCO=self.XCO,
        ).to(u.cm**(-2)).value
        self.refpixs["NCO"] = self.refpixs["m"]
        self.noisychans["NCO"] = self.noisychans["m"]
        self.sigma_noises["NCO"] = self.sigma_noises["m"]
        if self.verbose:
            print(f"CO column densities has been calculated\n")

    def calc_tau(self):
        """
        Computes the opacities of the model cube
        """
        if "NCO" not in self.cubes:
            self.calc_NCO()
        if self.verbose:
            print(f"\nComputing opacities...")
        self.cubes["tau"] = rt.tau_N(
            nu=rt.freq_caract_CO[self.rottrans],
            J=self.J,
            mu=0.112*u.D,
            Tex=self.Tex,
            dNdv=self.cubes["NCO"]*u.cm**(-2) / (self.abschanwidth*u.km/u.s),
        ).to("").value
        self.refpixs["tau"] = self.refpixs["m"]
        self.noisychans["tau"] = self.noisychans["m"]
        self.sigma_noises["tau"] = self.sigma_noises["m"]
        if self.verbose:
            print(f"Opacities has been calculated\n")

    def calc_I(self, opthin=False):
        """
        Calculates the intensity [Jy/beam] of the model cube.
        """
        if "tau" not in self.cubes:
            self.calc_tau()
        if self.verbose:
            print(f"\nComputing intensities...")
        func_I = rt.Inu_tau_thin if opthin else rt.Inu_tau
        ckI = "Ithin" if opthin else "I"
        self.cubes[ckI] = (func_I(
            nu=rt.freq_caract_CO[self.rottrans],
            Tex=self.Tex,
            Tbg=self.Tbg,
            tau=self.cubes["tau"],
        )*self.beamarea_sr).to(u.Jy).value
        self.refpixs[ckI] = self.refpixs["m"]
        self.noisychans[ckI] = self.noisychans["m"]
        self.sigma_noises[ckI] = self.sigma_noises["m"]
        if self.verbose:
            print(f"Intensities has been calculated\n")

    def calc_Ithin(self):
        """
        Computes the intensity [Jy/beam] of the model cube, taking into account
        the optically thin approximation 
        """
        self.calc_I(opthin=True)

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
            print(f"\nAdding source to {nck}...")
        self.cubes[nck] = np.copy(self.cubes[ck])
        value = value if value is not None else np.max(self.cubes[ck])
        if self.refpixs[ck][1]>=0 and self.refpixs[ck][0]>=0:
            self.cubes[nck][:, self.refpix[1], self.refpix[0]] = value
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = self.noisychans[ck]
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"A source has been added to {nck}, in pix [{self.refpixs[nck][0]:.2f}, {self.refpixs[nck][1]:.2f}] pix\n")

    def rotate(self, ck="m", forpv=False):
        """
        Rotates the cube an angle self.parot.

        Parameters
        -----------
        ck : str, optional
            Key of the cube to rotate 
        forpv : bool, optional
            If True, the image is rotated to calculate the PV along the bowshock
            axis 
        """
 
        nck = self._newck(ck, "r") if not forpv else self._newck(ck, "R")
        if self.verbose:
            if forpv:
                print(f"\nRotatng {nck} in order to compute the PV diagram...")
            else:
                print(f"\nRotatng {nck}...")
        # before allowing rotation of the model and not the cube
        # angle = -self.pa-90 if not forpv else self.pa+90
        # after allowing the model to be rotated
        angle = -self.parot if not forpv else self.papv + 90
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            self.cubes[nck][chan] = rotate(
                self.cubes[ck][chan],
                angle=angle,
                reshape=False,
                order=1
            )
        ang = angle * np.pi/180
        centerx = (self.nxs-1)/2
        centery = (self.nys-1)/2
        rp_center_x = self.refpixs[ck][0] - centerx
        rp_center_y = self.refpixs[ck][1] - centery
        self.refpixs[nck] = [
         +rp_center_x*np.cos(ang) + rp_center_y*np.sin(ang) + centerx,
         -rp_center_x*np.sin(ang) + rp_center_y*np.cos(ang) + centery
        ]
        self.noisychans[nck] = rotate(
            self.noisychans[ck],
            angle=angle,
            reshape=False,
            order=1,
        )
        self.sigma_noises[nck] = self.sigma_noises[ck]
        if self.verbose:
            print(f"{nck} has been rotated {angle} deg to compute the PV diagram\n")

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
            print(f"\nAdding noise to {nck}...")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])
        for chan in range(np.shape(self.cubes[ck])[0]):
            # sigma_noise = self.target_noise * 2 * np.sqrt(np.pi) \
            #          * np.sqrt(self.x_FWHM*self.y_FWHM) / 2.35
            sigma_noise = self.sigma_beforeconv if self.sigma_beforeconv is not None\
                else np.max(self.cubes[ck]) * self.maxcube2noise
            noise_matrix = np.random.normal(
                0, sigma_noise,
                size=np.shape(self.cubes[ck][chan])
                )
            self.cubes[nck][chan] = self.cubes[ck][chan] + noise_matrix
        self.refpixs[nck] = self.refpixs[ck]
        self.noisychans[nck] = noise_matrix
        self.sigma_noises[nck] = sigma_noise
        if self.verbose:
            print(f"Noise added to {nck}\n")

    def convolve(self, ck="m"):
        """
        Convolves the cube with the defined Gaussian kernel (self.bmaj,
        self.bmin, self.pabeam)
        
        Parameters
        -----------
        ck : str, optional
            Key of the cube to convolve
        """
 
        nck = self._newck(ck, "c")
        if self.verbose:
            print(f"\nConvolving {nck}... ")
        self.cubes[nck] = np.zeros_like(self.cubes[ck])

        if self.verbose:
            ts = []
            ut.progressbar_bowshock(0, self.nc,
                length=50, timelapsed=0, intervaltime=0)
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
                intervaltime = (tf-t0).total_seconds()
                ts.append(intervaltime)
                ut.progressbar_bowshock(
                    chan+1, self.nc, np.sum(ts), intervaltime, length=50)
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
{nck} has been convolved with a gaussian kernel with a size of [{self.x_FWHM:.2f}, {self.y_FWHM:.2f}] pix and with a PA of {self.pabeam:.2f}deg
"""
            )
            if "n" in nck:
                print(
f"""
The rms of the convolved image is {self.sigma_noises[nck]:.5} {self.bunits[self._q(nck)]}
""")

    def _useroutputcube2dostr(self, userdic):
        dictrad = {
            "mass": "m",
            "intensity": "I",
            "intensity_opthin": "Ithin",
            "CO_column_density": "NCO",
            "tot_column_density": "Ntot",
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
                if len(ops)>1:
                    ss = "".join([dictrad[s_user] for s_user in
                                  userdic[userkey] if s_user!="moments_and_pv"])
                    dostr = [f"{q}_{ss}"]
                if len(ops)==1:
                    dostr = [f"{q}"]
                dostrs += dostr
                self.listmompvs += dostr
            else:
                if len(ops)!=0:
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
            - "CO_column_density": Column density of the CO in cm-2.
            - "intensity": Intensity in Jy/beam.
            - "intensity_opthin": Intensity in Jy/beam, using the optically thin approximation.  
            - "tau": Opacities.
        
            The values of the dictionary are lists of strings indicating the
            operations to be performed over the cube. These are the available
            operations:
            
            - "add_source": Add a source at the reference pixel, just for spatial reference purposes.
            - "rotate": Rotate the whole spectral cube by an angle given by parot parameter.
            - "add_noise": Add gaussian noise, defined by maxcube2noise parameter.
            - "convolve": Convolve with a gaussian defined by the parameters bmaj, bmin, and pabeam.
            - "moments_and_pv": Computes the moments 0, 1, and 2, the maximum intensity and the PV diagram.
        
            The operations will be performed folowing the order of the strings
            in the list (from left to right). The list can be left empty if no
            operations are desired.
            
        Example:
        --------
        >>> cp = bs.CubeProcessing(...)
        >>> outcubes = {
        >>>    "intensity": ["add_noise", "convolve", "moments_and_pv"],
        >>>    "opacity": [],
        >>>    "CO_column_density": ["convolve"],
        >>>    "mass": [],
        >>> }
        >>> cp.calc(outcubes)

        will save 4 spectral cubes in fits format. The first one are the
        intensities with gaussian noise added, it will be convolved, and the
        moments and PV diagrams will be computed; the second cube will be the
        opacity; the third will be the CO_column_density, which will be
        convolved; and the forth cube will be the masses. The first spectral
        cube will be named I_nc.fits, the second tau.fits, the third NCO_c.fits,
        and the fourth m.fits.  
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
                    ck = q if i==0 else f"{q}_{ss[:i]}"
                    if self._newck(ck, s) not in self.cubes:
                        self.__getattribute__(self.dos[s])(ck=ck)
    
    def savecube(self, ck):
        """
        Saves the cube in fits format

        Parameters
        -----------
        ck : str
            Key of the cube to convolve
        """
        if self.coordcube == "offset":
            ctype1 = 'OFFSET'
            ctype2 = 'OFFSET'
            cunit1 = 'arcsec'
            cunit2 = 'arcsec'
            crval1 = 0
            crval2 = 0
            cdelt1 = self.arcsecpix
            cdelt2 = self.arcsecpix
        else:
            ctype1 = 'RA---SIN'
            ctype2 = 'DEC--SIN'
            cunit1 = 'deg'
            cunit2 = 'deg'
            crval1 = self.ra_source_deg
            crval2 = self.dec_source_deg
            cdelt1 = -self.arcsecpix / 3600
            cdelt2 = self.arcsecpix / 3600
        hdr = fits.Header()
        hdr["SIMPLE"] = True
        hdr["BITPIX"] = -32
        hdr["NAXIS"] = 3
        hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
        hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
        hdr["NAXIS3"] = np.shape(self.cubes[ck])[0]
        hdr["EXTEND"] = True
        hdr["BSCALE"] = 1
        hdr["BZERO"] = 0
        hdr["BMAJ"] = self.bmaj / 3600
        hdr["BMIN"] = self.bmin / 3600
        hdr["BPA"] = self.pabeam
        hdr["BTYPE"] = self.btypes[self._q(ck)]
        hdr["OBJECT"] = f'{self.modelname}'
        hdr["BUNIT"] = self.bunits[self._q(ck)]
        hdr["RADESYS"] = 'ICRS'
        hdr["LONPOLE"] = 1.800000000000E+02
        hdr["LATPOLE"] = 3.126777777778E+01
        hdr["PC1_1"] = 1
        hdr["PC2_1"] = 0
        hdr["PC1_2"] = 0
        hdr["PC2_2"] = 1
        hdr["PC1_3"] = 0
        hdr["PC2_3"] = 0
        hdr["PC3_1"] = 0
        hdr["PC3_2"] = 0
        hdr["PC3_3"] = 1
        hdr["CTYPE1"] = ctype1
        hdr["CRVAL1"] = crval1
        hdr["CDELT1"] = cdelt1
        hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
        hdr["CUNIT1"] = cunit1
        hdr["CTYPE2"] = ctype2
        hdr["CRVAL2"] = crval2
        hdr["CDELT2"] = cdelt2
        hdr["CRPIX2"] = self.refpixs[ck][1] + 1
        hdr["CUNIT2"] = cunit2
        hdr["CTYPE3"] = 'VRAD'
        hdr["CRVAL3"] = self.velchans[0]
        hdr["CDELT3"] = self.velchans[1] - self.velchans[0]
        hdr["CRPIX3"] = 1
        hdr["CUNIT3"] = 'km/s'
        hdr["PV2_1"] = 0
        hdr["PV2_2"] = 0
        hdr["RESTFRQ"] = 3.457380000000E+11
        hdr["SPECSYS"] = 'LSRK'
        hdr["ALTRVAL"] = 7.757120529450E+02
        hdr["ALTRPIX"] = -2.700000000000E+01
        hdr["VELREF"] = 257
        hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
        hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
        hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
        hdr["TIMESYS"] = 'UTC'
        hdr["OBSRA"] = 5.226562499999E+01
        hdr["OBSDEC"] = 3.126777777778E+01
        hdr["OBSGEO-X"] = 2.225142180269E+06
        hdr["OBSGEO-Y"] = -5.440307370349E+06
        hdr["OBSGEO-Z"] = -2.481029851874E+06
        hdr["DATE"] = f'{datetime.now().isoformat()}'
        hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

        self.hdrs[ck] = hdr
        hdu = fits.PrimaryHDU(self.cubes[ck])
        hdul = fits.HDUList([hdu])
        hdu.header = self.hdrs[ck]
        ut.make_folder(foldername=f'models/{self.modelname}/fits')
        hdul.writeto(f'models/{self.modelname}/fits/{ck}.fits', overwrite=True)
        if self.verbose:
            print(f'models/{self.modelname}/fits/{ck}.fits saved')

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
        >>> bscp = bs.CubeProcessing(bsc, ...)
        >>> outcubes = {
        >>>    "intensity": ["add_noise", "convolve", "moments_and_pv"],
        >>>    "opacity": [],
        >>>    "CO_column_density": ["convolve"],
        >>>    "mass": [],
        >>> }
        >>> bscp.savecubes(outputcubes)

        will save 4 spectral cubes in fits format. The first one are the
        intensities with gaussian noise added, it will be convolved, and the
        moments and PV diagrams will be computed; the second cube will be the
        opacity; the third will be the CO_column_density, which will be
        convolved; and the forth cube will be the masses. The first spectral
        cube will be named I_nc.fits, the second tau.fits, the third NCO_c.fits,
        and the fourth m.fits.  
        """
        cks = self._useroutputcube2dostr(userdic)
        for ck in cks:
            self.savecube(ck)

    def plot_channel(self, ck, chan, vmax=None, vmin=None,
                    cmap="inferno", savefig=None, add_beam=False, return_fig_axs=False):
        """
        Plots a channel map of a cube

        Parameters
        ----------
        ck : str
            Key of the cube to plot (see keys of self.cubes dictionary)
        chan : int
            Channel map to plot
        vmax : float, optional
            Maximum value of the colormap. If None (default), the maximum value of
            the channel is chosen.
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value of
            the channel is chosen.
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

    def plot_channels(self, ck, savefig=None, add_beam=False, 
                      return_fig_axs=False, **kwargs):
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
            Maximum value of the colormap. If None (default), the maximum value of
            the channel is chosen.
        vcenter : _type_, optional
            _description_, by default None
        vmin : float, optional
            Minimum value of the colormap. If None (default), the minimum value of
            the channel is chosen.
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
            If True, returns a tuple of the axes of the channel map and the
            colorbar. If False, does not return anything.
        
        Returns
        --------
        (fig, ax, cbax) : tuple of matplotlib.axes.Axes Axes of the channel map
            and the colorbar, only returns if return_fig_axs=True.
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

    def pvalongz(self, ck, halfwidth=0, savefits=False, filename=None):
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
        filename : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        pvimage : numpy.ndarray
            Position velocity diagram
        """
        pvimage = moments.pv(
            self.cubes[ck],
            int(self.refpixs[ck][1]),
            halfwidth=halfwidth,
            axis=1
        )
        if savefits:
            hdrpv = fits.Header()
            hdrpv["SIMPLE"] = True
            hdrpv["BITPIX"] = -32
            hdrpv["NAXIS"] = 2
            hdrpv["NAXIS1"] = np.shape(self.cubes[ck])[1]
            hdrpv["NAXIS2"] = np.shape(self.cubes[ck])[0]
            hdrpv["EXTEND"] = True
            hdrpv["BSCALE"] = 1
            hdrpv["BZERO"] = 0
            hdrpv["BTYPE"] = self.btypes[self._q(ck)]
            hdrpv["OBJECT"] = f'{self.modelname}'
            hdrpv["BUNIT"] = self.bunits[self._q(ck)]
            hdrpv["RADESYS"] = 'ICRS'
            hdrpv["LONPOLE"] = 1.800000000000E+02
            hdrpv["LATPOLE"] = 3.126777777778E+01
            hdrpv["PC1_1"] = 1
            hdrpv["PC2_1"] = 0
            hdrpv["PC1_2"] = 0
            hdrpv["PC2_2"] = 1
            hdrpv["CTYPE1"] = "OFFSET"
            hdrpv["CRVAL1"] = 0
            hdrpv["CDELT1"] = self.arcsecpix
            hdrpv["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdrpv["CUNIT1"] = "arcsec"
            hdrpv["CTYPE2"] = "VELOCITY"
            hdrpv["CRVAL2"] = self.velchans[0]
            hdrpv["CDELT2"] = self.chanwidth
            hdrpv["CRPIX2"] = 1
            hdrpv["CUNIT2"] = "km/s"
            hdrpv["PV2_1"] = 0
            hdrpv["PV2_2"] = 0
            hdrpv["RESTFRQ"] = 3.457380000000E+11
            hdrpv["SPECSYS"] = 'LSRK'
            hdrpv["ALTRVAL"] = 7.757120529450E+02
            hdrpv["ALTRPIX"] = -2.700000000000E+01
            hdrpv["VELREF"] = 257
            hdrpv["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdrpv["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdrpv["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdrpv["TIMESYS"] = 'UTC'
            hdrpv["OBSRA"] = 5.226562499999E+01
            hdrpv["OBSDEC"] = 3.126777777778E+01
            hdrpv["OBSGEO-X"] = 2.225142180269E+06
            hdrpv["OBSGEO-Y"] = -5.440307370349E+06
            hdrpv["OBSGEO-Z"] = -2.481029851874E+06
            hdrpv["DATE"] = f'{datetime.now().isoformat()}'
            hdrpv["ORIGIN"] = f'BOWSHOCKPY v{__version__}'
           
            hdu = fits.PrimaryHDU(pvimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdrpv
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_pv.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_pv.fits saved')
        return pvimage

    def sumint(self, ck, chan_range=None, savefits=False, filename=None):
        """
        Computes the image of the summation of pixels of the cube along the
        velocity axis
        
        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the summation. If None, the whole cube will be considered.
        savefits : boolean
            If True, save the PV-diagram in fits format.
        filename : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        sumint : numpy.ndarray 
            Image of the summation of the pixels of the cube along the velocty
            axis
 
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        sumintimage = moments.sumint(
            self.cubes[ck],
            chan_range=chan_range
        )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = 'OFFSET'
                ctype2 = 'OFFSET'
                cunit1 = 'arcsec'
                cunit2 = 'arcsec'
                crval1 = 0
                crval2 = 0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = 'RA---SIN'
                ctype2 = 'DEC--SIN'
                cunit1 = 'deg'
                cunit2 = 'deg'
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = fits.Header()
            hdr["SIMPLE"] = True
            hdr["BITPIX"] = -32
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            hdr["BSCALE"] = 1
            hdr["BZERO"] = 0
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f'{self.modelname}'
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["RADESYS"] = 'ICRS'
            hdr["LONPOLE"] = 1.800000000000E+02
            hdr["LATPOLE"] = 3.126777777778E+01
            hdr["PC1_1"] = 1
            hdr["PC2_1"] = 0
            hdr["PC1_2"] = 0
            hdr["PC2_2"] = 1
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1
            hdr["CUNIT2"] = cunit2
            hdr["PV2_1"] = 0
            hdr["PV2_2"] = 0
            hdr["RESTFRQ"] = 3.457380000000E+11
            hdr["SPECSYS"] = 'LSRK'
            hdr["ALTRVAL"] = 7.757120529450E+02
            hdr["ALTRPIX"] = -2.700000000000E+01
            hdr["VELREF"] = 257
            hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdr["TIMESYS"] = 'UTC'
            hdr["OBSRA"] = 5.226562499999E+01
            hdr["OBSDEC"] = 3.126777777778E+01
            hdr["OBSGEO-X"] = 2.225142180269E+06
            hdr["OBSGEO-Y"] = -5.440307370349E+06
            hdr["OBSGEO-Z"] = -2.481029851874E+06
            hdr["DATE"] = f'{datetime.now().isoformat()}'
            hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

            hdu = fits.PrimaryHDU(sumintimage)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_sumint.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_sumint.fits saved')
        return sumintimage

    def mom0(self, ck, chan_range=None, savefits=False, filename=None):
        """
        Computes the 0th order moment along the velocity axis
        
        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the summation. If None, the whole cube will be considered.
        savefits : boolean
            If True, save the PV-diagram in fits format.
        filename : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom0 : numpy.ndarray 
            Moment 0 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        mom0 = moments.mom0(
            self.cubes[ck],
            chan_vels=chan_vels,
            chan_range=chan_range,
            )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = 'OFFSET'
                ctype2 = 'OFFSET'
                cunit1 = 'arcsec'
                cunit2 = 'arcsec'
                crval1 = 0
                crval2 = 0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = 'RA---SIN'
                ctype2 = 'DEC--SIN'
                cunit1 = 'deg'
                cunit2 = 'deg'
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = fits.Header()
            hdr["SIMPLE"] = True
            hdr["BITPIX"] = -32
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            hdr["BSCALE"] = 1
            hdr["BZERO"] = 0
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f'{self.modelname}'
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["RADESYS"] = 'ICRS'
            hdr["LONPOLE"] = 1.800000000000E+02
            hdr["LATPOLE"] = 3.126777777778E+01
            hdr["PC1_1"] = 1
            hdr["PC2_1"] = 0
            hdr["PC1_2"] = 0
            hdr["PC2_2"] = 1
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1
            hdr["CUNIT2"] = cunit2
            hdr["PV2_1"] = 0
            hdr["PV2_2"] = 0
            hdr["RESTFRQ"] = 3.457380000000E+11
            hdr["SPECSYS"] = 'LSRK'
            hdr["ALTRVAL"] = 7.757120529450E+02
            hdr["ALTRPIX"] = -2.700000000000E+01
            hdr["VELREF"] = 257
            hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdr["TIMESYS"] = 'UTC'
            hdr["OBSRA"] = 5.226562499999E+01
            hdr["OBSDEC"] = 3.126777777778E+01
            hdr["OBSGEO-X"] = 2.225142180269E+06
            hdr["OBSGEO-Y"] = -5.440307370349E+06
            hdr["OBSGEO-Z"] = -2.481029851874E+06
            hdr["DATE"] = f'{datetime.now().isoformat()}'
            hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

            hdu = fits.PrimaryHDU(mom0)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_mom0.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom0.fits saved')
        return mom0

    def mom1(self, ck, chan_range=None, clipping=0, savefits=False, filename=None):
        """
        Computes the 1th order moment along the velocity axis
        
        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the summation. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter will be masked with 0 values.
        savefits : boolean, optional
            If True, save the PV-diagram in fits format.
        filename : str
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom1 : numpy.ndarray 
            Moment 1 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom1 = np.nan_to_num(
                moments.mom1(
                    cube_clipped,
                    chan_vels=chan_vels,
                    chan_range=chan_range
                )
            )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = 'OFFSET'
                ctype2 = 'OFFSET'
                cunit1 = 'arcsec'
                cunit2 = 'arcsec'
                crval1 = 0
                crval2 = 0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = 'RA---SIN'
                ctype2 = 'DEC--SIN'
                cunit1 = 'deg'
                cunit2 = 'deg'
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = fits.Header()
            hdr["SIMPLE"] = True
            hdr["BITPIX"] = -32
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            hdr["BSCALE"] = 1
            hdr["BZERO"] = 0
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f'{self.modelname}'
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["RADESYS"] = 'ICRS'
            hdr["LONPOLE"] = 1.800000000000E+02
            hdr["LATPOLE"] = 3.126777777778E+01
            hdr["PC1_1"] = 1
            hdr["PC2_1"] = 0
            hdr["PC1_2"] = 0
            hdr["PC2_2"] = 1
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1
            hdr["CUNIT2"] = cunit2
            hdr["PV2_1"] = 0
            hdr["PV2_2"] = 0
            hdr["RESTFRQ"] = 3.457380000000E+11
            hdr["SPECSYS"] = 'LSRK'
            hdr["ALTRVAL"] = 7.757120529450E+02
            hdr["ALTRPIX"] = -2.700000000000E+01
            hdr["VELREF"] = 257
            hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdr["TIMESYS"] = 'UTC'
            hdr["OBSRA"] = 5.226562499999E+01
            hdr["OBSDEC"] = 3.126777777778E+01
            hdr["OBSGEO-X"] = 2.225142180269E+06
            hdr["OBSGEO-Y"] = -5.440307370349E+06
            hdr["OBSGEO-Z"] = -2.481029851874E+06
            hdr["DATE"] = f'{datetime.now().isoformat()}'
            hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

            hdu = fits.PrimaryHDU(mom1)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_mom1.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom1.fits saved')
        return mom1

    def mom2(self, ck, chan_range=None, clipping=0, savefits=False, filename=None):
        """
        Computes the 2th order moment along the velocity axis
        
        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the summation. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter will be masked with 0 values.
        savefits : boolean, optional
            If True, save the PV-diagram in fits format.
        filename : str, optional
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom2 : numpy.ndarray 
            Moment 2 image of the cube
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom2 = np.nan_to_num(
                moments.mom2(
                    cube_clipped,
                    chan_vels=chan_vels,
                    chan_range=chan_range
                )
            )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = 'OFFSET'
                ctype2 = 'OFFSET'
                cunit1 = 'arcsec'
                cunit2 = 'arcsec'
                crval1 = 0
                crval2 = 0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = 'RA---SIN'
                ctype2 = 'DEC--SIN'
                cunit1 = 'deg'
                cunit2 = 'deg'
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = fits.Header()
            hdr["SIMPLE"] = True
            hdr["BITPIX"] = -32
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            hdr["BSCALE"] = 1
            hdr["BZERO"] = 0
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f'{self.modelname}'
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["RADESYS"] = 'ICRS'
            hdr["LONPOLE"] = 1.800000000000E+02
            hdr["LATPOLE"] = 3.126777777778E+01
            hdr["PC1_1"] = 1
            hdr["PC2_1"] = 0
            hdr["PC1_2"] = 0
            hdr["PC2_2"] = 1
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1
            hdr["CUNIT2"] = cunit2
            hdr["PV2_1"] = 0
            hdr["PV2_2"] = 0
            hdr["RESTFRQ"] = 3.457380000000E+11
            hdr["SPECSYS"] = 'LSRK'
            hdr["ALTRVAL"] = 7.757120529450E+02
            hdr["ALTRPIX"] = -2.700000000000E+01
            hdr["VELREF"] = 257
            hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdr["TIMESYS"] = 'UTC'
            hdr["OBSRA"] = 5.226562499999E+01
            hdr["OBSDEC"] = 3.126777777778E+01
            hdr["OBSGEO-X"] = 2.225142180269E+06
            hdr["OBSGEO-Y"] = -5.440307370349E+06
            hdr["OBSGEO-Z"] = -2.481029851874E+06
            hdr["DATE"] = f'{datetime.now().isoformat()}'
            hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

            hdu = fits.PrimaryHDU(mom2)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_mom1.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom2.fits saved')
        return mom2

    def mom8(self, ck, chan_range=None, clipping=0, savefits=False, filename=None):
        """
        Computes the maximum value of the cube along the velocity axis
        
        Parameters
        -----------
        ck : str
            Key of the cube to perfomr the PV-diagram.
        chan_range : list, optional
            Two element list with the last and first channels used to compute
            the summation. If None, the whole cube will be considered.
        clipping : float, optional
            Pixels with values smaller than the one given by clipping parameter
            will be masked with 0 values.
        savefits : boolean, optional
            If True, save the PV-diagram in fits format.
        filename : str, optional
            Full path name of the fits file. If None, it will be saved as
            models/{self.modelname}/fits/{ck}_pv.fits. If the path does not
            exist, it will be created.

        Returns
        --------
        mom8 : numpy.ndarray 
            Maximum value of the pixels of the cubes along the velocity axis 

        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        # chan_vels = self.velchans[chan_range[0]:chan_range[-1]]
        cube_clipped = np.copy(self.cubes[ck])
        clipping = clipping if clipping != 0 \
            else self.momtol_clipping * np.max(self.cubes[ck])
        cube_clipped[cube_clipped<clipping] = 0
        mom8 = np.nan_to_num(
                moments.mom8(
                    cube_clipped,
                    chan_range=chan_range,
                )
            )
        if savefits:
            if self.coordcube == "offset":
                ctype1 = 'OFFSET'
                ctype2 = 'OFFSET'
                cunit1 = 'arcsec'
                cunit2 = 'arcsec'
                crval1 = 0
                crval2 = 0
                cdelt1 = self.arcsecpix
                cdelt2 = self.arcsecpix
            else:
                ctype1 = 'RA---SIN'
                ctype2 = 'DEC--SIN'
                cunit1 = 'deg'
                cunit2 = 'deg'
                crval1 = self.ra_source_deg
                crval2 = self.dec_source_deg
                cdelt1 = -self.arcsecpix / 3600
                cdelt2 = self.arcsecpix / 3600
            hdr = fits.Header()
            hdr["SIMPLE"] = True
            hdr["BITPIX"] = -32
            hdr["NAXIS"] = 2
            hdr["NAXIS1"] = np.shape(self.cubes[ck])[2]
            hdr["NAXIS2"] = np.shape(self.cubes[ck])[1]
            hdr["EXTEND"] = True
            hdr["BSCALE"] = 1
            hdr["BZERO"] = 0
            hdr["BMAJ"] = self.bmaj / 3600
            hdr["BMIN"] = self.bmin / 3600
            hdr["BPA"] = self.pabeam
            hdr["BTYPE"] = self.btypes[self._q(ck)]
            hdr["OBJECT"] = f'{self.modelname}'
            hdr["BUNIT"] = self.bunits[self._q(ck)]
            hdr["RADESYS"] = 'ICRS'
            hdr["LONPOLE"] = 1.800000000000E+02
            hdr["LATPOLE"] = 3.126777777778E+01
            hdr["PC1_1"] = 1
            hdr["PC2_1"] = 0
            hdr["PC1_2"] = 0
            hdr["PC2_2"] = 1
            hdr["CTYPE1"] = ctype1
            hdr["CRVAL1"] = crval1
            hdr["CDELT1"] = cdelt1
            hdr["CRPIX1"] = self.refpixs[ck][0] + 1.
            hdr["CUNIT1"] = cunit1
            hdr["CTYPE2"] = ctype2
            hdr["CRVAL2"] = crval2
            hdr["CDELT2"] = cdelt2
            hdr["CRPIX2"] = self.refpixs[ck][1] + 1
            hdr["CUNIT2"] = cunit2
            hdr["PV2_1"] = 0
            hdr["PV2_2"] = 0
            hdr["RESTFRQ"] = 3.457380000000E+11
            hdr["SPECSYS"] = 'LSRK'
            hdr["ALTRVAL"] = 7.757120529450E+02
            hdr["ALTRPIX"] = -2.700000000000E+01
            hdr["VELREF"] = 257
            hdr["TELESCOP"] = f'BOWSHOCKPY v{__version__}'
            hdr["OBSERVER"] = f'BOWSHOCKPY v{__version__}'
            hdr["DATE-OBS"] = f'{datetime.now().isoformat()}'
            hdr["TIMESYS"] = 'UTC'
            hdr["OBSRA"] = 5.226562499999E+01
            hdr["OBSDEC"] = 3.126777777778E+01
            hdr["OBSGEO-X"] = 2.225142180269E+06
            hdr["OBSGEO-Y"] = -5.440307370349E+06
            hdr["OBSGEO-Z"] = -2.481029851874E+06
            hdr["DATE"] = f'{datetime.now().isoformat()}'
            hdr["ORIGIN"] = f'BOWSHOCKPY v{__version__}'

            hdu = fits.PrimaryHDU(mom8)
            hdul = fits.HDUList([hdu])
            hdu.header = hdr
            if filename is None:
                ut.make_folder(foldername=f'models/{self.modelname}/fits')
                filename = f'models/{self.modelname}/fits/{ck}_mom1.fits' 
            hdul.writeto(filename, overwrite=True)
            if self.verbose:
                print(f'models/{self.modelname}/fits/{ck}_mom8.fits saved')
        return mom8

    def plotpv(
            self, ck, halfwidth, ax=None, cbax=None,
            savefits=False, return_im=False, **kwargs
            ):
        """
        Plots the position velocity diagram.

        Parameters
        -----------
        pvimage : numpy.ndarray
            PV-diagram.
        chan_vels: list or numpy.ndarray
            1-dimensional array of the velocity corresponding to the channels.
        """
        ckpv = ck + "R"
        if ckpv not in self.cubes:
            self.rotate(ck, forpv=True)

        pvimage = self.pvalongz(
            ckpv,
            halfwidth=halfwidth,
            savefits=savefits,
            )
        rangex = np.array([
            -0.5-self.refpixs[ckpv][0],
            self.nxs-0.5-self.refpixs[ckpv][0]
            ]) * self.arcsecpix
        pl.plotpv(
            pvimage,
            rangex=rangex,
            chan_vels=self.velchans,
            ax=ax,
            cbax=cbax,
            cbarlabel=self._getunitlabel(ckpv),
            **kwargs,
            )
        if return_im:
            return pvimage

    def plotsumint(
            self, ck, chan_range=None, ax=None, cbax=None,
            savefits=False, return_im=False, add_beam=False,
              **kwargs):
        """
        Plots the sum of the pixels of the cubes along the velocity axis.

        Parameters
        -----------
        sumint : numpy.ndarray
            Sum of the pixels of the cubes along the velocity axis.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        sumint = self.sumint(
            ck, chan_range=chan_range, savefits=savefits,
            )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        pl.plotsumint(
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
            **kwargs,
            )
        if return_im:
            return sumint

    def plotmom0(
            self, ck, chan_range=None, ax=None, cbax=None,
            savefits=False, return_im=False, add_beam=False, 
            **kwargs):
        """
        Plots the moment 0.

        Parameters
        -----------
        mom0 : numpy.ndarray
            Moment 0 image.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        mom0 = self.mom0(
            ck, chan_range=chan_range, savefits=savefits,
            )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        pl.plotmom0(
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
            **kwargs,
            )
        if return_im:
            return mom0

    def plotmom1(
            self, ck, chan_range=None, mom1clipping=0, 
            ax=None, cbax=None, savefits=False, 
            return_im=False, add_beam=False, **kwargs
            ):
        """
        Plots the moment 1.

        Parameters
        -----------
        mom1 : numpy.ndarray
            Moment 1 image.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        clipping = float(mom1clipping.split("x")[0]) * self.sigma_noises[ck] \
            if mom1clipping !=0 else 0
        add_beam = add_beam if "c" in ck else False
        mom1 = self.mom1(
                ck,
                chan_range=chan_range,
                savefits=savefits,
                clipping=clipping,
            )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        pl.plotmom1(
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
            **kwargs,
            )
        if return_im:
            return mom1

    def plotmom2(
            self, ck, chan_range=None, mom2clipping=0, 
            ax=None, cbax=None, add_beam=False, savefits=False,
            return_im=False, **kwargs
            ):
        """
        Plots the moment 2.

        Parameters
        -----------
        mom2 : numpy.ndarray
            Moment 2 image.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        clipping = float(mom2clipping.split("x")[0]) * self.sigma_noises[ck]\
            if mom2clipping !=0 else 0
        add_beam = add_beam if "c" in ck else False
        mom2 = self.mom2(
                    ck,
                    chan_range=chan_range,
                    savefits=savefits,
                    clipping=clipping,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        pl.plotmom2(
            mom2,
            extent=extent,
            interpolation=None,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Velocity dispersion [km$^2$/s$^2$]",
            **kwargs,
            )
        if return_im:
            return mom2

    def plotmom8(
            self, ck, chan_range=None, ax=None, cbax=None,
            add_beam=False, savefits=False, return_im=False,
             **kwargs):
        """
        Plots the moment 8.

        Parameters
        -----------
        mom8 : numpy.ndarray
            Moment 8 image.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        """
        chan_range = chan_range if chan_range is not None else [0, self.nc]
        add_beam = add_beam if "c" in ck else False
        mom8 = self.mom8(
                    ck,
                    chan_range=chan_range,
                    savefits=savefits,
                )
        extent = np.array([
            -(-0.5-self.refpixs[ck][0]),
            -(self.nxs-0.5-self.refpixs[ck][0]),
            (-0.5-self.refpixs[ck][1]),
            (self.nys-0.5-self.refpixs[ck][1]),
            ]) * self.arcsecpix
        pl.plotmom8(
            mom8,
            extent=extent,
            ax=ax,
            cbax=cbax,
            add_beam=add_beam,
            bmin=self.bmin,
            bmaj=self.bmaj,
            pabeam=self.pabeam,
            cbarlabel="Peak " + self._getunitlabel(ck),
            **kwargs,
            )
        if return_im:
            return mom8

    def momentsandpv_and_params(
            self, ck, savefits=False, saveplot=False,
            mom1clipping=0, mom2clipping=0, verbose=True,
            chan_range=None, halfwidth_pv=0, add_beam=False,
            mom0values={v: None for v in ["vmax", "vcenter", "vmin"]},
            mom1values={v: None for v in ["vmax", "vcenter", "vmin"]},
            mom2values={v: None for v in ["vmax", "vcenter", "vmin"]},
            mom8values={v: None for v in ["vmax", "vcenter", "vmin"]},
            pvvalues={v: None for v in ["vmax", "vcenter", "vmin"]},
            ):
        """
        Computes the moments and position velocity diagram including also the
        main parameters of the model listed in the first ax
        
        Parameters
        -----------
        ck : str
            Key of the cube to which the moments and the position velocity diagram will be computed.
        savefits : bool
            If True, the moments and the position velocity diagram will be saved in fits format
        saveplot : bool
            If True, a plot of the moments and position velocity diagrams will be saved
        mom1clipping : float
            Clipping to in order to compute the moment 1. Pixels with values
            smaller than the one given by clipping parameter will be masked with
            0 values.
        add_beam : bolean, optional
            If True, plots a ellipse of the beam size in the left bottom corner.
        mom2clipping : float
            Clipping to in order to compute the moment 2. Pixels with values
            smaller than the one given by clipping parameter will be masked with
            0 values.
        mom8clipping : float
            Clipping to in order to compute the maximum value of the pixels along the velocity axis. Pixels with values smaller than the one given by clipping parameter will be masked with 0 values.
        mom0values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 0. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value of
            the moment image will be considered, respectively. Example:
            mom0values = {"vmax": None, "vcenter": None, "vmin": 0,}. 
        mom1values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 1. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value of
            the moment image will be considered, respectively. Example:
            mom1values = {"vmax": None, "vcenter": None, "vmin": 0,}. 
        mom2values : dict or None
            Dictionary with the maximum, central, and minimum value to show in
            the plot of the moment 2. If the dictionary value is None for vmax,
            vcenter, or vmin, then the maximum, central, or the minimum value of
            the moment image will be considered, respectively. Example:
            mom1values = {"vmax": None, "vcenter": None, "vmin": 0,}. 
        mom8values : dict or None
           Dictionary with the maximum, central, and minimum value to show in
           the plot of the maximum value along the velocity axis. If the
           dictionary value is None for vmax, vcenter, or vmin, then the
           maximum, central, or the minimum value of the moment image will be
           considered, respectively. Example: mom8values = {"vmax": None,
           "vcenter": None, "vmin": None,}. 
        """
 
        if verbose:
            print(
"""
\nComputing moments and the PV-diagram along the jet axis
"""
            )

        # ckpv = ck + "R"
        # if ckpv not in self.cubes:
        #     self.rotate(ck, forpv=True)

        fig = plt.figure(figsize=(14,10))
        gs = GridSpec(
            2, 3,
            width_ratios=[1]*3, # + [0.85]*2,
            height_ratios=[1]*2,
            hspace=0.3,
            wspace=0.25,
        )
        axs = {}
        cbaxs = {}
        gss = {}

        ik = "text"
        axs[ik] = plt.subplot(gs[0,0])
        axs[ik].set_axis_off()

        ik = "mom0"
        gss[ik] = gs[0,1].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "mom8"
        gss[ik] = gs[0,2].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "pv"
        gss[ik] = gs[1,0].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])


        ik = "mom1"
        gss[ik] = gs[1,1].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ik = "mom2"
        gss[ik] = gs[1,2].subgridspec(
                 2, 1,
                 height_ratios=[0.05, 1],
                 width_ratios=[1],
                 hspace=0.05,
             )
        axs[ik] = plt.subplot(gss[ik][1,0])
        cbaxs[ik] = plt.subplot(gss[ik][0,0])

        ak = "text"

        showtext = \
        fr"""
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
            axs["text"].text(0, 0.99-0.06*n, line, fontsize=12-self.nmodels,
                              transform=axs["text"].transAxes)

        add_beam = add_beam if "c" in ck else False
        ak = "mom0"
        self.plotmom0(
            ck=ck,
            chan_range=chan_range,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_im=False,
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
            return_im=False,
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
            return_im=False,
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
            return_im=False,
            **pvvalues,
            )

        ak = "mom8"
        self.plotmom8(
            ck=ck,
            chan_range=chan_range,
            ax=axs[ak],
            cbax=cbaxs[ak],
            savefits=savefits,
            return_im=False,
            add_beam=add_beam,
            **mom8values,
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



