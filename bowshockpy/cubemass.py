"""This module contains the tools that allows to compute the masses in a
spectral cube of a bowshock model"""

import warnings
from datetime import datetime

import numpy as np

import bowshockpy.plots as pl
import bowshockpy.utils as ut
from bowshockpy.modelproj import ObsModel
from bowshockpy.version import __version__

warnings.formatwarning = ut.formatwarning
warnings.filterwarnings("error", category=ut.UserError)
warnings.filterwarnings("always", category=UserWarning)


class BowshockCube(ObsModel):
    """
    Computes the spectral cube of the bowshock model

    Parameters
    -----------
    obsmodel : class instance
        Instance of ObsModel
    nphis : int
        Number of azimuthal angles phi to calculate the bowshock solution
    xpmax : float
        Physical size of the channel maps along the x axis [arcsec]
    vch0 : float
        Central velocity of the first channel map [km/s]
    vchf : float | None
        Central velocity of the last channel map [km/s]. If None, you should
        provide instead the channel width with chanwidth parameter. If a float
        is provided, chanwidth should be None, since chanwidth would be
        computed internally. Default is None.
    chanwidth : float | None
        Channel width of the spectral cube [km/s]. If positive, vch0<vchf,
        negative, vch0>vchf. If None, you should provide instead the central
        velocity of the last channel map. If a float is provided, vchf should
        be None, since vchf would be computed internally. Default is None.
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
    cic : bolean, optional
        Set to True to perform Cloud in Cell interpolation [1].
    vt : float | str, optional
        Thermal+turbulent line-of-sight velocity dispersion [km/s] If
        thermal+turbulent line-of-sight velocity dispersion is smaller than the
        instrumental spectral resolution, vt should be the spectral resolution.
        It can be also set to a integer times the channel width (e.g.,
        "2xchannel")
    tolfactor_vt : float, optional
        Neighbour channel maps around a given channel map with vch will stop
        being populated when their difference in velocity with respect to vch
        is higher than this factor times vt. The lower the factor, the quicker
        will be the code, but the total mass will be underestimated. If vt is
        not None, compare the total mass of the output cube with the mass
        parameter that the user has defined
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
    Many-Particle Physics, Vol. 739 (Springer), doi: 10.1007/978-3-540-74686-7.

    """

    def __init__(
        self,
        obsmodel,
        nphis,
        xpmax,
        vch0,
        vchf=None,
        chanwidth=None,
        nzs=200,
        nc=50,
        nxs=200,
        nys=200,
        refpix=[0, 0],
        cic=True,
        vt="2xchannel",
        tolfactor_vt=None,
        verbose=True,
        **kwargs,
    ):
        self.__dict__ = obsmodel.__dict__
        self.nphis = nphis
        self.xpmax = xpmax
        self.vch0 = vch0
        self.vchf = vchf
        self.chanwidth = chanwidth
        self.nzs = nzs
        self.nc = nc
        self.nxs = nxs
        self.nys = nys
        self.refpix = refpix
        self.cic = cic
        self.vt = vt
        self.tolfactor_vt = tolfactor_vt
        self.verbose = verbose
        self._calc_params_init()
        for kwarg in self.default_kwargs:
            kwarg_attr = (
                kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            )
            setattr(self, kwarg, kwarg_attr)

        self.nrs = 0
        self.rs = np.array([])
        self.dr = 0
        self.zs = np.array([])
        self.dzs = np.array([])

        self.phis = np.array([])
        self.dphi = 0

        self.vs = np.array([])
        self.velchans = np.array([])

        self._fromcube_mass = 0
        self.cube = np.array([])
        self.cube_sampling = np.array([])

    def _dimension_error(self, fromcube):
        warnings.warn(
            message=f"""
The provided cube into which the model is to be build has dimensions
{np.shape(fromcube)} but the dimensions of the desired model cube is
{(self.nc, self.nys, self.nxs)}. Please, provide a cube with the right
dimensions or do not provide any cube.
""",
            category=ut.UserError,
        )

    def _required_parameter_absent_error(self):
        warnings.warn(
            message="""
Both chanwidth and vchf parameters are None. Please, provide one value of type
float to one of them.
""",
            category=ut.UserError,
        )

    def _ambiguous_input_error(self):
        warnings.warn(
            message="""
Ambiguous input. The user provided values to both chanwidth and vchf. Only one
of them should be a float, the other should be None.
 """,
            category=ut.UserError,
        )

    def _outsidegrid_warning(self):
        warnings.warn(
            message="""
Part of the model lie outside the grid of the spectral cube! The model will be
truncated or not appearing at all in your spectral cube. This is due to at
least one of three reasons:
    - The image is too small. Try to make the image larger by increasing the
      number of pixels (parameters nxs and nys), or increase the physical size
      of the image (parameter xpmax).
    - The model is far away from the image center. Try to change the reference
      pixel where the physical center (the source) is found (parameter refpix).
    - The model is outside your velocity coverage. Try to change the range of
      velocity channels of the spectral cube (parameters vch0 and vchf,
      consider negative floats if the model is blueshifted).\n
""",
            category=UserWarning,
            stacklevel=1,
        )

    def _mass_consistency_warning(self, massloss):
        warnings.warn(
            message=rf"""
The integrated mass of the cube is {massloss:.1e} % less than the input total
mass of the bowshock. This can be due to several factors:
    - Part of the model lie outside the grid of the spectral cube. If this is
      not intended, try to solve it by making the maps larger, changing the
      reference pixel to center the model in the maps, or increasing the
      velocity coverage of the spectral cube.
    - The difference between the integrated mass of the cube and the input
      total mass of the bowshock model is due to numerical errors. If you think
      that the difference is too big, you can reduce it by increasing the
      number of points of the model (inceasng nzs or/and nphis parameters).
    - The masses corresponding to a channel maps are spread along the cube in
      the velocity axis following a Gaussian distribution, being sigma equal to
      vt parameter. This distribution is truncated at vt*tolfactor_vt in order
      to make the computation substatially faster, but it can result in an
      underestimation of the integrated mass of the spectral cube. Try to make
      tolfactor_vt larger.
""",
            category=UserWarning,
            stacklevel=1,
        )

    def _sampling_xy_warning(self):
        warnings.warn(
            message="""
It is possible that the model is not well sampled in the plane of sky given the
cube dimensions and the number of model points. You can ensure a better
sampling by increasing the number of model points (nzs parameter) or decreasing
the pixel size (nxs and nys parameters).
""",
            category=UserWarning,
        )

    def _sampling_v_warning(self):
        warnings.warn(
            message="""
It is possible that the model is not well sampled in the velocity direction
given the cube dimensions and the number of model points. You can ensure a
better sampling by increasing the number of model points in the azimuthal
direction (nphis parameter) or decreasing the pixel size (nxs and nys
parameters).
""",
            category=UserWarning,
        )

    def _sampling_phi_warning(self):
        warnings.warn(
            message="""
It is possible that the model is not well sampled in the plane of sky given the
cube dimensions and the number of azimuthal points of the model (nphis). You
can ensure a better sampling by increasing the number of model points in the
azimuthal direction (nphis parameter) or decreasing the pixel size (nxs and nys
parameters).
""",
            category=UserWarning,
        )

    def _calc_params_init(self):
        if self.chanwidth is None and self.vchf is None:
            self._required_parameter_absent_error()
        elif self.chanwidth is not None and self.vchf is not None:
            self._ambiguous_input_error()
        else:
            pass

        if self.chanwidth is None:
            self.chanwidth = (self.vchf - self.vch0) / (self.nc - 1)
        elif self.vchf is None:
            self.vchf = self.chanwidth * (self.nc - 1) + self.vch0
        self.abschanwidth = np.abs(self.chanwidth)

        self.vt = (
            self.vt
            if not isinstance(self.vt, str)
            else float(self.vt.split("x")[0]) * self.chanwidth
        )
        self.arcsecpix = self.xpmax / float(self.nxs)
        if self.refpix is None:
            if self.nxs % 2 == 0:
                xref = int(self.nxs / 2)
            else:
                xref = int((self.nxs - 1) / 2)
            if self.nys % 2 == 0:
                yref = int(self.nys / 2)
            else:
                yref = int((self.nys - 1) / 2)
            self.refpix = [xref, yref]

    def _cond_populatechan(self, diffv):
        """Truncates the Gaussian distribution of the emission along the
        velocity axis"""
        if self.tolfactor_vt is not None:
            return diffv < np.abs(self.vt) * self.tolfactor_vt
        return True

    def _wvzp(self, diffv, dmass):
        """
        Weight the masses across the velocity axis using a Gaussian
        distribution
        """
        normfactor = np.abs(self.chanwidth) / (np.sqrt(np.pi) * np.abs(self.vt))
        em = dmass * np.exp(-((diffv / self.vt) ** 2)) * normfactor
        return em

    def _docic(self, chan, diffv, xpix, ypix, dxpix, dypix, dmass):
        """Cloud In Cell method"""
        em = self._wvzp(diffv, dmass)
        self.cube[chan, ypix, xpix] += em * (1 - dxpix) * (1 - dypix)
        self.cube[chan, ypix, xpix + 1] += em * dxpix * (1 - dypix)
        self.cube[chan, ypix + 1, xpix] += em * (1 - dxpix) * dypix
        self.cube[chan, ypix + 1, xpix + 1] += em * dxpix * dypix

    def _dongp(self, chan, diffv, xpix, ypix, dmass):
        """Nearest Grid Point method"""
        em = self._wvzp(diffv, dmass)
        self.cube[chan, ypix, xpix] += em

    def _sampling(self, chan, xpix, ypix):
        """Keeps track of the sampling of the spectral cube"""
        self.cube_sampling[chan, ypix, xpix] += 1

    def _particle_in_cell(self, chan, diffv, xpix, ypix, dxpix, dypix, dmass):
        """Particle in cell scheme, weather CIC or NGP"""
        if self.cic:
            self._docic(chan, diffv, xpix, ypix, dxpix, dypix, dmass)
        else:
            self._dongp(chan, diffv, xpix, ypix, dmass)

    def _populatechan(self, chan, diffv, xpix, ypix, dxpix, dypix, dmass):
        """Populates the channel"""
        if self._cond_populatechan(diffv):
            self._particle_in_cell(chan, diffv, xpix, ypix, dxpix, dypix, dmass)
            if diffv < self.abschanwidth / 2:
                self._sampling(chan, xpix, ypix)

    def _check_mass_consistency(self):
        """Checks that the input total mass of the bowshock coincides with the
        total mass of the cube"""
        print("Checking total mass consistency...")
        intmass_cube = np.sum(self.cube)
        intmass_model = self.mass + self._fromcube_mass
        mass_consistent = np.isclose(intmass_cube, intmass_model)
        massloss = (intmass_model - intmass_cube) / self.mass * 100
        if mass_consistent:
            print(
                r"""
Mass consistency test passed: The input total mass of the bowshock model
coincides with the total mass of the cube.
"""
            )
        else:
            self._mass_consistency_warning(massloss)
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
            self._sampling_xy_warning()
        if vsamp:
            self._sampling_v_warning()
        if phisamp:
            self._sampling_phi_warning()

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

        self.phis = np.linspace(0, 2 * np.pi, self.nphis + 1)[:-1]
        self.dphi = self.phis[1] - self.phis[0]

        self.vs = np.array([self.vtot(zb) for zb in self.zs])
        self.velchans = np.linspace(self.vch0, self.vchf, self.nc)
        minvelchans = np.min(self.velchans)
        maxvelchans = np.max(self.velchans)

        if fromcube is None:
            self.cube = np.zeros((self.nc, self.nys, self.nxs))
        elif (fromcube is not None) and np.shape(fromcube) == (
            (self.nc, self.nys, self.nxs)
        ):
            self.cube = np.copy(fromcube)
            self._fromcube_mass = np.sum(fromcube)
        else:
            self._dimension_error(fromcube)

        self.cube_sampling = np.zeros((self.nc, self.nys, self.nxs))

        ci = np.cos(self.i)
        si = np.sin(self.i)
        cpa = np.cos(self.pa)
        spa = np.sin(self.pa)

        outsidegrid_warning = True
        ut.progressbar_bowshock(0, self.nzs, length=50, timelapsed=0, intervaltime=0)
        for iz, z in enumerate(self.zs):
            if self.verbose:
                t0 = datetime.now()

            if iz == 0:
                # Treat boundary of the outer parts of the bowshock wings
                intmass = self.intmass_analytical(self.rbf)
                intmass_halfdr = self.intmass_analytical(self.rbf - self.dr / 2)
                dmass = (intmass - intmass_halfdr) / self.nphis
            elif iz == len(self.zs) - 1:
                # Treat head boundary
                dmass = self.intmass_analytical(self.dr / 2) / self.nphis
            else:
                # Treat the rest of the bowshock
                dmass = self.dmass_func(z, self.dzs[iz], self.dphi)

            for phi in self.phis:
                # for phi in self.phis+self.dphi*np.random.rand():
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
                condition_inside_map = (
                    (xpix + 1 < self.nxs)
                    and (ypix + 1 < self.nys)
                    and (xpix > 0)
                    and (ypix > 0)
                )
                condition_inside_velcoverage = minvelchans <= vlsr <= maxvelchans
                if condition_inside_map and condition_inside_velcoverage:
                    dxpix = xpixcoord - xpix
                    dypix = ypixcoord - ypix
                    for chan, vchan in enumerate(self.velchans):
                        diffv = np.abs(vlsr - vchan)
                        self._populatechan(chan, diffv, xpix, ypix, dxpix, dypix, dmass)
                else:
                    if outsidegrid_warning:
                        self._outsidegrid_warning()
                        outsidegrid_warning = False
            if self.verbose:
                tf = datetime.now()
                intervaltime = (tf - t0).total_seconds()
                ts.append(intervaltime)
                ut.progressbar_bowshock(
                    iz + 1, self.nzs, np.sum(ts), intervaltime, length=50
                )
        _ = self._check_mass_consistency()
        self._check_sampling()

    def plot_channel(
        self,
        chan,
        vmax=None,
        vmin=None,
        cmap="inferno",
        savefig=None,
        return_fig_axs=False,
    ):
        """
        Plots a channel map of a cube

        Parameters
        ----------
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
        return_fig_axs : bool, optional
            If True, returns the figure, axes of the channel map, and the axes
            the colorbar.  If False, does not return anything.

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
        fig, axs, cbax = pl.plot_channel(
            cube=self.cube,
            chan=chan,
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            vmax=vmax,
            vmin=vmin,
            cmap=cmap,
            units=r"Mass / (pixel $\times$ channel) [Msun]",
            refpix=self.refpix,
            return_fig_axs=True,
        )
        if return_fig_axs:
            return fig, axs, cbax
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")

    def plot_channels(self, savefig=None, return_fig_axs=False, **kwargs):
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
        """
        fig, axs, cbax = pl.plot_channels(
            cube=self.cube,
            arcsecpix=self.arcsecpix,
            velchans=self.velchans,
            units=r"Mass / (pixel $\times$ channel) [Msun]",
            refpix=self.refpix,
            return_fig_axs=True,
            **kwargs,
        )
        if return_fig_axs:
            return fig, axs, cbax
        if savefig is not None:
            fig.savefig(savefig, bbox_inches="tight")
