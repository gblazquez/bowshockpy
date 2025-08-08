import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib import cm
from matplotlib import colors
from matplotlib import colormaps

from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture, RectangularAperture

from itertools import product
from matplotlib import ticker

import bowshockpy.utils as ut

class BowshockModelPlot():
    """
    Figure including the main parameters of the bowshock model, its morphology
    and kinematics, and the distribution of the surface density
    
    Parameters:
    -----------
    bsm : `~bowshockpy.models.BowshockModel` class instance
        Instance of the model to plot
    modelname : str, optional 
        Name of the model to include in the plot
    nzs : int, optional 
        Number of points used to compute the model solutions
    figsize: tuple, optional
        Tuple passed to `matplotib.pyplot.figure` to define the dimensions of
        the figure
    narrows: int, optional
        Number of arrows to show in order to indicate the velocity at each
        symmetrical half of the model.        
    v_arrow_ref: float, optional
        Velocity in km/s to use as reference in the reference arrow
    linespacing : float, optional
        Spacing between the text lines
    textbox_widthratio : float, optional
        Width ratio of the text ax to pass to GridSpec

    Attributes:
    -----------
    nzs : int
        Number of points used to compute the model solutions
    zs : numpy.ndarray
        Array of the z-coordinates of the model.
    dzs : numpy.ndarray
        Increment of z-coordinates between the points.
    Rs : numpy.ndarray
        Array with the radii of the model at each z-coordinate [km].
    thetas : numpy.ndarray
        Array of the polar angle of the position vector at each point of the
        model [radians].
    vs : numpy.ndarray
        Array of the total velocity for each point of the model [km/s].
    vrs : numpy.ndarray
        Array of the radial component of the velocity at each point of the model
        [km/s].
    vzs : numpy.ndarray
        Array of the z-coordinate component of the velocity at each point of the
        model [km/s].
    surfdenss : numpy.ndarray 
        Array of the surfance density of the shell at each z-coordinate [Msun
        km-2]
    surfdenss_gcm2 : numpy.ndarray
        Array of the surfance density of the shell at each z-coordinate [g cm-2]
    axs : dict
        Dictionary of matplotlib.axes.Axes in the figure
    cbaxs : dict
        Dictionary of matplotlib.axes.Axes of the colorbars in the figure
    """

    def __init__(
            self, bsm, modelname="none", nzs=200,
            figsize=(16,3), narrows=10, v_arrow_ref=100,
            linespacing=0.08, textbox_widthratio=0.7,
            ):
        self.mo = bsm
        self.modelname = modelname
        self.nzs = nzs
        self.nrs = nzs
        self.narrows = narrows 
        self.figsize = figsize
        self.v_arrow_ref = v_arrow_ref
        self.linespacing = linespacing
        self.textbox_widthratio = textbox_widthratio

        self.zs = np.array([])
        self.dzs = np.array([])
        self.Rs = np.array([])

        self.vrs = np.array([])
        self.vzs = np.array([])
        self.vs = np.array([])
        self.thetas = np.array([])
        self.Rs_arcsec = np.array([])
        self.zs_arcsec = np.array([])

        self.surfdenss = np.array([])
        self.surfdenss_gcm2 = np.array([])

        self.zs_arrows = np.array([])
        self.Rs_arrows = np.array([])
        self.zs_arrows_tip = np.array([])
        self.Rs_arrows_tip = np.array([])
        self.z_arrow_ref_tip = None
        self.R_arrow_ref_tip = None

        self.axs = {}
        self.cbaxs = {}
        self.maxsurfdens_plot = None
        self.minsurfdens_plot = None

        self._calc_solutions()
        self._calc_arrows()
        # self._create_axes()
        # self.plot()

    def _calc_solutions(self):
        # self.zsextended = self.zb_r(
        #     np.linspace(self.rbf, self.rj, self.nzs)
        # )
        # self.nrs = self.nzs
        self.rs = np.linspace(self.mo.rbf, 0, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.mo.zb_r(self.rs)
        self.dzs = self.mo.dz_func(self.mo.zb_r(self.rs), self.dr)

        self.vs = np.array([self.mo.vtot(zb) for zb in self.zs])
        self.Rs = np.array([self.mo.rb(zb) for zb in self.zs])
        self.vrs = np.array([self.mo.vr(zb) for zb in self.zs])
        self.vzs = np.array([self.mo.vz(zb) for zb in self.zs])
        self.vs = np.array([np.sqrt(vrr**2+vzz**2) for vrr, vzz in zip(self.vrs, self.vzs)])
        self.maxvs = np.max(self.vs)
        self.minvs = np.min(self.vs)
        self.thetas = np.array([np.arctan(self.Rs[i] / z)
                       for i, z in enumerate(self.zs)])

        self.Rs_arcsec = self.mo.km2arcsec(self.Rs)
        self.zs_arcsec = self.mo.km2arcsec(self.zs)
 
        self.surfdenss = np.array([self.mo.surfdens(zb) for zb in self.zs])
        self.surfdenss_gcm2 = self.mo.solMasskm2togcm2(self.surfdenss)

    def _calc_arrows(self):
        idx_arr = int(len(self.zs_arcsec)/self.narrows)
        self.larrow = 1 / np.max([np.max(self.vrs), np.max(self.vzs)])

        self.zs_arrows = self.zs_arcsec[::idx_arr]
        self.Rs_arrows = self.Rs_arcsec[::idx_arr]

        self.zs_arrows_tip = self.zs_arrows + self.larrow * self.vzs[::idx_arr]
        self.Rs_arrows_tip = self.Rs_arrows + self.larrow * self.vrs[::idx_arr]

    def _create_axes(self):
        nrow = 1
        ncol = 3
        wspace = 0.2
        hspace = 0.4
        width_ratios = [self.textbox_widthratio, 1, 1]
        height_ratios = [1] * nrow

        self.fig_model = plt.figure(figsize=self.figsize)
        gs = GridSpec(
            nrow, ncol,
            height_ratios=height_ratios,
            width_ratios=width_ratios,
            hspace=hspace, wspace=wspace
        )
        gss = {}
        gss[0] = gs[0, 1].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gss[1] = gs[0, 2].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        # gss[2] = gs[1, 1].subgridspec(
        #     2, 1,
        #     height_ratios=[0.05, 1],
        #     width_ratios=[1],
        #     hspace=0.05,
        # )

        self.axs["text"] = plt.subplot(gs[:, 0])
        self.axs[0] = plt.subplot(gss[0][1, 0])
        self.cbaxs[0] = plt.subplot(gss[0][0, 0])
        self.axs[1] = plt.subplot(gss[1][1, 0])
        self.cbaxs[1] = plt.subplot(gss[1][0, 0])
        self.axs["text"].set_axis_off()

    def plot(self):
        """
        Plots the 2D bowshock model
        """
        self._create_axes()
        showtext = \
            fr"""
            {self.modelname}

            $v_\mathrm{{j}} = {{{self.mo.vj:.2f}}}$ km/s
            $v_0 = {{{self.mo.v0:.2f}}}$ km/s
            $v_a = {{{self.mo.va:.2f}}}$ km/s
            $L_0 = {{{self.mo.L0_arcsec:.2f}}}$ arcsec
            $z_\mathrm{{j}} = {{{self.mo.zj_arcsec:.2f}}}$ arcsec
            $r_\mathrm{{b,f}} = {{{self.mo.rbf_arcsec:.2f}}}$ arcsec
            $m$ = ${{{self.mo.mass*10**4:.2f}}}\times10^{{-4}}$ M$_\odot$
            $t_\mathrm{{j}} = {{{self.mo.tj_yr:.2f}}}$ yr
            $\rho_a = {{{self.mo.rhoa_gcm3*10**20:.2f}}}\times 10^{{-20}}$ g cm$^{{-3}}$
            $\dot{{m}}_0 = {{{self.mo.mp0_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            $\dot{{m}}_{{a,f}} = {{{self.mo.mpamb_f_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 1.05-self.linespacing*n, line, fontsize=10,
                              transform=self.axs["text"].transAxes)

        """
        Deprojected shell Morph. and Kin., color velocity
        """
        cmap = "turbo_r"
        for i, zarcsec in enumerate(self.zs_arcsec):
            c = ut.get_color(
                [self.minvs, self.maxvs],
                self.vs[i],
                cmap,
            )
            self.axs[0].plot(
                zarcsec,
                self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
            self.axs[0].plot(
                zarcsec,
                -self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.Normalize(
                       vmax=self.maxvs,
                       vmin=self.minvs),
                   cmap=cmap,
               ),
               cax=self.cbaxs[0],
               orientation="horizontal",
        )

        for i in range(len(self.zs_arrows)):
            self.axs[0].annotate(
                "",
                xy=(self.zs_arrows_tip[i], self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))
            self.axs[0].annotate(
                "",
                xy=(self.zs_arrows_tip[i], -self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], -self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))

        xlims = [np.min(self.zs_arcsec),
                np.max(np.max(self.zs_arrows_tip))*1.1]
        ylims = [-np.max(self.Rs_arcsec)*1.3, np.max(self.Rs_arcsec)*1.3]
        self.axs[0].set_xlim(xlims)
        self.axs[0].set_ylim(ylims)
        larrow_scaled = self.larrow * self.v_arrow_ref
        self.z_arrow_ref = xlims[1]*0.97 - larrow_scaled
        self.R_arrow_ref = ylims[0] + np.diff(ylims)*0.05
        self.z_arrow_ref_tip = self.z_arrow_ref + larrow_scaled 
        self.R_arrow_ref_tip = self.R_arrow_ref + self.larrow * 0
        self.axs[0].annotate(
            "",
            xy=(self.z_arrow_ref_tip, self.R_arrow_ref_tip),
            xytext=(self.z_arrow_ref, self.R_arrow_ref),
            arrowprops=dict(arrowstyle="->"))

        self.axs[0].text(
            self.z_arrow_ref+0.0,
            self.R_arrow_ref+0.05,
            f"{self.v_arrow_ref:d} km/s"
        )

        self.axs[0].set_aspect("equal")
        self.axs[0].set_xlabel("Distance [arcsec]")
        self.axs[0].set_ylabel("Radius [arcsec]")

        self.cbaxs[0].tick_params(
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[0].set_xlabel(r"Speed [km/s]", )
        self.cbaxs[0].xaxis.set_label_position('top')


        """
        Deprojected shell Morph. and Kin., color density
        """
        self.minsurfdens_plot = np.percentile(self.surfdenss_gcm2[:-1], 0) 
        self.maxsurfdens_plot = np.percentile(self.surfdenss_gcm2[:-1], 70)
        norm = colors.LogNorm(
                    vmax=self.maxsurfdens_plot,
                    vmin=self.minsurfdens_plot,
        #           linthresh=self.maxsurfdens*0.99,
                    )
        cmap = "viridis"
        # we skip the point at the tip, there is a discontinuity and the surface
        # density is 0
        for i, zarcsec in enumerate(self.zs_arcsec[:-1]):
            c = ut.get_color(
                [self.minsurfdens_plot, self.maxsurfdens_plot],
                self.surfdenss_gcm2[i],
                cmap,
                customnorm=norm, 
            )
            self.axs[1].plot(
                zarcsec,
                self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
            self.axs[1].plot(
                zarcsec,
                -self.Rs_arcsec[i],
                color=c,
                marker="o",
            )
        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=norm,
                   cmap=cmap,
               ),
               cax=self.cbaxs[1],
               orientation="horizontal",
        )
        self.cbaxs[1].tick_params(
            axis="x", which="both", top=True, bottom=False
        )
#        self.cbaxs[1].set_xscale("log")

        for i in range(len(self.zs_arrows)):
            self.axs[1].annotate(
                "",
                xy=(self.zs_arrows_tip[i], self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))
            self.axs[1].annotate(
                "",
                xy=(self.zs_arrows_tip[i], -self.Rs_arrows_tip[i]),
                xytext=(self.zs_arrows[i], -self.Rs_arrows[i]),
                arrowprops=dict(arrowstyle="->"))

        xlims = [np.min(self.zs_arcsec),
                np.max(np.max(self.zs_arrows_tip))*1.1]
        ylims = [-np.max(self.Rs_arcsec)*1.3, np.max(self.Rs_arcsec)*1.3]
        self.axs[1].set_xlim(xlims)
        self.axs[1].set_ylim(ylims)
        larrow_scaled = self.larrow * self.v_arrow_ref
        self.z_arrow_ref = xlims[1]*0.97 - larrow_scaled
        self.R_arrow_ref = ylims[0] + np.diff(ylims)*0.05
        self.z_arrow_ref_tip = self.z_arrow_ref + larrow_scaled 
        self.R_arrow_ref_tip = self.R_arrow_ref + self.larrow * 0
        self.axs[1].annotate(
            "",
            xy=(self.z_arrow_ref_tip, self.R_arrow_ref_tip),
            xytext=(self.z_arrow_ref, self.R_arrow_ref),
            arrowprops=dict(arrowstyle="->"))

        self.axs[1].text(
            self.z_arrow_ref+0.0,
            self.R_arrow_ref+0.05,
            f"{self.v_arrow_ref:d} km/s"
        )

        self.axs[1].set_aspect("equal")
        self.axs[1].set_xlabel("Distance [arcsec]")
        #self.axs[1].set_ylabel("Radius [arcsec]")

        self.cbaxs[1].set_xlabel(r"Surface density [g cm$^{-2}$]", )
        self.cbaxs[1].xaxis.set_label_position('top')
        self.cbaxs[1].tick_params(
            axis="x",
             bottom=False, labelbottom=False,
             top=True, labeltop=True
         )
        # For some reason tick_params is not able to plot the ticks above if the
        # tick lables are in scientific notation. I have to do:
        self.cbaxs[1].xaxis.set_ticks_position('top')
        self.cbaxs[1].xaxis.set_label_position('top') 

    def savefig(self, figname=None, **kwargs):
        """
        Saves the plot of the bowhsock model.       

        Parameters
        ----------
        figname : str, optional
            Full path name of the figure. If None, the the full path name will
            be models/{self.modelname}/bowshock_model.pdf. If the folder tree does not exist, it will be created. 
        """
        if figname is None:
            ut.make_folder(f"models/{self.modelname}")
            figname = f"models/{self.modelname}/bowshock_model.pdf"
        self.fig_model.savefig(f"{figname}", bbox_inches="tight", **kwargs)



class BowshockObsModelPlot():
    """
    Figure including the main parameters of the bowshock model, its projected
    morphology, kinematics, and a PV diagram along the symmetry axis with the
    distribution of the surface density in color code.
    
    Parameters:
    -----------
    bsm : `~bowshockpy.models.BowshockModel` class instance
        Instance of the model to plot
    modelname : str, optional
        Name of the model to include in the plot
    nzs : int, optional
        Number of z coordinates used to compute the model solutions
    nphis : int, optional
        Number of phi coordinates used to compute the model solutions
    figsize: tuple, optional
        Tuple passed to `matplotib.pyplot.figure` to define the dimensions of
        the figure
    linespacing : float, optional
        Spacing between the text lines
    textbox_widthratio : float, optional
        Width ratio of the text ax to pass to GridSpec
    cmap : str, optional
        Colormap label
    minpointsize : float, optional
        Minsize of the points to plot
    maxpointsize : float, optional
        Minsize of the points to plot

    Attributes:
    -----------
    nrs : int
        Number of r coordinates used to compute the model solutions
    zs : numpy.ndarray
        Array of the z-coordinates of the model.
    dzs : numpy.ndarray
        Increment of z-coordinates between the points.
    Rs : numpy.ndarray
        Array with the radii of the model at each z-coordinate [km].
    thetas : numpy.ndarray
        Array of the polar angle of the position vector at each point of the
        model [radians].
    vs : numpy.ndarray
        Array of the total velocity for each point of the model [km/s].
    vrs : numpy.ndarray
        Array of the radial component of the velocity at each point of the model
        [km/s].
    vzs : numpy.ndarray
        Array of the z-coordinate component of the velocity at each point of the
        model [km/s].
    surfdenss : numpy.ndarray 
        Array of the surfance density of the shell at each z-coordinate [Msun
        km-2]
    surfdenss_gcm2 : numpy.ndarray
        Array of the surfance density of the shell at each z-coordinate [g cm-2]
    axs : dict
        Dictionary of matplotlib.axes.Axes in the figure
    cbaxs : dict
        Dictionary of matplotlib.axes.Axes of the colorbars in the figure
    """

    def __init__(
            self, bsmobs,
            modelname="none",
            nzs=150,
            nphis=150,
            figsize=(12,6), 
            linespacing=0.09,
            textbox_widthratio=0.8,
            cmap="turbo",
            minpointsize=0.1,
            maxpointsize=10,
            ):
        self.mo = bsmobs
        self.modelname = modelname
        self.nzs = nzs
        self.nphis = nphis
        self.nrs = nzs
        self.nphis = nphis
        self.figsize = figsize
        self.linespacing = linespacing
        self.textbox_widthratio = textbox_widthratio

        self.zs = np.array([])
        self.dzs = np.array([])
        self.Rs = np.array([])

        self.vrs = np.array([])
        self.vzs = np.array([])
        self.vs = np.array([])
        self.thetas = np.array([])
        self.Rs_arcsec = np.array([])
        self.zs_arcsec = np.array([])

        self.surfdenss = np.array([])
        self.surfdenss_gcm2 = np.array([])

        self.xps_phi90 = np.array([])
        self.xps_phi0 = np.array([])
        self.xps_phi180 = np.array([])
        self.vloss_phi0 = np.array([])
        self.vloss_phi90 = np.array([])
        self.vloss_phi180 = np.array([])
        self.xps_phi0_arcsec = np.array([])
        self.xps_phi90_arcsec = np.array([])
        self.xps_phi180_arcsec = np.array([])
        self.maxvlos = None
        self.minvlos = None

        self.axs = {}
        self.cbaxs = {}
        self.cmap = cmap
        self.minpointsize = minpointsize
        self.maxpointsize = maxpointsize
        self.minsurfdens_plot = None 
        self.maxsurfdens_plot = None 

        self._calc_solutions()
        # self._create_axes()
        # self.plot()

    def _calc_solutions(self):
        self.rs = np.linspace(self.mo.rbf, 0, self.nrs)
        self.dr = self.rs[0] - self.rs[1]
        self.zs = self.mo.zb_r(self.rs)
        self.dzs = self.mo.dz_func(self.mo.zb_r(self.rs), self.dr)

        self.vs = np.array([self.mo.vtot(zb) for zb in self.zs])
        self.Rs = np.array([self.mo.rb(zb) for zb in self.zs])
        self.vrs = np.array([self.mo.vr(zb) for zb in self.zs])
        self.vzs = np.array([self.mo.vz(zb) for zb in self.zs])
        self.vs = np.array([np.sqrt(vrr**2+vzz**2) for vrr, vzz in zip(self.vrs, self.vzs)])
        self.maxvs = np.max(self.vs)
        self.minvs = np.min(self.vs)
        self.thetas = np.array([np.arctan(self.Rs[i] / z)
                       for i, z in enumerate(self.zs)])

        self.Rs_arcsec = self.mo.km2arcsec(self.Rs)
        self.zs_arcsec = self.mo.km2arcsec(self.zs)
 
        self.surfdenss = np.array([self.mo.surfdens(zb) for zb in self.zs])
        self.surfdenss_gcm2 = self.mo.solMasskm2togcm2(self.surfdenss)

        self.xps_phi90 = np.array([self.mo.xp(zb,phi=np.pi/2) for zb in self.zs])
        self.xps_phi0 = np.array([self.mo.xp(zb,phi=0) for zb in self.zs])
        self.xps_phi180 = np.array([self.mo.xp(zb,phi=np.pi) for zb in self.zs])
        self.vloss_phi0 = -np.array([self.mo.vzp(zb,phi=0) for zb in self.zs])
        self.vloss_phi90 = -np.array([self.mo.vzp(zb,phi=np.pi/2) for zb in self.zs])
        self.vloss_phi180 = -np.array([self.mo.vzp(zb,phi=np.pi) for zb in self.zs])
        self.maxvlos = np.max([self.vloss_phi0, self.vloss_phi180])
        self.minvlos = np.min([self.vloss_phi0, self.vloss_phi180])

        self.xps_phi0_arcsec = self.mo.km2arcsec(self.xps_phi0)
        self.xps_phi90_arcsec = self.mo.km2arcsec(self.xps_phi90)
        self.xps_phi180_arcsec = self.mo.km2arcsec(self.xps_phi180)

        phi_0_1 = -np.pi/2 
        phi_f_1 = +np.pi/2
 
        nphis_half = int(self.nphis/2)
        self.phis_1 = np.linspace(phi_0_1, phi_f_1, nphis_half)[:-1]

        self.xp_zs_phis_1 = np.array([
            [self.mo.xp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_1]
        )
        self.yp_zs_phis_1 = np.array([
            [self.mo.yp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_1]
        )
        self.vlos_zs_phis_1 = np.array([
            [-self.mo.vzp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_1]
        )

        phi_0_2 = +np.pi/2 
        phi_f_2 = +np.pi*3/2
 
        self.phis_2 = np.linspace(phi_0_2, phi_f_2, nphis_half)[:-1]

        self.xp_zs_phis_2 = np.array([
            [self.mo.xp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_2]
        )
        self.yp_zs_phis_2 = np.array([
            [self.mo.yp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_2]
        )
        self.vlos_zs_phis_2 = np.array([
            [-self.mo.vzp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_2]
        )

        phi_0_3 = 0
        phi_f_3 = np.pi
        
        nphis_half = int(self.nphis/2)
        self.phis_3 = np.linspace(phi_0_3, phi_f_3, nphis_half)[:-1]
        
        self.xp_zs_phis_3 = np.array([
            [self.mo.xp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_3]
        )
        self.zp_zs_phis_3 = np.array([
            [self.mo.zp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_3]
        )
        self.vlos_zs_phis_3 = np.array([
            [-self.mo.vzp(zb, phi)
             for zb in self.zs]
            for phi in self.phis_3]
        )

    def _create_axes(self):

        nrow = 2
        ncol = 1
        wspace = 0.15
        hspace = 0.4
        width_ratios = [1,]
        height_ratios = [1, 1] 
        
        self.fig_model = plt.figure(figsize=self.figsize)
        gs = GridSpec(
            nrow, ncol,
            height_ratios=height_ratios,
            width_ratios=width_ratios,
            hspace=hspace, wspace=wspace
        )
        
        gss = {}
        gss[0] = gs[0, 0].subgridspec(
            1, 3,
            height_ratios=[1],
            width_ratios=[0.75, 1, 0.5],
            hspace=0.05,
        )
        gss[1] = gs[1, 0].subgridspec(
            1, 3,
            height_ratios=[1],
            width_ratios=[1, 1, 1],
            hspace=0.05,
        )
        
        gsss = {}
        gsss[0] = gss[0][0, 0].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gsss[1] = gss[0][0, 1:2].subgridspec(
            1, 2,
            height_ratios=[1],
            width_ratios=[1, 0.05],
            wspace=0.05,
        )
        gsss[2] = gss[1][0, 0].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gsss[3] = gss[1][0, 1].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        gsss[4] = gss[1][0, 2].subgridspec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.05,
        )
        
        self.axs["text"] = plt.subplot(gsss[0][:, 0])
        self.axs[0] = plt.subplot(gsss[1][0, 0])
        self.cbaxs[0] = plt.subplot(gsss[1][0, 1])
        self.axs[1] = plt.subplot(gsss[2][:, 0])
        self.axs[2] = plt.subplot(gsss[3][:, 0])
        self.axs[3] = plt.subplot(gsss[4][:, 0])
        self.axs["text"].set_axis_off()

    def plot(self):
        """
        Plots the 2D bowshock model
        """
        self._create_axes()
        showtext = \
            fr"""
            {self.modelname}
            $i = {{{self.mo.i*180/np.pi:.2f}}}^\circ$
            $v_\mathrm{{vsys}} = {{{self.mo.vsys:.2f}}}$ km/s
            $v_\mathrm{{j}} = {{{self.mo.vj:.2f}}}$ km/s
            $v_0 = {{{self.mo.v0:.2f}}}$ km/s
            $v_a = {{{self.mo.va:.2f}}}$ km/s
            $L_0 = {{{self.mo.L0_arcsec:.2f}}}$ arcsec
            $z_\mathrm{{j}} = {{{self.mo.zj_arcsec:.2f}}}$ arcsec
            $r_\mathrm{{b,f}} = {{{self.mo.rbf_arcsec:.2f}}}$ arcsec
            $m$ = ${{{self.mo.mass*10**4:.2f}}}\times10^{{-4}}$ M$_\odot$
            $t_\mathrm{{j}} = {{{self.mo.tj_yr:.2f}}}$ yr
            $\rho_a = {{{self.mo.rhoa_gcm3*10**20:.2f}}}\times 10^{{-20}}$ g cm$^{{-3}}$
            $\dot{{m}}_0 = {{{self.mo.mp0_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            $\dot{{m}}_{{a,f}} = {{{self.mo.mpamb_f_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 1.05-self.linespacing*n, line, fontsize=10,
                              transform=self.axs["text"].transAxes)

        
        """
        Projected shell Morph. and Kin.
        """
        
        # op controls the plotting order of the points 
        op = 1 if self.mo.i <= np.pi/2 else -1
        norm = colors.Normalize(
                    vmax=self.maxvlos + self.mo.vsys,
                    vmin=self.minvlos + self.mo.vsys,
                )
 
        range_point_sizes = np.linspace(
            self.maxpointsize,
            self.minpointsize,
            len(self.zs)
            )[::op]
        point_sizes = [[i]*len(self.vlos_zs_phis_1) for i in range_point_sizes]
        xps_arcsec_1 = self.mo.km2arcsec(self.xp_zs_phis_1.T[::op])
        yps_arcsec_1 = self.mo.km2arcsec(self.yp_zs_phis_1.T[::op])
        vlos_1 = self.vlos_zs_phis_1.T[::op] + self.mo.vsys
        self.axs[1].scatter(
            xps_arcsec_1,
            yps_arcsec_1,
            c=vlos_1,
            cmap=self.cmap,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            s=point_sizes,
        )
        xps_arcsec_2 = self.mo.km2arcsec(self.xp_zs_phis_2.T[::op])
        yps_arcsec_2 = self.mo.km2arcsec(self.yp_zs_phis_2.T[::op])
        vlos_2 = self.vlos_zs_phis_2.T[::op] + self.mo.vsys
        self.axs[1].scatter(
            xps_arcsec_2,
            yps_arcsec_2,
            c=vlos_2,
            cmap=self.cmap,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            s=point_sizes,
        )
        self.axs[1].set_aspect("equal")
        oxlim = self.axs[1].get_xlim()
        oylim = self.axs[1].get_ylim()
        self.axs[1].set_xlim([oxlim[0], oxlim[1]+0.6])
        self.axs[1].set_ylim([oylim[0], oylim[1]+0.3])
        text_obj = self.axs[1].text(
            0.64, 0.75, "Observer's\n view", transform=self.axs[1].transAxes)

        # text_obj = self.axs[1].text(
        #     0.75, 0.83, "Observer\n view", transform=self.axs[1].transAxes)
        # display_text_f = text_obj.get_window_extent().get_points()[-1]
        # xtext_f, ytext_f = \
        #     self.axs[1].transData.inverted().transform(display_text_f)
        # oxlim = self.axs[1].get_xlim()
        # self.axs[1].set_xlim([oxlim[0], np.max([xtext_f+1, oxlim[1]])])

        self.axs[1].set_xlabel("Projected length [arcsec]")
        self.axs[1].set_ylabel("Proj. radius [arcsec]")

        range_point_sizes = np.linspace(
            self.maxpointsize,
            self.minpointsize,
            len(self.zs)
            )[::-op]
        point_sizes = [[i]*len(self.vlos_zs_phis_1) for i in range_point_sizes]
        xps_arcsec_1 = self.mo.km2arcsec(self.xp_zs_phis_2.T[::-op])
        yps_arcsec_1 = self.mo.km2arcsec(self.yp_zs_phis_2.T[::-op])
        vlos_1 = self.vlos_zs_phis_2.T[::-op] + self.mo.vsys
        self.axs[2].scatter(
            xps_arcsec_1,
            yps_arcsec_1,
            c=vlos_1,
            cmap=self.cmap,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            s=point_sizes,
        )
        xps_arcsec_2 = self.mo.km2arcsec(self.xp_zs_phis_1.T[::-op])
        yps_arcsec_2 = self.mo.km2arcsec(self.yp_zs_phis_1.T[::-op])
        vlos_2 = self.vlos_zs_phis_1.T[::-op] + self.mo.vsys
        self.axs[2].scatter(
            xps_arcsec_2,
            yps_arcsec_2,
            c=vlos_2,
            cmap=self.cmap,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            s=point_sizes,
        )
        self.axs[2].set_aspect("equal")
        oxlim = self.axs[2].get_xlim()
        oylim = self.axs[2].get_ylim()
        self.axs[2].set_xlim([oxlim[0], oxlim[1]+0.6])
        self.axs[2].set_ylim([oylim[0], oylim[1]+0.3])
        text_obj = self.axs[2].text(
            0.8, 0.85, "Back", transform=self.axs[2].transAxes)
        # display_text_f = text_obj.get_window_extent().get_points()[-1]
        # xtext_f, ytext_f = \
        #     self.axs[2].transData.inverted().transform(display_text_f)
        # oxlim = self.axs[2].get_xlim()
        # self.axs[2].set_xlim([oxlim[0], np.max([xtext_f+1, oxlim[1]])])

        self.axs[2].set_xlabel("Projected length [arcsec]")

        range_point_sizes = np.linspace(
            self.maxpointsize,
            self.minpointsize,
            len(self.zs)
            )[::1]
        point_sizes = [[i]*len(self.vlos_zs_phis_1) for i in range_point_sizes]
        
        xps = self.xp_zs_phis_3.T[::1]
        zps = self.zp_zs_phis_3.T[::1]
        
        xps_rot = xps*np.sin(self.mo.i) + zps * np.cos(self.mo.i)
        zps_rot = -(xps*np.cos(self.mo.i) - zps * np.sin(self.mo.i))
        xps_rot_arcsec = self.mo.km2arcsec(xps_rot)
        zps_rot_arcsec = self.mo.km2arcsec(zps_rot)
        
        self.axs[0].scatter(
            xps_rot_arcsec,
            zps_rot_arcsec,
            c=self.vlos_zs_phis_3.T[::1] + self.mo.vsys,
            cmap=self.cmap,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            s=point_sizes,
        )
        
        self.axs[0].plot(0, 0, "*k")
        
        larrow = 0.5
        
        # Necessary to transform from display to physical coordinates
        self.axs[0].set_xlim(self.axs[0].get_xlim())
        self.axs[0].set_ylim(self.axs[0].get_ylim())
        display_coords_arrow = self.axs[0].transAxes.transform((0.98, 0.95))
        zarrow, Rarrow =  self.axs[0].transData.inverted().transform(display_coords_arrow)
        display_coords_text = self.axs[0].transAxes.transform((0.85, 0.83))
        xtext, ytext =  self.axs[0].transData.inverted().transform(display_coords_text)

        zarrow_tip = larrow*np.cos(self.mo.i) + zarrow
        Rarrow_tip = larrow*np.sin(self.mo.i) + Rarrow

        self.axs[0].annotate(
            "",
            xy=(zarrow_tip, Rarrow_tip,),
            xytext=(zarrow, Rarrow,),
            arrowprops=dict(arrowstyle="->")
            )

        text_obj = self.axs[0].text(xtext, ytext, "To observer")
        display_text_f = text_obj.get_window_extent().get_points()[-1]
        xtext_f, ytext_f = self.axs[0].transData.inverted().transform(display_text_f)

        oxlim = self.axs[0].get_xlim()
        oylim = self.axs[0].get_ylim()
        
        x_newlim = np.max([zarrow, zarrow_tip, oxlim[1], xtext_f])
        y_newlim = np.max([Rarrow, Rarrow_tip, oylim[1]])
        
        self.axs[0].set_xlim([oxlim[0], x_newlim*1.1])
        self.axs[0].set_ylim([oylim[0], y_newlim*1.1])
        
        self.axs[0].set_aspect("equal")
        self.axs[0].set_xlabel(r"Length [arcsec]")
        self.axs[0].set_ylabel(r"Radius [arcsec]")

        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=norm,
                   cmap=self.cmap,
               ),
               cax=self.cbaxs[0],
               orientation="vertical",
        )
        
        self.cbaxs[0].tick_params(
            which="both",
            left=False, labelleft=False,
            right=True, labelright=True
        )
        self.cbaxs[0].set_ylabel(r"Line-of-sight velocity [km/s]", )
        self.cbaxs[0].yaxis.set_label_position('right')



        """
        PV diagram: projected velocity
        """
        self.axs[3].scatter(
            self.xps_phi180_arcsec,
            self.vloss_phi180 + self.mo.vsys,
            marker="o",
            c=self.vloss_phi180 + self.mo.vsys,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            cmap=self.cmap,
        )
        self.axs[3].scatter(
            self.xps_phi0_arcsec,
            self.vloss_phi0 + self.mo.vsys,
            marker="o",
            c=self.vloss_phi0 + self.mo.vsys,
            vmax=self.maxvlos + self.mo.vsys,
            vmin=self.minvlos + self.mo.vsys,
            cmap=self.cmap,
            label="PV diagram\n along axis"
        )
        allvelsarray = np.array([
            self.vloss_phi0[:-1] + self.mo.vsys,
            self.vloss_phi180[:-1] + self.mo.vsys]).ravel()
        argmaxvelpv = np.argmax(np.abs(allvelsarray))
        if allvelsarray[argmaxvelpv]<0:
            self.axs[3].invert_yaxis()
        else:
            pass

        self.axs[3].set_box_aspect(1) 
        self.axs[3].set_xlabel("Projected length [arcsec]")
        self.axs[3].set_ylabel("Line-of-sight velocity [km/s]")
        # self.axs[3].text(0.1, 0.8, "PV diagram\n along axis",
        #                 transform=self.axs[3].transAxes)

        self.axs[3].legend(frameon=False, markerscale=0)
         
    def savefig(self, figname=None, **kwargs):
        """
        Saves the plot of the bowhsock model.       

        Parameters
        ----------
        figname : str, optional
            Full path name of the figure. If None, the the full path name will
            be models/{self.modelname}/bowshock_projected.pdf. If the folder tree does not exist, it will be created. 
        """
        if figname is None:
            ut.make_folder(f"models/{self.modelname}")
            figname = f"models/{self.modelname}/bowshock_projected.jpg"
        self.fig_model.savefig(f"{figname}", bbox_inches="tight", **kwargs)

def add_beam_to_ax(
        ax, bmin, bmaj, pabeam, extent, 
        add_beam_rectangle=True, beam_factor_rectangle=1.5,
        beam_color_rectangle="k", beam_edgecolor_rectangle="w",
        beam_fill_rectangle=False, beam_color="w", beam_fill=True
        ):
    # in radians
    pa = pabeam * np.pi/180. + np.pi/2
    # semi-major axis in pixels
    a = bmaj  
    # semi-minor axis in pixels
    b = bmin 
    xbpos = extent[0] - np.abs(extent[0]-extent[1])/20 - np.max([a,b])
    ybpos = extent[2] + np.abs(extent[2]-extent[3])/20 + np.max([a,b])
    geometry = EllipseGeometry(
        x0=xbpos, y0=ybpos, sma=a*0.5, eps=(1-b/a), pa=np.pi-pa)
    if add_beam_rectangle:
        aper_rectangle = RectangularAperture(
            (geometry.x0, geometry.y0),
            geometry.sma * 2 * beam_factor_rectangle,
            geometry.sma * 2 * beam_factor_rectangle)
        aper_rectangle.plot(
            ax=ax,
            color=beam_color_rectangle,
            ec=beam_edgecolor_rectangle,
            fc=beam_color_rectangle,
            fill=beam_fill_rectangle,
            linewidth=0.5,
            zorder=10000,
        )
    aper = EllipticalAperture(
        (geometry.x0, geometry.y0), geometry.sma,
         geometry.sma*(1 - geometry.eps),
         geometry.pa)
    aper.plot(
        ax,
        color=beam_color,
        fill=beam_fill,
        # linewidth=beam_linewidth,
        zorder=10000,
        )

def plot_channel(cube, chan, arcsecpix, velchans,
        vmax=None, vmin=None, cmap="inferno", units="Mass [Msun]",
        refpix=[0,0], markorigin=True, add_beam=False,
        bmin=None, bmaj=None, pabeam=None, return_fig_axs=False,):
    """
    Plots a channel map of a spectral cube

    Parameters
    ----------
    cube : numpy.ndarray()
        Spectral cube from which the channel is plotted
    chan : int
        Index of the channel to plot
    arcsecpix : float
        Arcseconds per pixel
    velchans : list or numpy.ndarray()
        Array with the velocities of each channel
    vmax : float, optional
        Maximum value of the colormap. If None (default), the maximum value of
        the channel is chosen.
    vmin : float, optional
        Minimum value of the colormap. If None (default), the minimum value of
        the channel is chosen.
    cmap : str, optional
        Label of the colormap, by default "inferno".
    units : str, optional
        Units of the values of the cube, by default "Mass [Msun]"
    refpix : list, optional
        Pixel of reference, by default [0,0]
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
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

    fig = plt.figure(figsize=(4,3.75))
    
    gs = GridSpec(1, 2,
        wspace=0.05, hspace=0.05,
        width_ratios=[1, 0.05],
        height_ratios=[1]
        )
    ax = plt.subplot(gs[0, 0])
    cbax = plt.subplot(gs[0, 1])
    
    vmax = vmax if vmax is not None else np.max(cube[chan])
    vmin = vmin if vmin is not None else np.min(cube[chan])
    norm = colors.Normalize(
        vmax=vmax,
        vmin=vmin
    )
    _, nys, nxs = np.shape(cube)
    
    extent = np.array([
        -(-0.5-refpix[0]),
        -(nxs-0.5-refpix[0]),
        (-0.5-refpix[1]),
        (nys-0.5-refpix[1]),
        ]) * arcsecpix
     
    im = ax.imshow(
        cube[chan],
        origin="lower",
        norm=norm,
        cmap=cmap,
        extent=extent,
    )
    if markorigin:
        ax.plot(0, 0, "w+")
    ax.set_aspect("equal")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.text(0.075, 0.9,
        f"v$_{{LSR}}$ = {velchans[chan]:.2f} km/s",
        color="w",
        transform=ax.transAxes,
        #fontsize=10,
        )

    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    cbar = plt.colorbar(im, label=units, cax=cbax)
    cbax.tick_params(
        axis="y", right=True, left=False,
        labelright=True, direction="in", color="w")
    
    ax.set_ylabel("Dec. [arcsec]")
    ax.set_xlabel("R.A. [arcsec]")

    if return_fig_axs:
        return fig, ax, cbax


def plot_channels(cube, arcsecpix, velchans, 
        ncol=4, nrow=4, figsize=None, wspace=0.05, hspace=0.0, vmax=None,
        vcenter=None, vmin=None, cmap="inferno", units="Mass [Msun]",
        xmajor_locator=1, xminor_locator=0.2, ymajor_locator=1,
        yminor_locator=0.2, refpix=[0,0], markorigin=True, add_beam=False,
        bmin=None, bmaj=None, pabeam=None, return_fig_axs=False):
    """
    Plots several channel map of a spectral cube.

    Parameters
    ----------
    cube : numpy.ndarray()
        Spectral cube from which the channel is plotted
    arcsecpix : float
        Arcseconds per pixel
    velchans : list or numpy.ndarray()
        Array with the velocities of each channel
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
    vcenter : float, optional
        Center value of the colormap. If None (default), the middle vale will be
        chosen.
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
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    return_fig_axs : bool, optional
        If True, returns the figure, axes of the channel map, and the axes the
        colorbar. If False, does not return anything.
        
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

    size_factor = 2.5
    beam_ax = (nrow-1) * ncol
    figsize = figsize if figsize is not None \
        else (ncol*size_factor, nrow*size_factor)
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrow, ncol+1, 
                  height_ratios=[1]*nrow,
                  width_ratios=[1]*ncol+[0.1])
    gs.update(wspace=wspace, hspace=hspace)
    
    axs = {}
    for n, (i,j) in enumerate(product(np.arange(nrow), np.arange(ncol))):
        axs[n] = plt.subplot(gs[i,j])
    cbax = plt.subplot(gs[:, ncol])
    
    nc, nys, nxs = np.shape(cube)
    nchanscube = nrow*ncol
    chans_plot = nc 
    selint = int(chans_plot/nchanscube)
    initchan = selint 
    vmin = vmin if vmin is not None else np.min(cube)
    vmax = vmax if vmax is not None else np.max(cube)
    vcenter = vcenter if vcenter is not None else (vmax + vmin) / 2.
    norm = colors.TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
        
    iter_grid = [
        i for i in product(
            [i for i in range(nrow)],
            [j for j in range(ncol)]
        )
    ]
    
    extent = np.array(
        [-(-0.5-refpix[0]),
         -(nxs-0.5-refpix[0]),
         (-0.5-refpix[1]),
         (nys-0.5-refpix[1]),
        ]) * arcsecpix
     
    for chan, (i, j) in enumerate(iter_grid):
        data = cube[initchan::selint][chan]
        im = axs[chan].imshow(
            data,
            origin="lower",
            extent=extent,
            norm=norm,
            cmap="inferno",
        )   
        if markorigin:
            axs[chan].plot(0, 0, "w+")
        axs[chan].set_aspect("equal")
       
        axs[chan].text(0.05, 0.9,
                s=f"{velchans[initchan::selint][chan]:.2f} km/s",
                color="w",
                transform=axs[chan].transAxes,
                fontsize=10)
        axs[chan].tick_params(
            which="both",
            top=True,
            right=True,
            direction="in",
            color="w",
        )
        axs[chan].xaxis.set_major_locator(
            ticker.MultipleLocator(xmajor_locator))    
        axs[chan].yaxis.set_major_locator(
            ticker.MultipleLocator(ymajor_locator))    
        axs[chan].xaxis.set_minor_locator(
            ticker.MultipleLocator(xminor_locator))    
        axs[chan].yaxis.set_minor_locator(
            ticker.MultipleLocator(yminor_locator))    
        if (j > 0) and (i < nrow-1):
            axs[chan].set_xticklabels([])
            axs[chan].set_yticklabels([])
        if (i == (nrow-1)) and (j > 0):
            axs[chan].set_yticklabels([])
        if (i < (nrow-1)) and (j == 0):
            axs[chan].set_xticklabels([])

    if add_beam:
        add_beam_to_ax(axs[beam_ax], bmin, bmaj, pabeam, extent,)

    cbar = plt.colorbar(
           cm.ScalarMappable(
               norm=norm,
               cmap=cmap,
           ),
           cax=cbax,
           orientation="vertical",
        )
    cbax.tick_params(
        axis="y", right=True, left=False,
        labelright=True, direction="in", color="w"
        )
    cbax.set_ylabel(units)
    if return_fig_axs:
        return fig, axs, cbax


def plotpv(pvimage, rangex, chan_vels, ax=None,
        cbax=None, vmax=None, vcenter=None, vmin=None,
        cmap="nipy_spectral", interpolation="bilinear",
        cbarlabel="Intensity [Jy/beam]", return_fig_axs=False):
    """
    Plots the Position-Velocity diagram
    
    Parameters
    ----------
    pvimage : numpy.ndarray
        Position velocity diagram to plot.
    rangex : list or numpy.ndarray
        2 element list or numpy.ndarray corresponding to the physical
        coordinates of the boundaries of the image in the spatial direction.
    chan_vels : list or numpy.ndarray 
        list or numpy.ndarray with the velocities corresponding to the channels
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes` instance in which the position velodity diagram is drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    cmap : str, optional
        Label of the colorbar
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the Position
        Velocity diagram
    return_fig_axs : bool, optional
        If True, returns the figure, axes of the channel map, and the axes the
        colorbar. If False, does not return anything.
 
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
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1,0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass
    vmax = vmax if vmax is not None else np.max(pvimage[~np.isnan(pvimage)])
    vmin = vmin if vmin is not None else np.min(pvimage[~np.isnan(pvimage)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    chanwidth = chan_vels[1] - chan_vels[0]
    im = ax.imshow(
        pvimage,
        origin="lower",
        extent=[
            rangex[0], rangex[1],
            chan_vels[0]-chanwidth/2,
            chan_vels[-1]+chanwidth/2
            ],
        norm=norm,
        cmap=cmap,
        interpolation=interpolation,
        )
    ax.set_aspect(np.abs(rangex[0]-rangex[-1]) / np.abs(chan_vels[0]-chan_vels[-1]) )
    ax.set_ylabel("Velocity [km/s]")
    ax.set_xlabel("Distance [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.set_xlabel(cbarlabel)
    cbax.xaxis.set_label_position("top")
    if return_fig_axs:
        return fig, ax, cbax


def plotsumint(sumint, ax=None, cbax=None, extent=None,
        vmax=None, vcenter=None, vmin=None, cmap="inferno",
        interpolation="bilinear", cbarlabel="Intensity", add_beam=False,
        bmin=None, bmaj=None, pabeam=None, markorigin=True,
        return_fig_axs=False):
    """
    Plots the sumation of all the pixels along the velocity axis
    
    Parameters
    ----------
    sumint : numpy.ndarray
        Image of the sumation of all the pixels along the velocty axis.
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the position velodity diagram is
        drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    cmap : str, optional
        Label of the colorbar
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the plot
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    return_fig_axs : bool, optional
        If True, returns the figure, axes of the channel map, and the axes the
        colorbar. If False, does not return anything.
 
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
 
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass
    vmax = vmax if vmax is not None else np.max(sumint[~np.isnan(sumint)])
    vmin = vmin if vmin is not None else np.min(sumint[~np.isnan(sumint)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        sumint,
        origin="lower",
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    if markorigin:
        ax.plot(0, 0, "w+")
    ax.set_aspect("equal")

    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    plt.colorbar(im, cax=cbax, orientation="horizontal")
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.set_xlabel(rf"$\sum\mathrm{{{cbarlabel}}}_i$")
    cbax.xaxis.set_label_position("top")
    if return_fig_axs:
        return fig, ax, cbax


def plotmom0(mom0, ax=None, cbax=None, extent=None,
        vmax=None, vcenter=None, vmin=None, cmap="inferno",
        interpolation="bilinear", cbarlabel="Moment 0 [Jy/beam km/s]",
        add_beam=False, bmin=None, bmaj=None, pabeam=None, markorigin=True,
        return_fig_axs=False):
    """
    Plots the moment 0 
    
    Parameters
    ----------
    mom0 : numpy.ndarray
        Image of the moment 0.
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    cmap : str, optional
        Label of the colorbar
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the plot
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    return_fig_axs : bool, optional
        If True, returns the figure, axes of the channel map, and the axes the
        colorbar. If False, does not return anything.
 
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
 
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    vmax = vmax if vmax is not None else np.max(mom0[~np.isnan(mom0)])
    vmin = vmin if vmin is not None else np.min(mom0[~np.isnan(mom0)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        mom0,
        origin="lower",
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.set_aspect("equal")
    if markorigin:
        ax.plot(0, 0, "w+")
    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")
    if return_fig_axs:
        return fig, ax, cbax


def plotmom1(
        mom1, ax=None, cbax=None, extent=None, vmin=None, vmax=None,
        vcenter=None, extend_cbar="max", bg="black", cmap_ref='jet_r',
        interpolation="bilinear", cbarlabel="Moment 1 [km/s]", add_beam=False,
        bmin=None, bmaj=None, pabeam=None, markorigin=True,
        return_fig_axs_velcmap=False,):
    """
    Plots the moment 1 
    
    Parameters
    ----------
    mom1 : numpy.ndarray
        Image of the moment 1.
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    extend_cbar : str,
        Extremum to extend the colorbar: "max" or "min"
    bg : str
        Color of the background image.
    cmap_ref : str, optional
        Label of the colorbar used for reference.
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the plot
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    return_fig_axs_velcmap : bool, optional
        If True, returns the figure, axes of the channel map, the axes the
        colorbar, and the colormap created to plot the moment 1. If False, does
        not return anything.
        
    Returns:
    -------- 
    fig : matplotlib.figure.Figure
        Figure instance, only return_fig_axs_velcmap=True
    ax : matplotlib.axes.Axes
        Axes of the channel map, only return_fig_axs_velcmap=True
    cbax : tuple of matplotlib.axes.Axes
        Axes of the channel map and the colorbar, only returns if
        return_fig_axs_velcmap=True.
    velcmap : matplotlib.colors.ListedColormap
        Colormap created to plot the moment 1
    """
 
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    if type(cmap_ref) is str:
        cmap = colormaps[cmap_ref]
    else:
        cmap = cmap_ref
    velcolors = cmap(np.linspace(0, 1, 256))
    if bg == "black":
        bgcolor = np.array([0/256, 0/256, 0/256, 1])
    elif bg == "white":
        bgcolor = np.array([256/256, 256/256, 256/256, 1])
    velcolors[:1, :] = bgcolor

    velcmap = ListedColormap(velcolors)

    if extend_cbar == "max":
        velcmap = ListedColormap(velcolors[::-1])

    vmin = vmin if vmin is not None else np.min(
        mom1[(~np.isnan(mom1)) & (~np.isclose(0,mom1,atol=1))]
        )
    vmax = vmax if vmax is not None else np.max(
        mom1[(~np.isnan(mom1)) & (~np.isclose(0,mom1,atol=1))]
        )
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vcenter=vcenter, vmax=vmax, vmin=vmin)
    im = ax.imshow(
        mom1,
        origin="lower",
        extent=extent,
        norm=norm,
        cmap=velcmap,
        interpolation=interpolation,
    )
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.set_aspect("equal")
    if markorigin:
        ax.plot(0, 0, "w+")
    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    plt.colorbar(im, cax=cbax, orientation="horizontal",
                 extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")
    if return_fig_axs_velcmap:
        return fig, ax, cbax, velcmap


def plotmom2(mom2, ax=None, cbax=None, extent=None,
        vmin=None, vmax=None, vcenter=None, extend_cbar="max", bg="black",
        cmap_ref='jet_r', cbarlabel="Moment 2 [km/s]", interpolation=None,
        add_beam=False, bmin=None, bmaj=None, pabeam=None, markorigin=True,
        return_fig_axs_velcmap=False):
    """
    Plots the moment 2 
    
    Parameters
    ----------
    mom2 : numpy.ndarray
        Image of the moment 2.
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    extend_cbar : str,
        Extremum to extend the colorbar: "max" or "min"
    bg : str
        Color of the background image.
    cmap_ref : str, optional
        Label of the colorbar used for reference.
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the plot
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]
    return_fig_axs_velcmap : bool, optional
        If True, returns the figure, axes of the channel map, the axes the
        colorbar, and the colormap created to plot the moment 2. If False, does
        not return anything.

    Returns:
    -------- 
    fig : matplotlib.figure.Figure
        Figure instance, only return_fig_axs_velcmap=True
    ax : matplotlib.axes.Axes
        Axes of the channel map, only return_fig_axs_velcmap=True
    cbax : tuple of matplotlib.axes.Axes
        Axes of the channel map and the colorbar, only returns if
        return_fig_axs_velcmap=True.
    velcmap : matplotlib.colors.ListedColormap
        Colormap created to plot the moment 2
    """
 
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    if type(cmap_ref) is str:
        cmap = colormaps[cmap_ref]
    else:
        cmap = cmap_ref
    velcolors = cmap(np.linspace(0, 1, 256))
    if bg == "black":
        bgcolor = np.array([0/256, 0/256, 0/256, 1])
    elif bg == "white":
        bgcolor = np.array([256/256, 256/256, 256/256, 1])
    velcolors[:1, :] = bgcolor
    velcmap = ListedColormap(velcolors)
    if extend_cbar == "max":
        velcmap = ListedColormap(velcolors[::-1])
    vmin = vmin if vmin is not None else np.min(
        mom2[(~np.isnan(mom2)) & (~np.isclose(0,mom2,atol=1))]
        )
    vmax = vmax if vmax is not None else np.max(
        mom2[(~np.isnan(mom2)) & (~np.isclose(0,mom2,atol=1))]
        )
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vcenter=vcenter, vmax=vmax, vmin=vmin)
    im = ax.imshow(
        mom2,
        origin="lower",
        extent=extent,
        norm=norm,
        cmap=velcmap,
        interpolation=interpolation,
    )
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.set_aspect("equal")
    if markorigin:
        ax.plot(0, 0, "w+")
    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    plt.colorbar(im, cax=cbax, orientation="horizontal",
                extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")
    if return_fig_axs_velcmap:
        return fig, ax, cbax, velcmap


def plotmom8(mom8, ax=None, cbax=None, extent=None,
        vmax=None, vcenter=None, vmin=None, cmap="inferno",
        interpolation="bilinear", cbarlabel="Moment 8", add_beam=False,
        bmin=None, bmaj=None, pabeam=None, markorigin=True,
        return_fig_axs=False):
    """
    Plots the maximum value of the pixels along the velocity axis
    
    Parameters
    ----------
    mom8 : numpy.ndarray
        Image of the maximum pixel along the velocty axis.
    ax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : matplotlib.axes.Axes, optional
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : float, optional
        Maximum value of the colorbar
    vcenter : float, optional
        Center value of the colorbar
    vmin : float, optional
        Minimum value of the colorbar
    cmap : str, optional
        Label of the colorbar
    interpolation : str, optional
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : str, optional
        String with information on the quantity represented in the plot
    markorigin : boolean, optional
        If True, a marker will be plot at [0,0]. Default True.
    add_beam : bolean, optional
        If True, plots a ellipse of the beam size in the left bottom corner.
    bmin : float, optional
        Beam minor axis [arcsec]
    bmaj : float, optional
        Beam major axis [arcsec]
    pabeam : float
        Beam position angle [degrees]

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
 
    if ax is None or cbax is None:
        fig = plt.figure(figsize=(5,5.5))
        gs = GridSpec(
            2, 1,
            height_ratios=[0.05, 1],
            width_ratios=[1],
            hspace=0.1,
            wspace=0.00,
        )
        ax = plt.subplot(gs[1, 0])
        cbax = plt.subplot(gs[0, 0])
    else:
        pass

    vmax = vmax if vmax is not None else np.max(mom8[~np.isnan(mom8)])
    vmin = vmin if vmin is not None else np.min(mom8[~np.isnan(mom8)])
    vcenter = vcenter if vcenter is not None else (vmax - vmin) / 2 + vmin
    norm = TwoSlopeNorm(vmax=vmax, vcenter=vcenter, vmin=vmin)
    im = ax.imshow(
        mom8,
        origin="lower",
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        extent=extent,
        )
    if extent is None:
        ax.set_ylabel("Dec. [pixel]")
        ax.set_xlabel("R.A. [pixel]")
    else:
        ax.set_ylabel("Dec. [arcsec]")
        ax.set_xlabel("R.A. [arcsec]")
    ax.minorticks_on()
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        color="w",
    )
    ax.set_aspect("equal")
    if markorigin:
        ax.plot(0, 0, "w+")
    if add_beam:
        add_beam_to_ax(ax, bmin, bmaj, pabeam, extent,)

    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")
    if return_fig_axs:
        return fig, ax, cbax