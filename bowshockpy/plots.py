import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib import colors

import bowshockpy.utils as ut



class BowshockModelPlot():
    """
    Figure including the main parameters of the bowshock model, its morphology
    and kinematics, and the distribution of the surface density

    Parameters:
    -----------
    bsm : `~bowshockpy.models.NarrowJet` class instance
        Instance of the model to plot
    modelname : optional, str
        Name of the model to include in the plot
    nzs : optional, int
        Number of points used to compute the model solutions
    figsize: optional, tuple
        Tuple passed to `matplotib.pyplot.figure` to define the dimensions of
        the figure
    narrows: optional, int
        Number of arrows to show in order to indicate the velocity at each
        symmetrical half of the model.        
    v_arrow_ref: optional, float
        Velocity in km/s to use as reference in the reference arrow

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
            figsize=(14,6), narrows=10, v_arrow_ref=100):
        self.mo = bsm
        self.modelname = modelname
        self.nzs = nzs
        self.nrs = nzs
        self.narrows = narrows 
        self.figsize = figsize
        self.v_arrow_ref = v_arrow_ref

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
        width_ratios = [1] * ncol
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
            $t_\mathrm{{j}} = {{{self.mo.tj_yr:.2f}}}$ yr
            $\rho_a = {{{self.mo.rhoa_gcm3*10**20:.2f}}}\times 10^{{-20}}$ g cm$^{{-3}}$
            $\dot{{m}}_0 = {{{self.mo.mp0_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            $\dot{{m}}_{{a,f}} = {{{self.mo.mpamb_f_solmassyr*10**6:.2f}}}\times10^{{-6}}$ M$_\odot$ yr$^{{-1}}$
            """

        self.axs["text"].set_axis_off()
        for n, line in enumerate(showtext.split("\n")):
             self.axs["text"].text(0, 1.1-0.09*n, line, fontsize=10,
                              transform=self.axs["text"].transAxes)

        """
        Deprojected shell Morph. and Kin.
        """
        for i, zarcsec in enumerate(self.zs_arcsec):
            c = ut.get_color(
                [self.minvs, self.maxvs],
                self.vs[i],
                "turbo_r",
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
                   cmap="turbo_r",
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
        Deprojected shell Morph. and Kin.
        """
        for i, zarcsec in enumerate(self.zs_arcsec):
            c = ut.get_color(
                [self.minvs, self.maxvs],
                self.vs[i],
                "viridis",
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
        self.minsurfdens = np.percentile(self.surfdenss_gcm2, 99) 
        self.maxsurfdens = np.percentile(self.surfdenss_gcm2, 85)
        cbar = plt.colorbar(
               cm.ScalarMappable(
                   norm=colors.LogNorm(
                       vmax=self.maxsurfdens,
                       vmin=self.minsurfdens,
#                       linthresh=self.maxsurfdens*0.99,
                       ),
                   cmap="viridis",
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

        self.cbaxs[1].tick_params(
            bottom=False, labelbottom=False,
            top=True, labeltop=True
        )
        self.cbaxs[1].set_xlabel(r"Surface density [g cm$^{-2}$]", )
        self.cbaxs[1].xaxis.set_label_position('top')

    def savefig(self, figname=None):
        """
        Saves the plot of the bowhsock model.       

        Parameters
        ----------
        figname : optional, str
            Full path name of the figure. If None, the the full path name will
            be models/{self.modelname}/bowshock_plot.pdf. If the folder tree does not exist, it will be created. 
        """
        if figname is None:
            ut.make_folder(f"models/{self.modelname}")
            figname = f"models/{self.modelname}/bowshock_plot.pdf"
        self.fig_model.savefig(f"{figname}", bbox_inches="tight")

