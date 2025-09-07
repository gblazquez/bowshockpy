"""This module contains the class that projects the bowshock model morphology
and kinematics"""

import numpy as np

import bowshockpy.plots as pl
from bowshockpy.models import BaseModel


class ObsModel(BaseModel):
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
        self.m = model
        self.i_deg = i_deg
        self.i = i_deg * np.pi / 180
        self.pa_deg = pa_deg
        self.pa = pa_deg * np.pi / 180
        self.vsys = vsys
        # for param in model.__dict__:
        #     setattr(self, param, getattr(model, param))
        # for kwarg in self.default_kwargs:
        #     kwarg_attr = (
        #         kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
        #     )
        #     setattr(self, kwarg, kwarg_attr)

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
        a = self.m.alpha(zb)
        return self.m.vtot(zb) * (
            np.cos(a) * np.cos(self.i) - np.sin(a) * np.cos(phi) * np.sin(self.i)
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
        return self.m.rb(zb) * np.cos(phi) * np.cos(self.i) + zb * np.sin(self.i)

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
        return self.m.rb(zb) * np.sin(phi)

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
        return -self.m.rb(zb) * np.cos(phi) * np.sin(self.i) + zb * np.cos(self.i)

    def get_obsmodelplot(
        self,
        modelname="none",
        nzs=150,
        nphis=150,
        figsize=(12, 6),
        linespacing=0.09,
        textbox_widthratio=0.8,
        cmap="turbo",
        minpointsize=0.1,
        maxpointsize=10,
        **kwargs,
    ):
        """
        Plot a figure including the main parameters of the bowshock model, its
        morphology and kinematics, and the distribution of the surface density

        Parameters
        -----------
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
        kwargs : optional
            Keyword arguments into `~bowshockpy.plot.BowshockObsModelPlot`

        Returns
        --------
        modelplot : `~bowshockpy.plot.BowshockObsModelPlot` class instance
            An instance of a class BowshockModelPlot, which contains
            information on the figure and the model data
        """
        modelplot = pl.BowshockObsModelPlot(
            self,
            modelname=modelname,
            nzs=nzs,
            nphis=nphis,
            figsize=figsize,
            linespacing=linespacing,
            textbox_widthratio=textbox_widthratio,
            cmap=cmap,
            minpointsize=minpointsize,
            maxpointsize=maxpointsize,
            **kwargs,
        )
        return modelplot
