import numpy as np

import bowshockpy.plots as pl
from bowshockpy.models import BowshockModel


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
        self.i = i_deg * np.pi / 180
        self.pa_deg = pa_deg
        self.pa = pa_deg * np.pi / 180
        self.vsys = vsys
        # for param in model.__dict__:
        #     setattr(self, param, getattr(model, param))
        for kwarg in self.default_kwargs:
            kwarg_attr = (
                kwargs[kwarg] if kwarg in kwargs else self.default_kwargs[kwarg]
            )
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
        return self.rb(zb) * np.cos(phi) * np.cos(self.i) + zb * np.sin(self.i)

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
        return -self.rb(zb) * np.cos(phi) * np.sin(self.i) + zb * np.cos(self.i)

    def get_obsmodelplot(self, **kwargs):
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
            An instance of a class BowshockModelPlot, which contains
            information on the figure and the model data
        """
        modelplot = pl.BowshockObsModelPlot(self, **kwargs)
        return modelplot
