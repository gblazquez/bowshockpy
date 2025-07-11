import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib import colormaps
from matplotlib import colors

from astropy.io import fits
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve

import os

from datetime import datetime

from bowshockpy._header_default import hdr_str_default

def print_example(nexample):
    """
    Prints one of the available examples of input file to run bowshockpy.

    Parameters:
    -----------
    nexample : str or int
        Number of the example to print. There are 4 examples:
            - Example 1: A redshfted bowshock
            - Example 2: A blueshifted bowshock
            - Example 3: A side-on bowshock
            - Example 4: Several bowshocks in one cube
    """
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    with open(f"example{nexample}.py", "w") as wr:
        with open(ROOT_DIR+f"/inputfiles/example{nexample}.py", "r") as re:
            for line in re:
                wr.write(line)

def list2str(a, precision=2):
    _list = [float(f'{i:.{precision}f}') for i in a]
    _str = str(_list) if len(_list)>1 else str(_list[0])
    return _str

def progressbar_bowshock(
        iteration, total, timelapsed, intervaltime,
        decimals = 1, length = 100, fill = '─', printend = "\r"
        ):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + ')' + ' ' * (length - filledLength)
    print(f'  0{bar}{percent}% | {timelapsed:.0f}/{intervaltime*total:.0f}s', end = printend)
    if iteration == total:
        print()

def create_hdr(**kwargs):
    hdr = fits.Header.fromstring(hdr_str_default)
    if len(kwargs) != 0:
         for kwarg in kwargs:
             hdr[kwarg] = kwargs[kwarg]
    return hdr

def make_folder(foldername=None):
    if not os.path.exists(foldername):
        os.makedirs(foldername)

def write_log(path, mode, pars, printtime=False):
    with open(path, mode) as f:
        if printtime:
            time = datetime.now().strftime("%d/%m/%y %H:%M:%S")
            f.writelines(f"###################\n# {time}\n###################\n")
        f.writelines("pars = {\n")
        for par in pars:
            f.writelines(f'    "{par}": "{pars[par]}",\n')
        f.writelines("}\n\n")

def mb_sa_gaussian_f(maja, mina):
    """
    Computes the solid angle of a gaussian main beam and θmaj and θmin as the
    half-power beam widths

    Parameters:
    -----------
    maja : astropy.units.Quantity
        Beam major axis in degrees or radians
    mina : astropy.units.Quantity
        Beam minor axis  in degrees or radians

    Returns:
    --------
    omega_M : astropy.units.sr
        Beam solid angle in stereoradians
    """
    omega_M = np.pi * maja * mina / (4 * np.log(2))
    return omega_M.to(u.sr)

def gaussconvolve(data, x_FWHM, y_FWHM, pa, return_kernel=False):
    """
    Convolves data with a Gaussian kernel

    Parameters:
    -----------
    data : numpy.ndarray
        Data to convolve
    x_FWHM : float
        Full width half maximum of the Gaussian kernel for the x direction
    y_FWHM : float
        Full width half maximum of the Gaussian kernel for the y direction
    pa : float
        Position angle in degrees
    return_kernel : optional, bool
        Whether to return the kernel or not

    Returns:
    --------
    data_conv : numpy.ndarray
        Convolved data
    kernel : numpy.ndarray
        Image of the Gaussian kernel. Is returned only if  return_kernel = True
        """
    x_stddev = x_FWHM / (2 * np.sqrt(2 * np.log(2)))
    y_stddev = y_FWHM / (2 * np.sqrt(2 * np.log(2)))
    # Gausskernel 0 and 1 entries are the FWHM, the third the PA
    kernel = Gaussian2DKernel(
        x_stddev=x_stddev,
        y_stddev=y_stddev,
        theta=pa*np.pi/180)
    data_conv = convolve(data, kernel)
    if return_kernel:
        return data_conv, kernel
    else:
        return data_conv

def get_color(vel_range, vel, cmap, norm="linear"):
    """
    Gets the color that corresponds in a colormap linearly interpolated taking
    into account the values at the limits.

    Parameters:
    -----------
    vel_range : list
        List with 2 elements defining the range of values to be represented by
        the colors
    vel : float
        Value to get the corresponding color from
    cmap : str
        Colormap label
    norm : optional, str
        Set "linear" for a linear scale, "log" for log scale.
    """
    cmapp = colormaps.get_cmap(cmap)
    if norm == "linear":
        norm = colors.Normalize(vmin=vel_range[0], vmax=vel_range[-1])
    elif norm == "log":
        norm = colors.LogNorm(vmin=vel_range[0], vmax=vel_range[-1])
    rgba = cmapp(norm(vel))
    color = colors.to_hex(rgba)
    return color

def plotpv(pvimage, rangex, chan_vels, ax=None, cbax=None,
        vmax=None, vcenter=None, vmin=None,
        cmap="nipy_spectral", interpolation="bilinear",
        cbarlabel="Intensity [Jy/beam]",):
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
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes` instance in which the position velodity diagram is drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    cmap : optional, str
        Label of the colorbar
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the Position
        Velocity diagram
    """
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5))
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

def plotsumint(sumint, ax=None, cbax=None, extent=None,
               vmax=None, vcenter=None, vmin=None,
               cmap="inferno", interpolation="bilinear",
               cbarlabel="Intensity",):
    """
    Plots the sumation of all the pixels along the velocity axis
    
    Parameters
    ----------
    sumint : numpy.ndarray
        Image of the sumation of all the pixels along the velocty axis.
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the position velodity diagram is
        drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    cmap : optional, str
        Label of the colorbar
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the plot
    """
 
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
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
    ax.set_aspect("equal")
    plt.colorbar(im, cax=cbax, orientation="horizontal")
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.set_xlabel(rf"$\sum\mathrm{{{cbarlabel}}}_i$")
    cbax.xaxis.set_label_position("top")

def plotmom0(mom0, ax=None, cbax=None, extent=None,
            vmax=None, vcenter=None, vmin=None,
            cmap="inferno",
            interpolation="bilinear", cbarlabel="Moment 0 [Jy/beam km/s]",):
    """
    Plots the moment 0 
    
    Parameters
    ----------
    mom0 : numpy.ndarray
        Image of the moment 0.
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    cmap : optional, str
        Label of the colorbar
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the plot
    """
 
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
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
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")


def plotmom1(mom1, ax=None, cbax=None, extent=None,
              vmin=None, vmax=None, vcenter=None,
              extend_cbar="max", return_velcmap=False,
              bg="black", cmap_ref='jet_r',
              interpolation="bilinear", cbarlabel="Moment 1 [km/s]"):
    """
    Plots the moment 1 
    
    Parameters
    ----------
    mom1 : numpy.ndarray
        Image of the moment 1.
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    extend_cbar : str,
        Extremum to extend the colorbar: "max" or "min"
    bg : str
        Color of the background image.
    cmap_ref : optional, str
        Label of the colorbar used for reference.
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the plot

    Returns:
    -------- 
    velcmap : matplotlib.colors.ListedColormap
        Colormap created to plot the moment 1
    """
 
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
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
    plt.colorbar(im, cax=cbax, orientation="horizontal",
                 extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")

    if return_velcmap:
        return velcmap
    else:
        pass

def plotmom2(mom2, ax=None, cbax=None, extent=None,
              vmin=None, vmax=None, vcenter=None,
              extend_cbar="max", return_velcmap=False,
              bg="black", cmap_ref='jet_r', cbarlabel="Moment 2 [km$^2$/s$^2$]",
              interpolation=None):
    """
    Plots the moment 2 
    
    Parameters
    ----------
    mom2 : numpy.ndarray
        Image of the moment 2.
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    extend_cbar : str,
        Extremum to extend the colorbar: "max" or "min"
    bg : str
        Color of the background image.
    cmap_ref : optional, str
        Label of the colorbar used for reference.
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the plot

    Returns:
    -------- 
    velcmap : matplotlib.colors.ListedColormap
        Colormap created to plot the moment 2
    """
 
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
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
    plt.colorbar(im, cax=cbax, orientation="horizontal",
                extend=extend_cbar, label=cbarlabel)
    cbax.tick_params(axis="x", top=True, bottom=False, labelbottom=False,
                     labeltop=True, direction="in")
    cbax.xaxis.set_label_position("top")
    if return_velcmap:
        return velcmap
    else:
        pass

def plotmom8(mom8, ax=None, cbax=None, extent=None,
            vmax=None, vcenter=None, vmin=None, cmap="inferno",
            interpolation="bilinear", cbarlabel="Moment 8"):
    """
    Plots the maximum value of the pixels along the velocity axis
    
    Parameters
    ----------
    mom8 : numpy.ndarray
        Image of the maximum pixel along the velocty axis.
    ax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the position velodity diagram is drawn.
    cbax : optional, matplotlib.axes.Axes
        The matplotlib.axes.Axes instance in which the color bar is drawn.
    extent : list
        Physical coordinates of the boundaries of the image.
    vmax : optional, float
        Maximum value of the colorbar
    vcenter : optional, float
        Center value of the colorbar
    vmin : optional, float
        Minimum value of the colorbar
    cmap : optional, str
        Label of the colorbar
    interpolation : optional, str
        Interpolation to pass to matplotlib.pyplot.imshow
    cbarlabel : optional, str
        String with information on the quantity represented in the plot
    """
 
    if ax is None or cbax is None:
        plt.figure(figsize=(5,5.5))
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
    plt.colorbar(im, cax=cbax, orientation="horizontal", label=cbarlabel)
    cbax.tick_params(
        axis="x", top=True, bottom=False,
        labelbottom=False, labeltop=True,
        direction="in",
        )
    cbax.xaxis.set_label_position("top")

