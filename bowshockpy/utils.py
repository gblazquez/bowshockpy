import os
from itertools import groupby

import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve
from matplotlib import colormaps, colors


def print_example(example):
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
    root_dir = os.path.dirname(os.path.abspath(__file__))
    with open(f"{example}", "w") as wr:
        with open(root_dir + f"/inputfiles/{example}", "r") as re:
            for line in re:
                wr.write(line)


def list2str(a, precision=2):
    """
    Converts a list to a str

    Parameters
    ----------
    a : list
        List to convert as string
    precision : int
        Number of decimals to display
    """
    _list = [float(f"{i:.{precision}f}") for i in a]
    _str = str(_list) if len(_list) > 1 else str(_list[0])
    return _str


def progressbar_bowshock(
    iteration,
    total,
    timelapsed,
    intervaltime,
    decimals=1,
    length=100,
    fill="─",
    printend="\r",
):
    """
    Bowshock-like progress bar

    Parameters
    ----------
    iteration : int
        Current iteraction
    total : int
        Total iteractions
    timelapsed : float
        Current time elapsed
    intervaltime : float
        Duration of an iteraction
    decimals : int
        Number of decimals to show
    length : float
        Length of the progress bar
    fill : str
        String to define the filled part of the progress bar
    printend : str
        End of the progress bar
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledlength = int(length * iteration // total)
    _bar = fill * filledlength + ")" + " " * (length - filledlength)
    print(
        f"  0{_bar}{percent}% | {timelapsed:.0f}/{intervaltime*total:.0f}s",
        end=printend,
    )
    if iteration == total:
        print()


def make_folder(foldername=None):
    """
    Makes a folder of the model

    Parameters
    ----------
    foldername : str
        Name of the folder of the model
    """
    if not os.path.exists(foldername):
        os.makedirs(foldername)


def mb_sa_gaussian_f(maja, mina):
    """
    Computes the solid angle of a Gaussian main beam

    Parameters:
    -----------
    maja : astropy.units.Quantity
        Beam major axis (FWHM) in degrees or radians
    mina : astropy.units.Quantity
        Beam minor axis (FWHM) in degrees or radians

    Returns:
    --------
    omega_m : astropy.units.sr
        Beam solid angle in stereoradians
    """
    omega_m = np.pi * maja * mina / (4 * np.log(2))
    return omega_m.to(u.sr)


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
        x_stddev=x_stddev, y_stddev=y_stddev, theta=pa * np.pi / 180
    )
    data_conv = convolve(data, kernel)
    if return_kernel:
        return data_conv, kernel
    return data_conv


def get_color(vel_range, vel, cmap, norm="linear", customnorm=None):
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
    customnorm : optional, str
        Custom norm from `matplotlib.colors`
    """
    cmapp = colormaps.get_cmap(cmap)
    if customnorm is not None:
        _norm = customnorm
    else:
        if norm == "log" and customnorm is None:
            _norm = colors.LogNorm(vmin=vel_range[0], vmax=vel_range[-1])
        else:
            _norm = colors.Normalize(vmin=vel_range[0], vmax=vel_range[-1])

    rgba = cmapp(_norm(vel))
    color = colors.to_hex(rgba)
    return color


class VarsInParamFile:
    """
    This class takes as attributes the keys and values of a dictionary

    Parameters
    ----------
    params : dict
        Input dictionary
    """

    def __init__(self, params):
        self.filename = params["__file__"]
        for key in params:
            if key.startswith("__") is False:
                setattr(self, key, params[key])


def allequal(inputlist):
    """
    Checks if all elements of an iterale object are equal

    Parameters
    ----------
    inputlist : list
        List object to check that all its elements are equal

    Returns
    -------
    boolean
        True if all elements are equal, False if they are not
    """
    if isinstance(inputlist[0], np.ndarray):
        _list = [list(i) for i in inputlist]
    else:
        _list = inputlist
    g = groupby(_list)
    return next(g, True) and not next(g, False)

