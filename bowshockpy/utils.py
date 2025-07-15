import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
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
    if norm == "linear" and customnorm is None:
        norm = colors.Normalize(vmin=vel_range[0], vmax=vel_range[-1])
    elif norm == "log" and customnorm is None:
        norm = colors.LogNorm(vmin=vel_range[0], vmax=vel_range[-1])
    elif customnorm is not None:
        norm = customnorm

    rgba = cmapp(norm(vel))
    color = colors.to_hex(rgba)
    return color

class VarsInParamFile():
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
