"""This module contains the functions that perform the statistical moments to
the spectral cubes"""

import numpy as np


def sumint(cube, chan_range):
    """
    Sums all pixels of a cube along the 0 axis.

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array.
    chan_range : list
        Two element list with the last and first channels used to compute the
        summation.

    Returns
    -------
    sumint_im : np.ndarray
        2 dimensional array
    """
    sumint_im = np.sum(cube[np.arange(*chan_range), :, :], axis=0)
    return sumint_im


def mom0(cube, chan_vels, chan_range):
    """
    Computes the moment of order 0 of a cube along the 0 axis.

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array.
    chan_vels: list or numpy.ndarray
        1-dimensional array of the velocity corresponding to the channels.
    chan_range : list
        Two element list with the last and first channels used to compute the
        moment.

    Return:
    -------
    mom0_im : numpy.ndarray
        Moment of order 0

    References:
    -----------
    https://www.aoc.nrao.edu/~kgolap/casa_trunk_docs/CasaRef/image.moments.html
    """
    dv = np.abs(chan_vels[1] - chan_vels[0])
    mom0_im = dv * sumint(cube, chan_range)
    return mom0_im


def sumIixvi(cube, chan_vels, chan_range, exp=1):
    """
    Computes the summation: sum (I_i * v_i)**exp

    Parameters
    ----------
    cube : numpy.array
        Input 3 dimensional array
    chan_vels : list or numpy.array
        1-dimensional array of the velocity corresponding to the channels.
    chan_range : list
        Two element list with the last and first channels used to compute the moment.
    exp : int, optional
        Exponent, by default 1

    Returns
    -------
    _type_
        _description_
    """
    Iixvi = np.array(
        [
            cube[i, :, :] * chan_vels[i] ** exp
            for i in np.arange(chan_range[0], chan_range[1])
        ]
    )
    return np.sum(Iixvi, axis=0)


def mom1(cube, chan_vels, chan_range):
    """
    Computes the moment of order 1 (intensity weighted mean velocity field) of
    a cube along the 0 axis.

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array.
    chan_vels: list or numpy.ndarray
        1-dimensional array of the velocity corresponding to the channels.
    chan_range : list
        Two element list with the last and first channels used to compute the
        moment.

    Return:
    -------
    mom1_im : numpy.ndarray
        Moment of 1 order

    References:
    -----------
    https://www.aoc.nrao.edu/~kgolap/casa_trunk_docs/CasaRef/image.moments.html

    But note that the equation for the moment 1 is wrong in this reference
    (units should be km/s, and not an adimensional value). The moment1 is  the
    intensity intensity weighted mean velocity = Sigma I_i v_i/Sigma I_i
    """
    mom1_im = sumIixvi(cube, chan_vels, chan_range) / sumint(cube, chan_range)
    return mom1_im


def mom2(cube, chan_vels, chan_range):
    """
    Computes the moment 2 (the intensity weighted dispersion) of a cube along
    the 0 axis

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array.
    chan_vels: list or numpy.ndarray
        1-dimensional array of the velocity corresponding to the channels.
    chan_range : list
        Two element list with the last and first channels used to compute
        moment.

    Return:
    -------
    disp : numpy.ndarray
        Moment of 2 order

    References:
    -----------
    https://www.aoc.nrao.edu/~kgolap/casa_trunk_docs/CasaRef/image.moments.html

    But note that the equation for the moment 2 is wrong in this reference
    (units should be km/s, no sqrt(km/s)). The moment2 is intensity weighted
    dispersion = [Sigma I_i v^2_i/Sigma I_i]**(1/2)

    """
    disp = np.sqrt(
        sumIixvi(cube, chan_vels, chan_range, exp=2) / sumint(cube, chan_range)
        - mom1(cube, chan_vels, chan_range) ** 2
    )
    return disp


def mom8(cube, chan_range):
    """
    Computes the moment the maximum value of a cube along the along the 0 axis

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array
    chan_range : list
        Two element list with the last and first channels used to compute the
        moment

    Return:
    -------
    mom8_im : numpy.ndarray
        Moment 8 (maximum value along the 0 axis)

    References:
    -----------
    https://www.aoc.nrao.edu/~kgolap/casa_trunk_docs/CasaRef/image.moments.html
    """
    mom8_im = np.max(cube[chan_range[0] : chan_range[1], :, :], axis=0)
    return mom8_im


def pv(cube, xpv, halfwidth, axis=1):
    """
    Computes the Position-Velocity diagram of a cube

    Parameters
    ----------
    cube : numpy.ndarray
        Input 3 dimensional array.
    xpv : int
        Pixel along which the PV-diagram is computed.
    halfwidth : int
        Number of pixels around xpv that will be taking into account to compute
        the PV-diagram.

    Return:
    -------
    pv_im : numpy.ndarray
        PV-diagram.
    """

    pixarray = np.array(
        [
            *np.arange(xpv - halfwidth, xpv),
            xpv,
            *np.arange(xpv + 1, xpv + halfwidth + 1),
        ]
    )
    selected_data = cube[:, pixarray, :]
    pv_im = np.mean(selected_data, axis=axis)
    return pv_im
