import numpy as np

def sumint(cube, chan_range):
    """
    chan0                       chanf
      |   data  |  data  |  data  |
      0         1        2        3
    """
    return np.sum(cube[np.arange(*chan_range),:,:], axis=0)

def mom0(cube, chan_vels, chan_range):
    dv = np.abs(chan_vels[1]-chan_vels[0])
    return dv * sumint(cube, chan_range)

def sumIixvi(cube, chan_vels, chan_range, exp=1):
    Iixvi = np.array([cube[i, :, :]*chan_vels[i]**exp
                      for i in np.arange(chan_range[0], chan_range[1])])
    return np.sum(Iixvi, axis=0)

def mom1(cube, chan_vels, chan_range):
    return sumIixvi(cube, chan_vels, chan_range) / sumint(cube, chan_range)
    # or
    # dv = np.abs(chan_vels[1]-chan_vels[0])
    # return sumIixvi(cube, chan_vels, chan_range) * dv / mom0(cube, chan_vels, chan_range)

def mom2(cube, chan_vels, chan_range):
    # vimmom12 = np.array([(chan_vels[i] - mom1(cube,chan_vels,chan_range))**2
    #                      for i in np.arange(chan_range[0], chan_range[1])])
    # Iixdisp2 = np.array([cube[i, :, :]*(chan_vels[i] - mom1(cube,chan_vels,chan_range))**2
    #               for i in np.arange(chan_range[0], chan_range[1])])
    # sum_Iixdisp2 = np.sum(Iixdisp2, axis=0)
    # return np.sqrt(sum_Iixdisp2 / mom0(cube,chan_vels,chan_range))
    disp = np.sqrt(sumIixvi(cube, chan_vels, chan_range, exp=2)
                   /sumint(cube, chan_range) - mom1(cube, chan_vels, chan_range)**2)

    return disp

def mom8(cube, chan_range):
    return np.max(cube[chan_range[0]:chan_range[1], :, :], axis=0)

def pv(cube, xpv, halfwidth, axis=1):
    pixarray = np.array([*np.arange(xpv-halfwidth, xpv),
                         xpv,
                         *np.arange(xpv+1, xpv+halfwidth+1)])
    selected_data = cube[:, pixarray, :]
    return np.mean(selected_data, axis=axis)
