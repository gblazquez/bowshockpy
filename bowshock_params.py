import numpy as np

from datetime import datetime

"""
Use this file to define all the the parameters needed to run bowshock.py.

For more information about the physical meaning of some of these parameters, see
Tabone et al. [2018)
"""

"""
MODEL OUTPUTS
"""
# Name of the model folder
modelname = f"L1448_{datetime.now().strftime('%y%m%d_%H%M%S')}"

# Plot 2D bowshock model [True/False]
bs2Dplot = True

# List of the output cube quantites and operations performed onto the cubes.
# The string should follow this format: {quantity}_{operations}
# Available quantities:
#     m: Mass
#     NCO: CO column density
#     tau: Opacity
#     I: Intensity
#     Ithin: Intensity taking into account the optically thin approximation.
#
# Available operations:
#     s: Add source
#     r: Rotate
#     n: Add gaussian noise
#     c: Convolve with beam
# Operations can be concatenated. Some examples of elements that could be
# included in the list are:
#    "m": Compute the masses in every pixel and channel.
#    "tau_r": Compute the opacities in every pixel and channel, and rotate the model.
#    "I_rnc": Compute the intensities in every pixel and channel, rotate, add
#    noise, and convolve.
# The list can be left empty if no cube is desired
# Example of outcubes:
# outcubes = ["m", "m_r", "I_rc", "tau_rc", "NCO_rc", "Ithin_rc"]
outcubes = ["I_srnc", ]


# Verbose messages about the computation? [True/False]
verbose = True

"""
OBSERVER PARAMETERS
"""

# Source distance to the observer [pc]
distpc = 300

# Jet inclination angle with respect to the line of sight [degrees]
i = 180-70

# Systemic velocity of the source [km/s]
vsys = + 5

# Source coordinates [deg, deg]
ra_source_deg, dec_source_deg = 51.41198333, 30.73479833

# Reference pixel [[int, int] or None]
# The x'y' pixel coordinates of the source. The x'y' is the plane of the sky,
# being x' the projection of the symmetry axis z onto the plane of the sky. In
# the model output cube, previous to any rotation defined by the PA angle of the
# x' axis, x' is the abscisa axis and y' is the ordinate axis.
refpix = None

# Add source to the image at the reference pixel? [True/False]
add_source = True

# Noise of the map
maxcube2noise = 10


"""
BOWSHOCK PARAMETERS
"""

# Number of bowshocks to model
nbowshocks = 5

# Excitation temperature [K]
Tex = 100

# Background temperature [K]
Tbg = 2.7

# Mean molecular mass per H molecule
muH2 = 2.8

# CO abundance
XCO = 8.5 * 10**(-5)


# The individual bowshock parameters must end in _{bowshock_number}. For example, the jet
# velocity for the third bowshock is vj_3

"""
bowshock 1 [red]
"""
# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_1 = 0

# Characteristic length scale [arcsec]
L0_1 = 0.5

# Distance between the working surface and the source [arcsec]
zj_1 = 1. / np.sin(i*np.pi/180)

# Jet velocity
vj_1 = (53-vsys) / (-np.cos(i*np.pi/180))

# Ambient (or wind) velocity [km/s]
vw_1 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_1 = 20

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs_1 = 1 

# Total mass of the bowshock [solar masses]
mass_1 = 0.00031

"""
bowshock 2 [red]
"""
# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_2 = 0

# Characteristic length scale [arcsec]
L0_2 = 0.5

# Distance between the working surface and the source [arcsec]
zj_2 = 3.38 / np.sin(i*np.pi/180)

# Jet velocity
vj_2 = (56-vsys) / (-np.cos(i*np.pi/180))

# Ambient (or wind) velocity [km/s]
vw_2 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_2 = 20

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs_2 = 1

# Total mass of the bowshock [solar masses]
mass_2 = 0.00031 * 1.25

"""
bowshock 3 [red]
"""
# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_3 = 0

# Characteristic length scale [arcsec]
L0_3 = 0.5

# Distance between the working surface and the source [arcsec]
zj_3 = 5.6 / np.sin(i*np.pi/180)

# Jet velocity
vj_3 = (60-vsys) / (-np.cos(i*np.pi/180))

# Ambient (or wind) velocity [km/s]
vw_3 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_3 = 20

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs_3 = 1

# Total mass of the bowshock [solar masses]
mass_3 = 0.00031 * 1.5


"""
bowshock 4 [red]
"""
# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_4 = 0

# Characteristic length scale [arcsec]
L0_4 = 0.5

# Distance between the working surface and the source [arcsec]
zj_4 = 8 / np.sin(i*np.pi/180)

# Jet velocity
vj_4 = (63-vsys) / (-np.cos(i*np.pi/180))

# Ambient (or wind) velocity [km/s]
vw_4 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_4 = 20

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs_4 = 1 

# Total mass of the bowshock [solar masses]
mass_4 = 0.00031 * 1.75


"""
bowshock 5 [red]
"""
# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj_5 = 0

# Characteristic length scale [arcsec]
L0_5 = 0.5

# Distance between the working surface and the source [arcsec]
zj_5 = 10.4 / np.sin(i*np.pi/180)

# Jet velocity
vj_5 = (65-vsys) / (-np.cos(i*np.pi/180))

# Ambient (or wind) velocity [km/s]
vw_5 = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0_5 = 20

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs_5 = 1

# Total mass of the bowshock [solar masses]
mass_5 = 0.00031 *2




"""
SPECTRAL CUBE PARAMETERS
"""

# Number of points to model
nzs = 500

# Number of azimuthal angle phi to calculate the bowshock solution
nphis = 500

# Number of spectral channel maps
nc = 50

# Central velocity of the first channel map [km/s]
vch0 = 0

# Central velocity of the last channel map [km/s]
vchf = +100

# Number of pixels in the x and y axes
nxs, nys = (250, 250)

# Reference pixel where the physical center (the source) is found
refpix = [50, 125]

# Size of the channel maps along the x axis [arcsec]
xpmax = 15

# Position angle of the jet axis [degrees]
pa = 161.5

# Beam size [arcsec]
bmaj, bmin = (0.420, 0.287)

# Beam position angle [degrees]
pabeam = -17.2 

# Thermal+turbulent line-of-sight velocity dispersion [km/s] If
# thermal+turbulent line-of-sight velocity dispersion is smaller than the
# instrumental spectral resolution, vt should be the spectral resolution.
# It can be also set to a integer times the channel width (e.g., "2xchannel")
vt = "2xchannel"

# Cloud in Cell interpolation? [True/False]
CIC = True

# Neighbour channel maps around a given channel map with vch will stop being
# populated when their difference in velocity with respect to vch is higher than
# this factor times vt. The lower the factor, the quicker will be the code, but
# the total mass will be underestimated. If vt is not None, compare the total
# mass of the output cube with the 'mass' parameter that the user has defined
tolfactor_vt = 3
