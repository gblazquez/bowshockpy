"""
Use this file to define all the the parameters needed to run bowshock.py.

For more information about the physical meaning of some of these parameters, see
Tabone et al. [2018) 
"""

"""
MODEL OUTPUTS
"""
# Name of the model folder
modelname = "bs1_prueba3"

# Plot 2D bowshock model [True/False]
bs2Dplot = True

# List of the output cube quantites and operations performed onto the cubes.
# The string should follow this format: {quantity}_{operations}
# Available quantities:
#     m: Mass
#     NCO: CO column density
#     tau: Opacity
#     I: Intensity
#     Ithin: Intensity taking into account the optically thin approximation. The
#
# Available operations:
#     s: Add source
#     r: Rotate
#     n: Add gaussian noise
#     c: Convolve with beam
# Operations can be concatenated, e.g, "I_rc" rotates and then convolves the
# Intensity cube. List can be left empty if no cube is desired
# Example:
# outcubes = ["m", "m_r", "I_rc", "tau_rc", "NCO_rc", "Ithin_rc"]
outcubes = ["I_rnc", ]


# Verbose messages about the computation? [True/False]
verbose = True 


"""
BOWSHOCK PARAMETERS
"""

# Jet radius. Set this parameter to zero, the channel maps generator
# are not yet generalized for jet radius>0 [arcsec]
rj = 0

# Characteristic length scale [arcsec]
L0 = 0.5

# Distance between the working surface and the source [arcsec]
zj = 7.31

# Jet velocity
vj = 108

# Ambient (or wind) velocity [km/s]
vw = 0

# Velocity at which the material is ejected from the internal working surface [km/s]
v0 = 17

# Final radius of the bowshock [arcsec]. Set None if you want to end the
# bowshock model at the theoretical final radius (see Tabone et al. 2018)
rbf_obs = None

# Total mass of the bowshock [solar masses]
mass = 0.00031

# Excitation temperature [K]
Tex = 100

# Background temperature [K]
Tbg = 2.7

# Mean molecular mass per H molecule
muH2 = 2.8

# CO abundance
XCO = 8.5 * 10**(-5)


"""
OBSERVER PARAMETERS
"""

# Source distance to the observer [pc]
distpc = 300

# Jet inclination angle with respect to the line of sight [degrees]
i = 20

# Systemic velocity of the source [km/s]
vsys = 0

# Source coordinates [deg, deg]
ra_source_deg, dec_source_deg = 52.26570167, 31.26771556

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
vchf = -150

# Number of pixels in the x and y axes 
nxs, nys = (250, 250)

# Reference pixel where the physical center (the source) is found
refpix = [-25, 125]

# Size of the channel maps along the x axis [arcsec]
xpmax = 3

# Position angle of the jet axis [degrees]
pa = 120

# Beam size [arcsec]
bmaj, bmin = (0.173, 0.091)

# Beam position angle [degrees] 
pabeam = 30

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