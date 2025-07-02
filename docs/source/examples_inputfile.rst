
Examples of input file
=======================

Examples of input files.


Example 1
-------------------

This example of input file generates one redshifted bowshock

.. code-block:: python 
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example1"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "intensity_opthin": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": ["convolve"],
      "opacity": [],
      "CO_column_density": ["convolve"],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 300
  
  vsys = + 0
  
  ra_source_deg, dec_source_deg = 51.41198333, 30.73479833
  
  
  """
  BOWSHOCK PARAMETERS
  """
  
  nbowshocks = 1
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = "3-2"
  
  XCO = 8.5 * 10**(-5)
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 180-45
  
  L0_1 = 0.7
  
  zj_1 = 3.5
  
  vj_1 = 73
  
  va_1 = 0
  
  v0_1 = 5
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00031 * 1.5
  
  pa_1 = -20
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 100
  
  nphis = 500
  
  nc = 50
  
  vch0 = 30
  
  vchf = 63
  
  nxs, nys = (200, 200)
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  CIC = True
  
  tolfactor_vt = 3
  
  refpix = [80, 30]
  
  parot = 0
  
  sigma_beforeconv = 0.1
  
  maxcube2noise = 0.07
  
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
  mom0values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom1values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom2values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom8values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  pvvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }

Example 2
-------------------

This example of input file generates two blueshifted bowshock

.. code-block:: python 
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example2"
  
  bs2Dplot = True
  
  outcubes = {
      "intensity": ["add_noise", "convolve", "moments_and_pv"],
      "intensity_opthin": ["add_noise", "convolve", "moments_and_pv"],
      "opacity": ["convolve"],
      "opacity": [],
      "CO_column_density": ["convolve"],
      "mass": [],
      }
  
  verbose = True
  
  """
  OBSERVER PARAMETERS
  """
  
  distpc = 300
  
  vsys = + 0
  
  ra_source_deg, dec_source_deg = 51.41198333, 30.73479833
  
  
  """
  BOWSHOCK PARAMETERS
  """
  
  nbowshocks = 2
  
  Tex = 100
  
  Tbg = 2.7
  
  muH2 = 2.8
  
  J = "3-2"
  
  XCO = 8.5 * 10**(-5)
  
  
  """
  bowshock 1 [redshifted]
  """
  
  i_1 = 180-45
  
  L0_1 = 0.7
  
  zj_1 = 3.5
  
  vj_1 = 73
  
  va_1 = 0
  
  v0_1 = 5
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00031 * 1.5
  
  pa_1 = -20
  
  """
  bowshock 2 [redshifted]
  """
  
  i_2 = 180-45
  
  L0_2 = 0.8
  
  zj_2 = 4.5
  
  vj_2 = 80
  
  va_2 = 0
  
  v0_2 = 7
  
  rbf_obs_2 = 1
  
  mass_2 = 0.00035 * 1.5
  
  pa_2 = -20
  
  
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 100
  
  nphis = 500
  
  nc = 50
  
  vch0 = 30
  
  vchf = 70
  
  nxs, nys = (200, 200)
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  CIC = True
  
  tolfactor_vt = 3
  
  refpix = [80, 30]
  
  parot = 0
  
  sigma_beforeconv = 0.1
  
  maxcube2noise = 0.07
  
  
  
  """
  MOMENTS AND PV PARAMETERS
  """
  
  savefits = True
  
  saveplot = True
  
  mom1clipping = "5xsigma"
  
  mom2clipping = "4xsigma"
  
  mom0values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom1values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom2values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  mom8values = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
  
  pvvalues = {
      "vmax": None,
      "vcenter": None,
      "vmin": None,
  }
