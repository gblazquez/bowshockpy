
Examples of input file
=======================

Examples of input files.


Example 1
-------------------

This example of input file generates two blueshifted bowshock

.. code-block:: python
  
  """
  MODEL OUTPUTS
  """
  modelname = f"example3"
  
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
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
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
  bowshock 1 [blueshifted]
  """
  
  i_1 = 25
  
  L0_1 = 0.8
  
  zj_1 = 3.5
  
  vj_1 = 80
  
  va_1 = 0
  
  v0_1 = 10
  
  rbf_obs_1 = 1.1
  
  mass_1 = 0.00015
  
  pa_1 = +40
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 100
  
  nphis = 500
  
  nc = 50
  
  vch0 = -25
  
  vchf = -80
  
  nxs, nys = (200, 200)
  
  xpmax = 4
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  CIC = True
  
  tolfactor_vt = 3
  
  refpix = [125, 75]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.03
  
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
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
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
  
  mass_1 = 0.00015
  
  pa_1 = -20
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 100
  
  nphis = 500
  
  nc = 50
  
  vch0 = 35
  
  vchf = 65
  
  nxs, nys = (200, 200)
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  CIC = True
  
  tolfactor_vt = 3
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.05
  
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

Example 3
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
  
  distpc = 400
  
  vsys = + 5
  
  ra_source_deg, dec_source_deg = 84.095, -6.7675
  
  
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
  
  i_1 = 180-55
  
  L0_1 = 0.7
  
  zj_1 = 3
  
  vj_1 = 73
  
  va_1 = 0
  
  v0_1 = 4
  
  rbf_obs_1 = 1
  
  mass_1 = 0.00015
  
  pa_1 = -20
  
  """
  bowshock 2 [redshifted]
  """
  
  i_2 = 180-55
  
  L0_2 = 0.8
  
  zj_2 = 4
  
  vj_2 = 77
  
  va_2 = 0
  
  v0_2 = 4
  
  rbf_obs_2 = 1
  
  mass_2 = 0.00020 
  
  pa_2 = -20
  
  
  
  
  """
  SPECTRAL CUBE PARAMETERS
  """
  
  nzs = 100
  
  nphis = 500
  
  nc = 50
  
  vch0 = 30
  
  vchf = 57
  
  nxs, nys = (200, 200)
  
  xpmax = 5
  
  papv = pa_1
  
  bmaj, bmin = (0.420, 0.287)
  
  pabeam = -17.2
  
  vt = "2xchannel"
  
  CIC = True
  
  tolfactor_vt = 3
  
  refpix = [80, 30]
  
  coordcube = "sky"
  
  parot = 0
  
  sigma_beforeconv = 0.05
  
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
