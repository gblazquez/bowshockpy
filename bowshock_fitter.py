import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, Normalize, ListedColormap
from matplotlib import cm

import astropy.units as u
from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord, concatenate
from astropy.wcs import WCS

import shelve

from datetime import datetime

import subprocess
result = subprocess.run(['whoami'], stdout=subprocess.PIPE)
user = result.stdout.decode().strip()

import sys
sys.path
sys.path.append('/home/guille/py_envs/radio310/lib/python3.10')

import bowshock as bs

today = datetime.now().strftime("%y%m%d")

rk = "prueba"
model = "nj"
ps = {
 'foldername': f"aj_prueba_{today}",
 'rj': 0.00,
 'L0': 0.50,
 'zj': None,
 'vj': None,
 'vw': 0,
 'v0': 17,
 'i': 20 * np.pi / 180,
 'vsys': 8.5,
 'vhead': -93,
 'xhead': 2.5,
 'xorigin': 0,
'distpc': 300,
'm0rate/rhow': None,
'm0/rhow': None,
'tj': None,
'tdyn': None,
'rbf': None,
'nzs':300
}

bsf = bs.BowshockFitter(
    ps,
    rk=rk,
    model=model,
)
# plt.show()
plt.show(block=False)
