"""Initialization of BowshockPy"""

from .cubemass import MassCube
from .cubeproc import CubeProcessing
from .genbow import generate_bowshock
from .inputfiles import *
from .modelproj import ObsModel
from .models import BaseModel, BowshockModel
from .moments import mom0, mom1, mom2, maxintens, pv, sumint
from .utils import gaussconvolve, get_color, mb_sa_gaussian_f, print_example
from .version import __version__
