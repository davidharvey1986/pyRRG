from .measure_moms import *
from .mmm import *
from .main import *
import pkg_resources  # part of setuptools
__version__ = pkg_resources.require("pyRRG")[0].version
