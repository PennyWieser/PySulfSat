__author__ = 'Penny Wieser'


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings as w
import scipy.optimize as optimize
from scipy.special import erf


# This has the core calculations, e.g., molar fractions, cation fractions etc.
from PySulfSat.core_calcs import *
# This has the loading excel function
from PySulfSat.import_data import *
# This has the actual SCSS calculations
from PySulfSat.scss_calcs2 import *
# This has stuff for performing s6+ corrections
from PySulfSat.s6_corrections import *
from PySulfSat.scas_calc import *
from PySulfSat.sulf_mass_balance import *
from PySulfSat.cali_plots import *
from PySulfSat.mantle_melting import *


# version
from ._version import __version__





