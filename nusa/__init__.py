# -*- coding: utf-8 -*-
"""
NuSA 0.2.0 (Numerical Structural Analysis in Python)
Author: Pedro Jorge De Los Santos
E-mail: delossantosmfq@gmail.com
Blog: jorgedelossantos.github.io // &  numython.github.io
License: MIT License
"""

__version__ = "0.2.0"
__author__ = "P.J. De Los Santos"
__email__ = "delossantosmfq@gmail.com"

from nusa.core import *
from nusa.element import *
from nusa.model import *
from _experimental import *
from nusa.mesh import *

import matplotlib as mpl
mpl.rc("figure", facecolor="#F9F9F9", titleweight="bold")
mpl.rc("axes", facecolor="#FFFFFF")
mpl.rc("font", size=10)

# Print NUSA info
#~ print("NuSA %s\n%s\n"%(__version__,10*"_"))
