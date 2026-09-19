"""NuSA: Numerical Structural Analysis in Python."""

from .version import __version__
from .core import Element, Model, Node
from .element import Bar, Beam, LinearTriangle, Spring, Truss
from .model import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    SpringModel,
    TrussModel,
)

__author__ = "P.J. De Los Santos"
__email__ = "delossantosmfq@gmail.com"

__all__ = [
    "__version__",
    "Model",
    "Element",
    "Node",
    "Spring",
    "Bar",
    "Truss",
    "Beam",
    "LinearTriangle",
    "SpringModel",
    "BarModel",
    "TrussModel",
    "BeamModel",
    "LinearTriangleModel",
]
