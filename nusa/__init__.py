"""NuSA: Numerical Structural Analysis in Python."""

from .version import __version__
from .analysis import LinearStaticAnalysis, solve
from .result import StaticResult
from .reporting import simple_report
from .post import element_field, nodal_field
from .visualization import (
    plot_deformed_shape,
    plot_element_field,
    plot_moment_diagram,
    plot_nodal_field,
    plot_shear_diagram,
)
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
    "solve",
    "LinearStaticAnalysis",
    "StaticResult",
    "simple_report",
    "element_field",
    "nodal_field",
    "plot_deformed_shape",
    "plot_nodal_field",
    "plot_element_field",
    "plot_moment_diagram",
    "plot_shear_diagram",
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
