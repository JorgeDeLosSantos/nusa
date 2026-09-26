"""Regression tests for the explicit top-level NuSA public API."""

import importlib

import matplotlib as mpl
import nusa


EXPECTED_PUBLIC_API = {
    "__version__",
    "solve",
    "LinearStaticAnalysis",
    "StaticResult",
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
}


def test_top_level_public_api_is_explicit():
    assert set(nusa.__all__) == EXPECTED_PUBLIC_API

    namespace = {}
    exec("from nusa import *", namespace)
    exported = {name for name in namespace if not name.startswith("__")}
    assert exported == EXPECTED_PUBLIC_API - {"__version__"}
    assert namespace["__version__"] == nusa.__version__


def test_peripheral_modules_do_not_leak_into_top_level_namespace():
    for name in (
        "Modeler",
        "read_file",
        "read_msh",
        "ModelFromFiles",
        "read_model",
        "_read_truss_model",
        "_read_spring_model",
    ):
        assert not hasattr(nusa, name)


def test_mesh_remains_available_as_explicit_submodule():
    mesh = importlib.import_module("nusa.mesh")

    assert hasattr(mesh, "Modeler")


def test_removed_legacy_modules_are_not_importable():
    legacy_modules = (
        "nusa.io",
        "nusa._experimental",
        "nusa._3d_experimental",
        "nusa.graph",
        "nusa.templates",
        "nusa.lib",
        "nusa._lib",
    )

    for module_name in legacy_modules:
        try:
            importlib.import_module(module_name)
        except ModuleNotFoundError:
            pass
        else:
            raise AssertionError(f"{module_name} should not be importable")


def test_importing_nusa_does_not_override_matplotlib_rcparams():
    keys = ("figure.facecolor", "axes.facecolor", "font.size")
    original = {key: mpl.rcParams[key] for key in keys}
    sentinel = {
        "figure.facecolor": "pink",
        "axes.facecolor": "yellow",
        "font.size": 13.0,
    }

    try:
        mpl.rcParams.update(sentinel)
        importlib.reload(nusa)
        assert {key: mpl.rcParams[key] for key in keys} == sentinel
    finally:
        mpl.rcParams.update(original)
