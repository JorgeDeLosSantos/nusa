"""Smoke tests for executable FEM examples."""

from pathlib import Path
import runpy

import numpy as np

from nusa import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    SpringModel,
    StaticResult,
    TrussModel,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def _load_example(relative_path):
    return runpy.run_path(str(EXAMPLES / relative_path))


def test_fem_examples_do_not_use_wildcard_imports():
    directories = ("spring", "bar", "truss", "beam", "linear_triangle")

    for directory in directories:
        for file_path in (EXAMPLES / directory).rglob("*.py"):
            text = file_path.read_text(encoding="utf-8")
            assert "import *" not in text, f"Wildcard import found in {file_path}"


def test_spring_example_executes_and_returns_result():
    namespace = _load_example("spring/spring_01.py")
    result = namespace["test1"]()

    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))


def test_bar_example_executes_and_returns_result():
    namespace = _load_example("bar/bar_1.py")
    result = namespace["test1"]()

    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))


def test_truss_example_builds_model_and_solves_to_result():
    namespace = _load_example("truss/truss_01.py")
    model = namespace["build_model"]()

    assert isinstance(model, TrussModel)
    result = model.solve()
    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))


def test_beam_example_executes_and_returns_result():
    namespace = _load_example("beam/beam_2.py")
    result = namespace["test2"]()

    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))


def test_linear_triangle_example_builds_model_and_solves_to_result():
    namespace = _load_example(
        "linear_triangle/simple_triangle/simple_triangle.py"
    )
    model = namespace["build_model"]()

    assert isinstance(model, LinearTriangleModel)
    result = model.solve()
    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))
