"""Smoke tests for executable FEM examples."""

from pathlib import Path
import runpy

import numpy as np

from nusa import (
    BarModel,
    BeamModel,
    LinearTriangleModel,
    SpringModel,
    TrussModel,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def _load_example(relative_path):
    return runpy.run_path(str(EXAMPLES / relative_path))


def test_fem_examples_do_not_use_wildcard_imports():
    directories = ("spring", "bar", "truss", "beam", "linear_triangle")

    for directory in directories:
        for path in (EXAMPLES / directory).rglob("*.py"):
            text = path.read_text(encoding="utf-8")
            assert "import *" not in text, f"Wildcard import found in {path}"


def test_spring_example_executes_and_returns_solved_model():
    namespace = _load_example("spring/spring_01.py")
    model = namespace["test1"]()

    assert isinstance(model, SpringModel)
    assert model._is_assembled is True
    assert np.isfinite(model.nodes[2].ux)


def test_bar_example_executes_and_returns_solved_model():
    namespace = _load_example("bar/bar_1.py")
    model = namespace["test1"]()

    assert isinstance(model, BarModel)
    assert model._is_assembled is True
    assert np.isfinite(model.nodes[1].ux)


def test_truss_example_executes_and_returns_solved_model():
    namespace = _load_example("truss/truss_01.py")
    model = namespace["build_model"]()

    assert isinstance(model, TrussModel)
    assert model._is_assembled is True
    assert np.all(np.isfinite([model.nodes[0].ux, model.nodes[0].uy]))


def test_beam_example_executes_and_returns_solved_model():
    namespace = _load_example("beam/beam_2.py")
    model = namespace["test2"]()

    assert isinstance(model, BeamModel)
    assert model._is_assembled is True
    assert np.isfinite(model.nodes[1].uy)


def test_linear_triangle_example_executes_and_returns_solved_model():
    namespace = _load_example("linear_triangle/simple_triangle/simple_triangle.py")
    model = namespace["build_model"]()

    assert isinstance(model, LinearTriangleModel)
    assert model._is_assembled is True
    assert np.all(np.isfinite([model.nodes[1].ux, model.nodes[1].uy]))
