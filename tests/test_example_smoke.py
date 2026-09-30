"""Smoke tests for executable FEM examples."""

from pathlib import Path
import runpy

import matplotlib.pyplot as plt
import numpy as np
import pytest

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


@pytest.mark.parametrize(
    ("relative_path", "function_name"),
    [
        ("spring/simple_case.py", "simple_case"),
        ("spring/spring_02.py", "test2"),
        ("spring/spring_03.py", "test3"),
        ("bar/bar_2.py", "test2"),
        ("beam/beam_1.py", "test1"),
        ("beam/beam_3.py", "test3"),
        ("beam/beam_4.py", "test4"),
        ("beam/beam_5.py", "test5"),
    ],
)
def test_additional_0_4_examples_execute(relative_path, function_name, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    namespace = _load_example(relative_path)

    result = namespace[function_name]()

    assert isinstance(result, StaticResult)
    assert np.all(np.isfinite(result.displacements))
