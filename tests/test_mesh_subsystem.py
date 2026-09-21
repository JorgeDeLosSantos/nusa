"""Regression tests for the optional mesh subsystem."""

from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from nusa import _mesh
from nusa.mesh import Modeler


def test_rectangle_geometry_generates_closed_surface_code():
    modeler = Modeler()

    loop, surface = modeler.add_rectangle((0.0, 0.0), (2.0, 1.0), esize=0.25)
    code = modeler.geom.get_code()

    assert loop == "1"
    assert surface == "1"
    assert code.count("Point(") == 4
    assert code.count("Line(") == 4
    assert "Line Loop(1)" in code
    assert "Plane Surface(1)" in code


def test_polygon_requires_at_least_three_points():
    modeler = Modeler()

    with pytest.raises(ValueError, match="at least three"):
        modeler.add_poly((0.0, 0.0), (1.0, 0.0))


def test_rectangle_requires_nonzero_area():
    modeler = Modeler()

    with pytest.raises(ValueError, match="nonzero area"):
        modeler.add_rectangle((0.0, 0.0), (0.0, 1.0))


def test_circle_requires_positive_radius():
    modeler = Modeler()

    with pytest.raises(ValueError, match="positive"):
        modeler.add_circle((0.0, 0.0), 0.0)


def test_circle_uses_four_quarter_arcs():
    modeler = Modeler()

    loop, surface = modeler.add_circle((1.0, 2.0), 0.5, esize=0.1)
    code = modeler.geom.get_code()

    assert loop == "1"
    assert surface == "1"
    assert code.count("Circle(") == 4
    assert "Line Loop(1) = {10001,10002,10003,10004};" in code


def test_subtract_surfaces_uses_outer_and_inner_loops():
    modeler = Modeler()
    outer = modeler.add_rectangle((0.0, 0.0), (2.0, 2.0))
    inner = modeler.add_circle((1.0, 1.0), 0.25)

    loops, surface = modeler.subtract_surfaces(outer, inner)
    code = modeler.geom.get_code()

    assert loops == (outer[0], inner[0])
    assert surface == "3"
    assert "Delete { Surface{" in code
    assert f"Plane Surface({surface}) = {{{outer[0]},{inner[0]}}};" in code
    assert not hasattr(modeler, "substract_surfaces")


def test_generate_mesh_stores_mesh_and_forwards_gmsh_options(monkeypatch):
    modeler = Modeler()
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    elements = np.array([[0, 1, 2]], dtype=int)
    calls = {}

    def fake_generate_mesh(*, verbose=False, gmsh_executable="gmsh"):
        calls["verbose"] = verbose
        calls["gmsh_executable"] = gmsh_executable
        return nodes, elements

    monkeypatch.setattr(modeler.geom, "generate_mesh", fake_generate_mesh)

    returned_nodes, returned_elements = modeler.generate_mesh(
        verbose=True,
        gmsh_executable="custom-gmsh",
    )

    assert calls == {"verbose": True, "gmsh_executable": "custom-gmsh"}
    assert returned_nodes is nodes
    assert returned_elements is elements
    assert modeler.nc is nodes
    assert modeler.ec is elements
    np.testing.assert_allclose(modeler.x, [0.0, 1.0, 0.0])
    np.testing.assert_allclose(modeler.y, [0.0, 0.0, 1.0])


def test_generate_mesh_from_file_uses_triangle_reader(monkeypatch, tmp_path):
    modeler = Modeler()
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    elements = np.array([[0, 1, 2]], dtype=int)
    mesh_path = tmp_path / "mesh.msh"

    monkeypatch.setattr(
        _mesh,
        "read_triangle_mesh",
        lambda filename: (nodes, elements),
    )

    returned_nodes, returned_elements = modeler.generate_mesh_from_file(mesh_path)

    assert returned_nodes is nodes
    assert returned_elements is elements
    assert modeler.nc is nodes
    assert modeler.ec is elements


def test_generate_mesh_from_file_reads_triangle_cells_with_meshio(tmp_path):
    import meshio

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    elements = np.array([[0, 1, 2]], dtype=int)
    mesh_path = tmp_path / "triangle.vtu"
    meshio.write_points_cells(mesh_path, nodes, [("triangle", elements)])

    modeler = Modeler()
    loaded_nodes, loaded_elements = modeler.generate_mesh_from_file(mesh_path)

    np.testing.assert_allclose(loaded_nodes, nodes)
    np.testing.assert_array_equal(loaded_elements, elements)


def test_generate_mesh_from_file_rejects_mesh_without_triangles(tmp_path):
    import meshio

    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ])
    lines = np.array([[0, 1]], dtype=int)
    mesh_path = tmp_path / "lines.vtu"
    meshio.write_points_cells(mesh_path, nodes, [("line", lines)])

    modeler = Modeler()
    with pytest.raises(ValueError, match="does not contain triangle"):
        modeler.generate_mesh_from_file(mesh_path)


def test_plot_mesh_requires_generated_or_loaded_mesh():
    modeler = Modeler()

    with pytest.raises(RuntimeError, match="Generate or load"):
        modeler.plot_mesh()


@pytest.mark.skipif(shutil.which("gmsh") is None, reason="Gmsh is not installed")
def test_real_gmsh_rectangle_mesh_smoke():
    modeler = Modeler()
    modeler.add_rectangle((0.0, 0.0), (1.0, 1.0), esize=0.5)

    nodes, elements = modeler.generate_mesh()

    assert nodes.ndim == 2
    assert nodes.shape[0] >= 3
    assert nodes.shape[1] >= 2
    assert elements.ndim == 2
    assert elements.shape[0] >= 1
    assert elements.shape[1] == 3
    assert elements.min() >= 0
    assert elements.max() < nodes.shape[0]


def test_simple_gmsh_reports_missing_executable(monkeypatch):
    geometry = _mesh.SimpleGMSH()
    geometry.add_point((0.0, 0.0, 0.0))

    monkeypatch.setattr(_mesh.shutil, "which", lambda executable: None)

    with pytest.raises(RuntimeError, match="gmsh --version"):
        geometry.generate_mesh(gmsh_executable="missing-gmsh")


def test_resolve_gmsh_command_uses_direct_executable(monkeypatch):
    monkeypatch.setattr(
        _mesh.shutil,
        "which",
        lambda executable: "/usr/bin/gmsh",
    )

    command = _mesh._resolve_gmsh_command(
        "gmsh",
        ["-2", "model.geo"],
        windows=False,
    )

    assert command == ["/usr/bin/gmsh", "-2", "model.geo"]


def test_resolve_gmsh_command_wraps_windows_batch_launcher(monkeypatch):
    monkeypatch.setattr(
        _mesh.shutil,
        "which",
        lambda executable: r"C:\Users\test\anaconda3\Scripts\gmsh.bat",
    )
    monkeypatch.setenv("COMSPEC", r"C:\Windows\System32\cmd.exe")

    command = _mesh._resolve_gmsh_command(
        "gmsh",
        ["-2", "model.geo"],
        windows=True,
    )

    assert command == [
        r"C:\Windows\System32\cmd.exe",
        "/d",
        "/c",
        r"C:\Users\test\anaconda3\Scripts\gmsh.bat",
        "-2",
        "model.geo",
    ]


def test_simple_gmsh_reports_gmsh_failure(monkeypatch):
    geometry = _mesh.SimpleGMSH()
    geometry.add_point((0.0, 0.0, 0.0))

    monkeypatch.setattr(_mesh.shutil, "which", lambda executable: "/usr/bin/gmsh")

    def failed_run(command, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=1,
            cmd=command,
            stderr="bad geometry",
        )

    monkeypatch.setattr(subprocess, "run", failed_run)

    with pytest.raises(RuntimeError, match="bad geometry"):
        geometry.generate_mesh()


def test_simple_gmsh_uses_temporary_msh2_output(monkeypatch):
    geometry = _mesh.SimpleGMSH()
    geometry.add_point((0.0, 0.0, 0.0))
    monkeypatch.setattr(_mesh.shutil, "which", lambda executable: "/usr/bin/gmsh")
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    elements = np.array([[0, 1, 2]], dtype=int)
    seen = {}

    def successful_run(command, **kwargs):
        seen["command"] = command
        output_path = Path(command[command.index("-o") + 1])
        output_path.write_text("placeholder", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="")

    def fake_reader(filename):
        assert Path(filename).exists()
        return nodes, elements

    monkeypatch.setattr(subprocess, "run", successful_run)
    monkeypatch.setattr(_mesh, "read_triangle_mesh", fake_reader)

    result_nodes, result_elements = geometry.generate_mesh()

    assert "-2" in seen["command"]
    assert "-format" in seen["command"]
    assert seen["command"][seen["command"].index("-format") + 1] == "msh2"
    assert result_nodes is nodes
    assert result_elements is elements
