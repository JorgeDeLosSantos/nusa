# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************

from pathlib import Path
import subprocess
import tempfile


def _load_meshio():
    try:
        import meshio
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "NuSA mesh utilities require meshio. "
            'Install them with: pip install "nusa[mesh]"'
        ) from exc
    return meshio


def read_triangle_mesh(filename):
    """Read mesh points and triangular connectivity with meshio."""
    meshio = _load_meshio()
    mesh = meshio.read(filename)

    try:
        triangles = mesh.cells_dict["triangle"]
    except KeyError as exc:
        raise ValueError(
            f"Mesh file {str(filename)!r} does not contain triangle cells"
        ) from exc

    return mesh.points, triangles


class SimpleGMSH:
    """Minimal Gmsh .geo builder used by :class:`nusa.mesh.Modeler`."""

    def __init__(self):
        self.ID_POINT = 0
        self.ID_LINE = 0
        self.ID_CIRCLE = 10000
        self.ID_LINE_LOOP = 0
        self.ID_PLANE_SURFACE = 0
        self.GMSH_CODE = []

    def add_point(self, coords, esize=0.1):
        self.ID_POINT += 1
        x, y = coords[0], coords[1]
        name = str(self.ID_POINT)
        self.GMSH_CODE.append(
            f"Point({name}) = {{ {x},{y},0,{esize} }};"
        )
        return name

    def add_line(self, p0, p1):
        self.ID_LINE += 1
        name = str(self.ID_LINE)
        self.GMSH_CODE.append(f"Line({name}) = {{ {p0},{p1} }};")
        return name

    def add_circle(self, p0, p1, p2=None):
        self.ID_CIRCLE += 1
        name = str(self.ID_CIRCLE)
        if p2 is None:
            self.GMSH_CODE.append(f"Circle({name}) = {{ {p1},{p0},{p1} }};")
        else:
            self.GMSH_CODE.append(f"Circle({name}) = {{ {p1},{p0},{p2} }};")
        return name

    def add_line_loop(self, *lines):
        self.ID_LINE_LOOP += 1
        name = str(self.ID_LINE_LOOP)
        self.GMSH_CODE.append(
            f"Line Loop({name}) = {{{','.join(lines)}}};"
        )
        return name

    def add_plane_surface(self, *loops):
        self.ID_PLANE_SURFACE += 1
        name = str(self.ID_PLANE_SURFACE)
        self.GMSH_CODE.append(
            f"Plane Surface({name}) = {{{','.join(loops)}}};"
        )
        return name

    def delete_surfaces(self, *surfaces):
        self.GMSH_CODE.append(
            f"Delete {{ Surface{{ {','.join(surfaces)} }}; }}"
        )

    def get_code(self):
        return "\n".join(self.GMSH_CODE)

    def generate_mesh(self, verbose=False, gmsh_executable="gmsh"):
        """Generate a 2D triangular mesh using the external Gmsh executable."""
        with tempfile.TemporaryDirectory(prefix="nusa-gmsh-") as tmpdir:
            tmpdir = Path(tmpdir)
            geo_path = tmpdir / "model.geo"
            msh_path = tmpdir / "model.msh"
            geo_path.write_text(self.get_code(), encoding="utf-8")

            command = [
                gmsh_executable,
                "-2",
                str(geo_path),
                "-format",
                "msh2",
                "-o",
                str(msh_path),
            ]

            try:
                result = subprocess.run(
                    command,
                    check=True,
                    text=True,
                    capture_output=not verbose,
                )
            except FileNotFoundError as exc:
                raise RuntimeError(
                    f"Gmsh executable {gmsh_executable!r} was not found on PATH"
                ) from exc
            except subprocess.CalledProcessError as exc:
                details = (exc.stderr or exc.stdout or "").strip()
                message = "Gmsh failed while generating the 2D mesh"
                if details:
                    message += f": {details}"
                raise RuntimeError(message) from exc

            if verbose and result.stdout:
                print(result.stdout, end="")

            if not msh_path.exists():
                raise RuntimeError("Gmsh completed without creating a mesh file")

            return read_triangle_mesh(msh_path)
