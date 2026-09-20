# ***********************************
#  Author: Pedro Jorge De Los Santos
#  E-mail: delossantosmfq@gmail.com
#  Blog: numython.github.io
#  License: MIT License
# ***********************************

import nusa._mesh as msh


class Modeler:
    """Small 2D geometry and triangular-mesh helper built around Gmsh."""

    def __init__(self):
        self.geom = msh.SimpleGMSH()

    def add_rectangle(self, p0, p1, esize=0.1):
        """Add an axis-aligned rectangular surface."""
        x0, y0 = p0[:2]
        x1, y1 = p1[:2]
        if x0 == x1 or y0 == y1:
            raise ValueError("Rectangle corners must define a nonzero area")

        xa, ya = x1, y0
        xb, yb = x0, y1
        p1_id = self.geom.add_point((x0, y0, 0), esize)
        p2_id = self.geom.add_point((xa, ya, 0), esize)
        p3_id = self.geom.add_point((x1, y1, 0), esize)
        p4_id = self.geom.add_point((xb, yb, 0), esize)
        l1 = self.geom.add_line(p1_id, p2_id)
        l2 = self.geom.add_line(p2_id, p3_id)
        l3 = self.geom.add_line(p3_id, p4_id)
        l4 = self.geom.add_line(p4_id, p1_id)
        loop = self.geom.add_line_loop(l1, l2, l3, l4)
        surface = self.geom.add_plane_surface(loop)
        return loop, surface

    def add_poly(self, *points, esize=0.1):
        """Add a polygonal surface from three or more 2D points."""
        if len(points) < 3:
            raise ValueError("A polygon requires at least three points")

        point_ids = [
            self.geom.add_point((point[0], point[1], 0), esize)
            for point in points
        ]
        lines = [
            self.geom.add_line(point_ids[index], point_ids[(index + 1) % len(point_ids)])
            for index in range(len(point_ids))
        ]
        loop = self.geom.add_line_loop(*lines)
        surface = self.geom.add_plane_surface(loop)
        return loop, surface

    def add_circle(self, p0, r, esize=0.1):
        """Add a circular surface."""
        if r <= 0:
            raise ValueError("Circle radius must be positive")

        xc, yc = p0[:2]
        center = self.geom.add_point((xc, yc, 0), esize)
        point = self.geom.add_point((xc + r, yc, 0), esize)
        circle = self.geom.add_circle(center, point)
        loop = self.geom.add_line_loop(circle)
        surface = self.geom.add_plane_surface(loop)
        return loop, surface

    def subtract_surfaces(self, outer, inner):
        """Create a plane surface with an inner hole."""
        outer_loop, outer_surface = outer
        inner_loop, inner_surface = inner
        self.geom.delete_surfaces(outer_surface, inner_surface)
        surface = self.geom.add_plane_surface(outer_loop, inner_loop)
        return (outer_loop, inner_loop), surface

    def _store_mesh(self, nodes, elements):
        self.nc = nodes
        self.ec = elements
        self.x = nodes[:, 0]
        self.y = nodes[:, 1]
        return nodes, elements

    def generate_mesh(self, verbose=False, gmsh_executable="gmsh"):
        """Generate a triangular mesh from the current geometry."""
        nodes, elements = self.geom.generate_mesh(
            verbose=verbose,
            gmsh_executable=gmsh_executable,
        )
        return self._store_mesh(nodes, elements)

    def generate_mesh_from_file(self, filename):
        """Load a triangular mesh file using meshio."""
        nodes, elements = msh.read_triangle_mesh(filename)
        return self._store_mesh(nodes, elements)

    def plot_mesh(self):
        """Plot the most recently generated or loaded triangular mesh."""
        if not hasattr(self, "nc") or not hasattr(self, "ec"):
            raise RuntimeError("Generate or load a mesh before calling plot_mesh()")

        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Polygon

        fig = plt.figure()
        ax = fig.add_subplot(111)
        patches = []
        for element in self.ec:
            points = [(self.nc[node, 0], self.nc[node, 1]) for node in element]
            patches.append(Polygon(points))

        collection = PatchCollection(
            patches,
            color="#25CDCD",
            edgecolor="#435959",
            alpha=0.8,
            lw=0.5,
        )
        ax.add_collection(collection)
        x0, x1, y0, y1 = self._rect_region()
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        plt.show()

    def _rect_region(self):
        x0, x1 = min(self.x), max(self.x)
        y0, y1 = min(self.y), max(self.y)
        dx = x1 - x0
        dy = y1 - y0
        kx = dx / 10.0 if dx else 0.1
        ky = dy / 10.0 if dy else 0.1
        return x0 - kx, x1 + kx, y0 - ky, y1 + ky
