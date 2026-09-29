"""Finite-element problem model families."""

from .core import Model


class SpringModel(Model):
    """One-dimensional spring model."""

    displacement_dofs = ("ux",)
    force_dofs = ("fx",)

    def __init__(self, name="Spring Model 01"):
        super().__init__(name=name, mtype="spring")
        self.dof = 1

    def add_force(self, node, force):
        values = self._validated_component_vector(
            force, self.force_dofs, "force"
        )
        self._record_applied_forces(node, **values)

    def add_constraint(self, node, **constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        if values:
            self._record_prescribed_displacements(node, **values)


class BarModel(Model):
    """One-dimensional axial bar model."""

    displacement_dofs = ("ux",)
    force_dofs = ("fx",)

    def __init__(self, name="Bar Model 01"):
        super().__init__(name=name, mtype="bar")
        self.dof = 1

    def add_force(self, node, force):
        values = self._validated_component_vector(
            force, self.force_dofs, "force"
        )
        self._record_applied_forces(node, **values)

    def add_constraint(self, node, **constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        if values:
            self._record_prescribed_displacements(node, **values)


class TrussModel(Model):
    """Two-dimensional truss model."""

    displacement_dofs = ("ux", "uy")
    force_dofs = ("fx", "fy")

    def __init__(self, name="Truss Model 01"):
        super().__init__(name=name, mtype="truss")
        self.dof = 2

    def add_force(self, node, force):
        values = self._validated_component_vector(
            force, self.force_dofs, "force"
        )
        self._record_applied_forces(node, **values)

    def add_constraint(self, node, **constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        if values:
            self._record_prescribed_displacements(node, **values)


class BeamModel(Model):
    """Euler-Bernoulli beam model."""

    displacement_dofs = ("uy", "ur")
    force_dofs = ("fy", "m")

    def __init__(self, name="Beam Model 01"):
        super().__init__(name=name, mtype="beam")
        self.dof = 2

    def add_force(self, node, force):
        values = self._validated_component_vector(force, ("fy",), "force")
        self._record_applied_forces(node, **values)

    def add_moment(self, node, moment):
        values = self._validated_component_vector(moment, ("m",), "moment")
        self._record_applied_forces(node, **values)

    def add_constraint(self, node, **constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        if values:
            self._record_prescribed_displacements(node, **values)


class LinearTriangleModel(Model):
    """Two-dimensional constant-strain triangle model."""

    displacement_dofs = ("ux", "uy")
    force_dofs = ("fx", "fy")

    def __init__(self, name="LT Model 01"):
        super().__init__(name=name, mtype="triangle")
        self.dof = 2

    def add_force(self, node, force):
        values = self._validated_component_vector(
            force, self.force_dofs, "force"
        )
        self._record_applied_forces(node, **values)

    def add_constraint(self, node, **constraint):
        values = self._validated_named_components(
            constraint, self.displacement_dofs, "constraint"
        )
        if values:
            self._record_prescribed_displacements(node, **values)
