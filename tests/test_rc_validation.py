"""Release-candidate validation tests for physical and assembly invariants."""

import numpy as np

from nusa import (
    Bar,
    BarModel,
    Beam,
    BeamModel,
    LinearTriangle,
    LinearTriangleModel,
    Node,
    Spring,
    SpringModel,
    Truss,
    TrussModel,
)


def _assert_symmetric_stiffness(model):
    model.assemble()
    np.testing.assert_allclose(
        model.stiffness_matrix,
        model.stiffness_matrix.T,
        rtol=0.0,
        atol=1e-12,
    )


def test_global_stiffness_is_symmetric_for_all_public_models():
    spring = SpringModel("spring symmetry")
    s1, s2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    spring.add_nodes([s1, s2])
    spring.add_element(Spring((s1, s2), 100.0))
    _assert_symmetric_stiffness(spring)

    bar = BarModel("bar symmetry")
    b1, b2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    bar.add_nodes([b1, b2])
    bar.add_element(Bar((b1, b2), E=200.0, A=3.0))
    _assert_symmetric_stiffness(bar)

    truss = TrussModel("truss symmetry")
    t1, t2 = Node((0.0, 0.0)), Node((3.0, 4.0))
    truss.add_nodes([t1, t2])
    truss.add_element(Truss((t1, t2), E=200.0, A=2.0))
    _assert_symmetric_stiffness(truss)

    beam = BeamModel("beam symmetry")
    bm1, bm2 = Node((0.0, 0.0)), Node((2.0, 0.0))
    beam.add_nodes([bm1, bm2])
    beam.add_element(Beam((bm1, bm2), E=200.0, I=4.0))
    _assert_symmetric_stiffness(beam)

    triangle = LinearTriangleModel("triangle symmetry")
    c1 = Node((0.0, 0.0))
    c2 = Node((1.0, 0.0))
    c3 = Node((0.0, 1.0))
    triangle.add_nodes([c1, c2, c3])
    triangle.add_element(
        LinearTriangle((c1, c2, c3), E=1000.0, nu=0.25, t=0.5)
    )
    _assert_symmetric_stiffness(triangle)


def test_spring_and_bar_reactions_balance_applied_loads():
    spring = SpringModel("spring equilibrium")
    s1, s2 = Node((0.0, 0.0)), Node((0.0, 0.0))
    spring.add_nodes([s1, s2])
    spring.add_element(Spring((s1, s2), 100.0))
    spring.add_constraint(s1, ux=0.0)
    spring.add_force(s2, (25.0,))
    spring.solve()

    assert np.isclose(
        spring.applied_loads.sum() + spring.reactions.sum(),
        0.0,
    )

    bar = BarModel("bar equilibrium")
    b1, b2, b3 = Node((0.0, 0.0)), Node((1.0, 0.0)), Node((2.0, 0.0))
    bar.add_nodes([b1, b2, b3])
    bar.add_elements([
        Bar((b1, b2), E=100.0, A=1.0),
        Bar((b2, b3), E=100.0, A=1.0),
    ])
    bar.add_constraint(b1, ux=0.0)
    bar.add_force(b3, (40.0,))
    bar.solve()

    assert np.isclose(
        bar.applied_loads.sum() + bar.reactions.sum(),
        0.0,
    )


def test_truss_and_triangle_reactions_balance_applied_force_components():
    truss = TrussModel("truss equilibrium")
    n1 = Node((0.0, 0.0))
    n2 = Node((1.0, 1.0))
    n3 = Node((2.0, 0.0))
    truss.add_nodes([n1, n2, n3])
    truss.add_elements([
        Truss((n1, n2), E=1000.0, A=1.0),
        Truss((n2, n3), E=1000.0, A=1.0),
        Truss((n1, n3), E=1000.0, A=1.0),
    ])
    truss.add_constraint(n1, ux=0.0, uy=0.0)
    truss.add_constraint(n3, ux=0.0, uy=0.0)
    truss.add_force(n2, (10.0, -30.0))
    truss.solve()

    applied = truss.applied_loads.reshape(-1, 2).sum(axis=0)
    reactions = truss.reactions.reshape(-1, 2).sum(axis=0)
    np.testing.assert_allclose(applied + reactions, [0.0, 0.0], atol=1e-10)

    triangle = LinearTriangleModel("triangle equilibrium")
    c1 = Node((0.0, 0.0))
    c2 = Node((1.0, 0.5))
    c3 = Node((0.0, 1.0))
    triangle.add_nodes([c1, c2, c3])
    triangle.add_element(
        LinearTriangle((c1, c2, c3), E=1000.0, nu=0.25, t=0.5)
    )
    triangle.add_constraint(c1, ux=0.0, uy=0.0)
    triangle.add_constraint(c3, ux=0.0, uy=0.0)
    triangle.add_force(c2, (12.0, -7.0))
    triangle.solve()

    applied = triangle.applied_loads.reshape(-1, 2).sum(axis=0)
    reactions = triangle.reactions.reshape(-1, 2).sum(axis=0)
    np.testing.assert_allclose(applied + reactions, [0.0, 0.0], atol=1e-10)


def test_beam_reactions_satisfy_force_and_moment_equilibrium():
    model = BeamModel("beam equilibrium")
    n1 = Node((0.0, 0.0))
    n2 = Node((3.0, 0.0))
    n3 = Node((5.0, 0.0))
    model.add_nodes([n1, n2, n3])
    model.add_elements([
        Beam((n1, n2), E=200.0, I=4.0),
        Beam((n2, n3), E=200.0, I=4.0),
    ])
    model.add_constraint(n1, uy=0.0, ur=0.0)
    model.add_force(n2, (-10.0,))
    model.add_moment(n3, (6.0,))
    model.solve()

    applied_force = sum(model.applied_load(node)["fy"] for node in model.nodes)
    reaction_force = sum(model.reaction(node)["fy"] for node in model.nodes)
    assert np.isclose(applied_force + reaction_force, 0.0, atol=1e-10)

    applied_moment = sum(
        model.applied_load(node)["m"] + node.x * model.applied_load(node)["fy"]
        for node in model.nodes
    )
    reaction_moment = sum(
        model.reaction(node)["m"] + node.x * model.reaction(node)["fy"]
        for node in model.nodes
    )
    assert np.isclose(applied_moment + reaction_moment, 0.0, atol=1e-10)
