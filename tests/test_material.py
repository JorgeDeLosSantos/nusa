import numpy as np
import pytest

from nusa import Material


def test_material_with_required_property():
    material = Material(E=210e9)
    assert material.E == 210e9
    assert material.nu is None
    assert material.name is None


def test_material_with_optional_properties():
    material = Material(E=210e9, nu=0.3, name='Steel')
    assert material.E == 210e9
    assert material.nu == 0.3
    assert material.name == 'Steel'


def test_material_normalizes_numeric_properties_to_float():
    material = Material(E=210_000_000_000, nu=0)
    assert isinstance(material.E, float)
    assert isinstance(material.nu, float)


@pytest.mark.parametrize('E', [0, -1, np.nan, np.inf, -np.inf, 'invalid', True])
def test_material_rejects_invalid_E(E):
    with pytest.raises(ValueError):
        Material(E=E)


@pytest.mark.parametrize('nu', [-1.0, 0.5, np.nan, np.inf, -np.inf, 'invalid', True])
def test_material_rejects_invalid_nu(nu):
    with pytest.raises(ValueError):
        Material(E=210e9, nu=nu)


@pytest.mark.parametrize('nu', [-0.99, 0.0, 0.499])
def test_material_accepts_valid_poisson_ratio(nu):
    material = Material(E=210e9, nu=nu)
    assert material.nu == float(nu)


def test_material_rejects_non_string_name():
    with pytest.raises(TypeError):
        Material(E=210e9, name=123)


def test_material_properties_are_read_only():
    material = Material(E=210e9, nu=0.3, name='Steel')
    with pytest.raises(AttributeError):
        material.E = 70e9
    with pytest.raises(AttributeError):
        material.nu = 0.25
    with pytest.raises(AttributeError):
        material.name = 'Aluminum'


def test_material_repr():
    material = Material(E=210e9, nu=0.3, name='Steel')
    assert repr(material) == "Material(E=210000000000.0, nu=0.3, name='Steel')"
