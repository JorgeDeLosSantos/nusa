import numpy as np
import pytest

from nusa import Section


def test_section_with_area():
    section = Section(A=1e-3)
    assert section.A == 1e-3
    assert section.I is None


def test_section_with_second_moment():
    section = Section(I=2e-6)
    assert section.A is None
    assert section.I == 2e-6


def test_section_with_area_and_second_moment():
    section = Section(A=1e-3, I=2e-6, name='Main section')
    assert section.A == 1e-3
    assert section.I == 2e-6
    assert section.name == 'Main section'


def test_section_requires_at_least_one_geometric_property():
    with pytest.raises(ValueError):
        Section()


@pytest.mark.parametrize('A', [0, -1, np.nan, np.inf, -np.inf, 'invalid', True])
def test_section_rejects_invalid_area(A):
    with pytest.raises(ValueError):
        Section(A=A)


@pytest.mark.parametrize('I', [0, -1, np.nan, np.inf, -np.inf, 'invalid', True])
def test_section_rejects_invalid_second_moment(I):
    with pytest.raises(ValueError):
        Section(I=I)


def test_section_normalizes_numeric_properties_to_float():
    section = Section(A=1, I=2)
    assert isinstance(section.A, float)
    assert isinstance(section.I, float)


def test_section_rejects_non_string_name():
    with pytest.raises(TypeError):
        Section(A=1e-3, name=123)


def test_section_properties_are_read_only():
    section = Section(A=1e-3, I=2e-6, name='Section 01')
    with pytest.raises(AttributeError):
        section.A = 2e-3
    with pytest.raises(AttributeError):
        section.I = 3e-6
    with pytest.raises(AttributeError):
        section.name = 'Other'


def test_section_repr_omits_undefined_properties():
    section = Section(A=0.001)
    assert repr(section) == 'Section(A=0.001)'
