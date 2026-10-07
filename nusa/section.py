'''Cross-section definitions for NuSA finite-element models.'''

import math


def _optional_positive_finite(value, name):
    '''Return None or a finite positive scalar.'''
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a finite positive scalar')
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be a finite positive scalar') from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{name} must be a finite positive scalar')
    return value


class Section:
    '''Cross-section properties used by structural elements.

    Parameters
    ----------
    A : float, optional
        Cross-sectional area.
    I : float, optional
        Second moment of area.
    name : str, optional
        Descriptive metadata for the section.
    '''

    def __init__(self, A=None, I=None, name=None):
        if name is not None and not isinstance(name, str):
            raise TypeError('Section name must be a string or None')

        self._A = _optional_positive_finite(A, 'Cross-sectional area A')
        self._I = _optional_positive_finite(I, 'Second moment of area I')
        if self._A is None and self._I is None:
            raise ValueError('Section requires at least one geometric property: A or I')
        self._name = name

    @property
    def A(self):
        '''Cross-sectional area, or None when unspecified.'''
        return self._A

    @property
    def I(self):
        '''Second moment of area, or None when unspecified.'''
        return self._I

    @property
    def name(self):
        '''Optional descriptive name.'''
        return self._name

    def __repr__(self):
        parts = []
        if self.A is not None:
            parts.append(f'A={self.A!r}')
        if self.I is not None:
            parts.append(f'I={self.I!r}')
        if self.name is not None:
            parts.append(f'name={self.name!r}')
        return f"Section({', '.join(parts)})"
