'''Material definitions for NuSA finite-element models.'''

import math


def _positive_finite(value, name):
    '''Return a finite positive scalar.'''
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a finite positive scalar')
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be a finite positive scalar') from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{name} must be a finite positive scalar')
    return value


def _poisson_ratio(value):
    '''Return a physically admissible isotropic Poisson ratio.'''
    message = 'nu must be a finite scalar in the range -1 < nu < 0.5'
    if isinstance(value, bool):
        raise ValueError(message)
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(message) from exc
    if not math.isfinite(value) or not (-1.0 < value < 0.5):
        raise ValueError(message)
    return value


class Material:
    '''Linear-elastic isotropic material properties.

    Parameters
    ----------
    E : float
        Young modulus. Must be finite and strictly positive.
    nu : float, optional
        Poisson ratio. When provided, must satisfy -1 < nu < 0.5.
    name : str, optional
        Descriptive metadata for the material.
    '''

    def __init__(self, E, nu=None, name=None):
        if name is not None and not isinstance(name, str):
            raise TypeError('Material name must be a string or None')

        self._E = _positive_finite(E, 'Young modulus E')
        self._nu = None if nu is None else _poisson_ratio(nu)
        self._name = name

    @property
    def E(self):
        '''Young modulus.'''
        return self._E

    @property
    def nu(self):
        '''Poisson ratio, or None when unspecified.'''
        return self._nu

    @property
    def name(self):
        '''Optional descriptive name.'''
        return self._name

    def __repr__(self):
        parts = [f'E={self.E!r}']
        if self.nu is not None:
            parts.append(f'nu={self.nu!r}')
        if self.name is not None:
            parts.append(f'name={self.name!r}')
        return f"Material({', '.join(parts)})"
