"""Physical constants."""

import pint

from . import dim
from .dim import D, Dimension
from .system import Units


class C:
    """Pint constant names."""

    gas = "molar_gas_constant"
    boltzmann = "boltzmann_constant"


def quantity(const: str) -> pint.Quantity:
    """Get physical constant pint Quantity.

    :param const: Physical constant name
    :return: Pint quantity
    """
    const = const.strip().lower()
    return pint.Quantity(const)


def dimension(const: str) -> Dimension:
    """Get physical constant dimensions.

    :param const: Physical constant name
    :return: Dimension
    """
    # Define dimension for each constant
    dim_dct = {
        C.gas: D.energy_per_substance / D.temperature,
        C.boltzmann: D.energy / D.temperature,
    }

    # Look up dimension
    const = const.strip().lower()
    if const not in dim_dct:
        msg = f"Unknown constant: {const}"
        raise ValueError(msg)
    return dim_dct[const]


def value(const: str, units: Units) -> float:
    """Determine physical constant value in unit system.

    :param const: Physical constant name
    :param units: Unit system
    :return: Value
    """
    unit = dim.unit(units, dimension(const))
    return pint.Quantity(const).m_as(unit)
