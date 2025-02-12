"""Unit system."""

import enum
import functools
import itertools
from collections.abc import Mapping, Sequence
from typing import Any, TypeAlias

import more_itertools as mit
import numpy
import pint
from numpy.typing import NDArray

from ..util.type_ import Frozen, Unit_

UR = pint.UnitRegistry()


# Model for specifying units
class Units(Frozen):
    """Unit system."""

    # Core units
    time: Unit_ = pint.Unit("s")
    temperature: Unit_ = pint.Unit("K")
    length: Unit_ = pint.Unit("cm")
    substance: Unit_ = pint.Unit("mol")
    pressure: Unit_ = pint.Unit("atm")
    energy: Unit_ = pint.Unit("cal")

    # Derived units
    @functools.cached_property
    def energy_per_substance(self) -> pint.Unit:
        """Energy per substance unit."""
        return pint.Unit(self.energy / self.substance)

    @functools.cached_property
    def concentration(self) -> pint.Unit:
        """Concentration unit."""
        return pint.Unit(self.substance / self.length**3)

    def rate_constant(self, order: int) -> pint.Unit:
        """Rate constant unit.

        :param order: Reaction order
        """
        if order == 1:
            return pint.Unit(self.time**-1)

        return pint.Unit(self.concentration ** (1 - order) * self.time**-1)


# Alias for unit-convertible data
UnitData: TypeAlias = str | pint.Unit
UnitsData: TypeAlias = Mapping[str, UnitData] | Units


# Internal unit system (defaults of above Units class)
UNITS = Units()


# Named unit dimensions
class Dimension(enum.StrEnum):
    """Unit dimensions."""

    time = "time"
    temperature = "temperature"
    length = "length"
    substance = "substance"
    pressure = "pressure"
    energy = "energy"
    energy_per_substance = "energy_per_substance"
    concentration = "concentration"
    rate_constant = "rate_constant"


# Check that all dimensions are defined in `Units`
assert all(
    hasattr(UNITS, d) for d in map(str, Dimension)
), f"Invalid dimension name: {Dimension}"


# Constructors
def from_unit_sequence(unit_seq: Sequence[pint.Unit]) -> Units:
    """Construct unit sytem from sequence of composite units.

    :param unit_seq: Units sequence
    :param strict: Whether to raise an error if units are inconsistent
    :return: Unit system
    """
    # Split into components and chain them
    unit_iter = itertools.chain.from_iterable(
        map(pint.util.to_units_container, unit_seq)
    )
    unit_pool = list(map(pint.Unit, unit_iter))

    # For each dimension, gather matching components
    data = {}
    # for dim in UNITS.model_dump():
    for dim, unit0 in UNITS:
        # Separate out compatible units from the pool
        unit_pool, units = map(
            list, mit.partition(pint.Unit(unit0).is_compatible_with, unit_pool)
        )
        if units:
            # Pop the first one, make sure they match, and add it
            unit = pint.Unit(units.pop())
            assert all(map(unit.is_compatible_with, units)), f"{unit} !~ {units}"
            data[dim] = unit

    return Units.model_validate(data)


# Properties
def gas_constant_value(units: Units) -> float:
    """Determine gas constant value in unit system.

    :param units: Unit system
    :return: Gas constant value
    """
    return pint.Quantity("molar_gas_constant").m_as(
        units.energy_per_substance / units.temperature
    )


def dimension_unit(units: Units, dim: str | Dimension, **kwargs) -> pint.Unit:
    """Determine the unit for a dimension.

    :param units: Unit system
    :param dim: Dimension
    :param **kwargs: Extra arguments for unit determination
    :return: Unit
    """
    dim = Dimension(dim)
    if dim == Dimension.rate_constant:
        order = kwargs.get("order")
        assert order is not None, f"Missing 'order': {kwargs}"
        return units.rate_constant(order)
    return getattr(units, dim)


def convert_dimension_value(
    units: Units, new_units: Units, dim: str | Dimension, val: Any, **kwargs
) -> NDArray[numpy.float64]:
    """Convert dimension value to new units.

    :param units: Units sytem
    :param new_units: New units system
    :param dim: Dimension
    :param val: Value
    :param **kwargs: Extra arguments for unit determination
    :return: New value
    """
    unit = dimension_unit(units, dim, **kwargs)
    new_unit = dimension_unit(new_units, dim, **kwargs)
    new_val = pint.Quantity(val, unit).m_as(new_unit)
    return numpy.array(new_val, dtype=numpy.float64)


# # Helpers
# def dimension_compatible_unit(dim: str) -> pint.Unit | None:
#     """Get compatible units for a dimension.

#     :param dim: Dimension name
#     :return: Compatible units
#     """
#     return next(iter(UR.get_compatible_units(f"[{dim}]")), None)
