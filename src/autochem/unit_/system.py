"""Unit system."""

import functools
import itertools
from collections.abc import Mapping, Sequence
from typing import Annotated, TypeAlias

import more_itertools as mit
import pint
import pydantic
from pydantic_core import core_schema

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
    def volume(self) -> pint.Unit:
        """Volume unit."""
        return pint.Unit(self.length**3)

    @functools.cached_property
    def concentration(self) -> pint.Unit:
        """Concentration unit."""
        return pint.Unit(self.substance / self.volume)

    def rate_constant(self, order: int) -> pint.Unit:
        """Rate constant unit.

        :param order: Reaction order
        """
        if order == 1:
            return pint.Unit(self.time**-1)

        return pint.Unit(self.concentration ** (1 - order) * self.time**-1)


# Annotated type for use in pydantic models
Units_ = Annotated[
    pydantic.SkipValidation[Units],
    pydantic.BeforeValidator(lambda x: Units.model_validate(x)),
    pydantic.PlainSerializer(lambda x: Units.model_validate(x).model_dump()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(pydantic.BaseModel)),
    ),
]


# Alias for unit-convertible data
UnitData: TypeAlias = str | pint.Unit
UnitsData: TypeAlias = Mapping[str, UnitData] | Units


# Internal unit system (defaults of above Units class)
UNITS = Units()


# Constructors
def from_unit_sequence(unit_seq: Sequence[pint.Unit]) -> Units:
    """Construct unit sytem from sequence of composite units.

    :param unit_seq: Units sequence
    :param strict: Whether to raise an error if units are inconsistent
    :return: Unit system
    """
    # Split into components and chain them
    unit_iter = itertools.chain.from_iterable(
        map(pint.util.to_units_container, unit_seq),
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
