"""Defines abstract base class for automatically handling model units."""

import abc
import functools
from collections.abc import Sequence
from typing import ClassVar

import numpy
from numpy.typing import NDArray

from ..util.type_ import Frozen
from . import system
from .system import UNITS, Dimension, Units, UnitsData


class UnitManager(Frozen, abc.ABC):

    _dimensions: ClassVar[dict[str, Dimension]]

    def __init__(self, units: UnitsData | None = None, **kwargs):
        if units is not None:
            units0 = Units.model_validate(units)
            for key, dim in self.__class__._dimensions.items():
                val0 = kwargs.get(key)
                kwargs[key] = system.convert_dimension_value(
                    units0, UNITS, dim, val0, **kwargs
                )

        super().__init__(**kwargs)


def manage_units(arg_dims: Sequence[Dimension], ret_dim: Dimension):
    """Transform function into a unit managing function.

    Converts arguments to internal units, calls the function, then converts the return
    value back to desired units.

    :param arg_dims: Argument dimensions
    :param ret_dim: Return dimensions
    :return: Function decorator
    """

    # def manage_units_(func0: UnitMethod) -> UnitMethod:
    def manage_units_(func0):

        @functools.wraps(func0)
        def func(
            self, *args, units: UnitsData | None = None, **kwargs
        ) -> float | NDArray[numpy.float64]:
            # If no units were specified, return as-is
            if units is None:
                return func0(self, *args, units=units, **kwargs)

            # Process units
            units0 = Units.model_validate(units)

            # Convert arguments
            args_ = tuple(
                system.convert_dimension_value(
                    units0, UNITS, dim, arg, **self.model_dump()
                )
                for dim, arg in zip(arg_dims, args, strict=True)
            )

            # Call function
            ret = func0(self, *args_, units=units, **kwargs)

            # Convert return
            ret_ = system.convert_dimension_value(
                UNITS, units0, ret_dim, ret, **self.model_dump()
            )
            return ret_

        return func

    return manage_units_
