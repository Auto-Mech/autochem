"""Unit dimensions."""

import functools
import operator
from collections.abc import Mapping
from typing import TypeAlias

import numpy as np
import pint
from numpy.typing import ArrayLike, NDArray

from .system import UNITS, Units


class Dimension:
    """Dimension class."""

    unit_exponents: dict[str, int]
    log: bool = False

    def __init__(
        self,
        arg: "str | Mapping[str, int] | Dimension | None" = None,
        *,
        log: bool | None = None,
        **kwargs: object,
    ) -> None:
        """Initialize dimension object."""
        if arg is not None:
            assert not kwargs, f"Invalid arguments: {arg}, {kwargs}"

        match arg:
            case Dimension():
                self.unit_exponents = arg.unit_exponents
                self.log = log or arg.log or False
            case str():
                self.unit_exponents = {arg: 1}
                self.log = log or False
            case Mapping():
                self.log = log or False
                self.unit_exponents = dict(arg)
            case _:
                assert arg is None, f"Cannot combine arg={arg} and kwargs={kwargs}"
                self.log = kwargs.pop("log", log or False)
                self.unit_exponents = kwargs

        assert all(hasattr(UNITS, k) for k in self.unit_exponents), (
            f"Invalid names: {self.unit_exponents}"
        )
        self.unit_exponents = {k: int(v) for k, v in self.unit_exponents.items()}

    def items(self) -> list[tuple[str, int]]:
        """Get items of dimension."""
        return sorted(self.unit_exponents.items(), key=lambda t: t[::-1], reverse=True)

    def __repr__(self) -> str:
        """Represent as string."""
        data_str = ", ".join(f"{k}={v}" for k, v in self.items())
        if self.log:
            data_str += ", log=True"
        return f"{self.__class__.__name__}({data_str})"

    def __pow__(self, power: int) -> "Dimension":
        """Raise dimension to a power.

        :param power: Power
        :return: New dimension
        """
        if self.log:
            msg = "Cannot raise logarithmic dimension to a power"
            raise ValueError(msg)
        unit_exponents = {k: v * power for k, v in self.unit_exponents.items()}
        return self.__class__(unit_exponents)

    def __mul__(self, other: "Dimension") -> "Dimension":
        """Multiply two dimensions together.

        :param other: Other dimension
        :return: New dimension
        """
        if self.log or other.log:
            msg = "Cannot multiply logarithmic dimensions"
            raise ValueError(msg)
        unit_exponents = {
            k: self.unit_exponents.get(k, 0) + other.unit_exponents.get(k, 0)
            for k in set(self.unit_exponents) | set(other.unit_exponents)
        }
        return self.__class__(unit_exponents)

    __rmul__ = __mul__

    def __truediv__(self, other: "Dimension") -> "Dimension":
        """Divide one dimension by another.

        :param other: Other dimension
        :return: New dimension
        """
        if self.log or other.log:
            msg = "Cannot divide logarithmic dimensions"
            raise ValueError(msg)
        data = {
            k: self.unit_exponents.get(k, 0) - other.unit_exponents.get(k, 0)
            for k in set(self.unit_exponents) | set(other.unit_exponents)
        }
        return self.__class__(data)


# Alias for dimension-convertible data
DimensionData: TypeAlias = str | Mapping[str, int] | Dimension


class D:
    """Easy access to common dimensions."""

    time = Dimension("time")
    temperature = Dimension("temperature")
    length = Dimension("length")
    substance = Dimension("substance")
    pressure = Dimension("pressure")
    energy = Dimension("energy")
    energy_per_substance = Dimension("energy_per_substance")
    volume = Dimension("volume")
    concentration = Dimension("concentration")
    rate_constant = Dimension("rate_constant")


def unit(units: Units, dim: DimensionData, **kwargs: object) -> pint.Unit:
    """Determine the unit for a dimension.

    :param units: Unit system
    :param dim: Dimension
    :param **kwargs: Extra arguments for unit determination
    :return: Unit
    """
    dim = Dimension(dim)

    def _unit(dim_name: str) -> pint.Unit:
        if dim_name == "rate_constant":
            order = kwargs.get("order")
            assert order is not None, f"Missing 'order': {kwargs}"
            return units.rate_constant(order)
        return getattr(units, dim_name)

    return functools.reduce(operator.mul, (_unit(k) ** v for k, v in dim.items()))


unit_ = unit


def conversion_factor(
    units: Units, new_units: Units, dim: DimensionData, **kwargs: object
) -> float:
    """Convert dimension value to new units.

    For log dimensions, this is the conversion factor for the argument of the log
    function. The conversion is performed by taking the log of this factor and adding it
    to the original value.

    :param units: Units sytem
    :param new_units: New units system
    :param dim: Dimension
    :param **kwargs: Extra arguments for unit determination
    :return: New value
    """
    unit = unit_(units, dim, **kwargs)
    new_unit = unit_(new_units, dim, **kwargs)
    return pint.Quantity(1, unit).m_as(new_unit)


def convert(
    units: Units,
    new_units: Units,
    dim: DimensionData,
    val: ArrayLike,
    **kwargs: object,
) -> NDArray[np.float64]:
    """Convert dimension value to new units.

    :param units: Units sytem
    :param new_units: New units system
    :param dim: Dimension
    :param val: Value
    :param **kwargs: Extra arguments for unit determination
    :return: New value
    """
    dim = Dimension(dim)
    factor = conversion_factor(units=units, new_units=new_units, dim=dim, **kwargs)
    val = np.array(val, dtype=np.float64)
    return val * factor if not dim.log else val + np.log(factor)


def log(dim: DimensionData) -> Dimension:
    """Create logarithmic dimension.

    :param dim: Dimension
    :return: Logarithmic dimension
    """
    return Dimension(dim, log=True)
