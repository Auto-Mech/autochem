"""Rate constant models."""

import abc
import warnings
from collections.abc import Mapping, Sequence
from numbers import Number
from typing import Annotated, ClassVar, Literal, Self

import altair as alt
import more_itertools as mit
import numpy as np
import pint
import pydantic
from numpy.polynomial import chebyshev
from numpy.typing import ArrayLike, NDArray
from pydantic import BeforeValidator, model_validator
from pydantic_core import core_schema

from .. import unit_
from ..unit_ import UNITS, C, D, Dimension, UnitManager, Units, UnitsData, const
from ..util import arrh, chemkin, func, mess, plot
from ..util.type_ import NDArray_, Scalable, Scalers, SubclassTyped
from . import blend
from .blend import BlendingFunction_


class Key:
    """Attribute keys."""

    # Independent variables
    T = "T"
    P = "P"
    # Dependent variables
    k = "k"


class BaseRate(UnitManager, Scalable, SubclassTyped, abc.ABC):
    """Abstract base class for rate constants."""

    order: int = 1

    @property
    def unit(self) -> pint.Unit:
        """Rate unit."""
        return UNITS.rate_constant(self.order)

    def __init__(self, units: UnitsData | None = None, **kwargs: object) -> None:
        """Instantiate rate."""
        super().__init__(units=units, **kwargs)

    @abc.abstractmethod
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float128]:
        """Evaluate rate constant.

        :param T: Temperature(s)
        :param P: Pressure(s)
        :param units: Input units and desired output units
        :return: Value(s)
        """

    @property
    def plot_mark(self) -> str:
        """Plot mark to use in altair."""
        return plot.Mark.line

    def plot_data(
        self,
        T: float | tuple[float, float] = (400, 1250),  # noqa: N803
        P: float | tuple[float, float] = 1,  # noqa: N803
        units: UnitsData | None = None,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Display as an Arrhenius plot.

        :param T: Temperature or temperature range
        :param P: Pressure or pressure range
        :param units: Units
        :return: Chart
        """
        if isinstance(T, Sequence) and isinstance(P, Number):
            P_ = P
            T_ = np.linspace(*T, 1000)  # noqa: N806
            x_data = T_
        elif isinstance(T, Number) and isinstance(P, Sequence):
            P_ = np.linspace(*P, 1000)  # noqa: N806
            T_ = T
            x_data = P_
        else:
            msg = f"3-dimensional plotting not yet implemented:\nT={T}\nP={P}"
            raise ValueError(msg)

        y_data = self(T_, P_, units=units)
        return x_data, y_data  # type: ignore

    def display(  # noqa: PLR0913
        self,
        T_range: tuple[float, float] = (400, 1250),  # noqa: N803
        P: float = 1,  # noqa: N803
        units: UnitsData | None = None,
        label: str | None = None,
        color: str | None = None,
        x_label: str = "1000/𝑇",  # noqa: RUF001
        y_label: str = "𝑘",  # noqa: RUF001
    ) -> alt.Chart:
        """Display as an Arrhenius plot.

        :param T_range: Temperature range
        :param P: Pressure
        :param units: Units
        :param x_label: X-axis label
        :param y_label: Y-axis label
        :return: Chart
        """
        T, k = self.plot_data(T=T_range, P=P, units=units)  # noqa: N806
        return plot.arrhenius(
            ks=[k],
            T=T,
            order=self.order,
            units=units,
            labels=[label] if label else None,
            colors=[color] if color else None,
            x_label=x_label,
            y_label=y_label,
            mark=self.plot_mark,
        )


def nan_array_to_none(arr: ArrayLike | None) -> NDArray | None:
    """Replace an array of NaNs with None.

    :param arr: Array or None
    :return: Array or None
    """
    return None if arr is None or np.all(np.isnan(arr)) else arr


def drop_invalid_rates(arr: ArrayLike | None) -> NDArray | None:
    """Replace negative rates with None.

    Drop invalid rates up to one past the last negative temperature.

    :param arr: Array or None
    :return: Array or None
    """
    if arr is None:
        return None
    arr = np.array(arr, dtype=float, copy=True)
    # Define mask to select up to one past the last negative value per row
    axis = 0
    shape = arr.shape
    row_idxs = np.expand_dims(np.arange(shape[axis]), axis=1)
    neg_idx = np.max(np.where(arr < 0, row_idxs, -2), axis=axis)
    mask = row_idxs <= (np.expand_dims(neg_idx, axis=axis) + 1)
    # Set those values to nan
    arr[mask] = np.nan
    return arr


def multiply_if_not_none(obj1: ArrayLike, obj2: ArrayLike) -> NDArray[np.float64]:
    if obj1 is None:
        return obj1
    if obj2 is None:
        return obj2
    return np.multiply(obj1, obj2).tolist()


class Rate(BaseRate):
    """Rate data."""

    T: list[float]
    P: list[float]
    k_data: Annotated[NDArray_, BeforeValidator(drop_invalid_rates)]
    k_high: Annotated[list[float] | None, BeforeValidator(nan_array_to_none)] = None

    # Private attributes
    type_: ClassVar[str] = "data"
    _scalers: ClassVar[Scalers] = {
        "k_data": np.multiply,
        "k_high": multiply_if_not_none,
    }
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "T": D.temperature,
        "P": D.pressure,
        "k_data": D.rate_constant,
        "k_high": D.rate_constant,
    }

    @model_validator(mode="after")
    def sort_temperatures(self) -> Self:
        idxs = np.argsort(self.T)
        self.T = np.take(self.T, idxs).tolist()
        self.k_data = self.k_data[idxs]
        self.k_high = (
            None if self.k_high is None else np.take(self.k_high, idxs).tolist()
        )
        return self

    def __truediv__(self, other: "Rate" | ArrayLike) -> Self:
        """Scalar division.

        :param c: Scalar value to divide by
        :return: Scaled object
        """
        if not isinstance(other, Rate):
            return super().__truediv__(other)

        k_data = self.k_data / other.k_data
        k_high = None
        if self.k_high is not None and other.k_high is not None:
            k_high = np.true_divide(self.k_high, other.k_high).tolist()

        return self.model_copy(deep=True, update={"k_data": k_data, "k_high": k_high})

    @property
    def plot_mark(self) -> str:
        """Plot mark to use in altair."""
        return plot.Mark.point

    def plot_data(
        self,
        T: float | tuple[float, float] = (400, 1250),  # noqa: N803
        P: float | tuple[float, float] = 1,  # noqa: N803
        units: UnitsData | None = None,
    ) -> tuple[NDArray, NDArray]:
        """Display as an Arrhenius plot.

        :param T: Temperature or temperature range
        :param P: Pressure or pressure range
        :param units: Units
        :return: Chart
        """
        if isinstance(T, Sequence) and isinstance(P, Number):
            T_ = self.T
            P_ = P
            (i_,) = np.where(np.greater_equal(T_, T[0]) & np.less_equal(T_, T[1]))
            x_data = np.take(T_, i_)
        elif isinstance(T, Number) and isinstance(P, Sequence):
            T_ = T
            P_ = self.P
            (i_,) = np.where(np.greater_equal(P_, P[0]) & np.less_equal(P_, P[1]))
            x_data = np.take(P_, i_)
        else:
            msg = f"3-dimensional plotting not yet implemented:\nT={T}\nP={P}"
            raise ValueError(msg)

        y_data = np.take(self(T_, P_, units=units), i_)
        return x_data, y_data

    @unit_.manage_units([D.temperature, D.pressure], D.rate_constant)
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray:
        """Evaluate rate constant."""
        interp_ = arrh.interpolator(self.T, self.P, self.k_data)
        return interp_(T, P)

    def __add__(self, other: "Rate") -> "Rate":
        """Add rates."""
        assert self.order == other.order, f"{self} !~ {other}"
        T, ixT1, ixT2 = np.intersect1d(self.T, other.T, return_indices=True)  # noqa: N806
        P, ixP1, ixP2 = np.intersect1d(self.P, other.P, return_indices=True)  # noqa: N806
        k_data1 = self.k_data[np.ix_(ixT1, ixP1)]
        k_data2 = other.k_data[np.ix_(ixT2, ixP2)]
        k_data = np.add(k_data1, k_data2)

        k_high = None
        if self.k_high is not None and other.k_high is not None:
            k_high1 = np.array(self.k_high)[ixT1]
            k_high2 = np.array(other.k_high)[ixT2]
            k_high = np.add(k_high1, k_high2)

        return self.__class__(order=self.order, T=T, P=P, k_data=k_data, k_high=k_high)

    def without_nan(self) -> "Rate":
        """Return a copy of the rate without temperatures giving NaNs.

        :return: Rate
        """
        k_data = self.k_data
        k_all = (
            k_data if self.k_high is None else np.column_stack((k_data, self.k_high))
        )
        not_nan = np.all(np.isfinite(k_all), axis=0)
        return self.__class__(
            order=self.order,
            T=np.array(self.T)[not_nan],
            P=self.P,
            k_data=self.k_data[:, not_nan],
        )

    def fill_nan(self, nan: float = 0.0) -> "Rate":
        """Return a copy of the rate with NaNs replaced with a value.

        :return: Rate
        """
        k_data = np.nan_to_num(self.k_data, nan=nan)
        k_high = (
            None
            if self.k_high is None
            else np.nan_to_num(self.k_high, nan=nan).tolist()
        )
        return self.__class__(
            order=self.order, T=self.T, P=self.P, k_data=k_data, k_high=k_high
        )

    def high_pressure_values(self) -> NDArray[np.float128]:
        """Return high-pressure rate values."""
        k_max = np.array(self(T=self.T, P=max(self.P)))
        return np.array(self.k_high) if self.k_high else k_max

    def is_pressure_dependent(
        self, T: Sequence[float] | None = None, tol: float = 0.2
    ) -> bool:
        """Determine whether or not the rate is pressure dependent.

        If data fittable data is only available for one pressure, the reaction
        will be treated as pressure-dependent.

        :param tol: Threshold for determining pressure dependence
        :return: `True` if it is, otherwise `False`
        """
        if self.is_empty():
            return False

        T_ = self.T if T is None else T
        P_ = self.fittable_pressures()
        data = self(T=T_, P=P_)

        # Mask of valid (non-NaN) entries
        mask = ~np.isnan(data)
        dim0 = np.arange(data.shape[0])
        dim1 = np.arange(data.shape[1])

        # Indices of first and last valid values per row
        idx_lo = np.where(mask, dim1, np.inf).argmin(axis=1)
        idx_hi = np.where(mask, dim1, -np.inf).argmax(axis=1)

        if np.all(idx_lo == idx_hi):
            return True

        # Extract first and last values, handling rows with all-NaN
        k_lo = np.where(mask.any(axis=1), data[dim0, idx_lo], np.nan)
        k_hi = np.where(mask.any(axis=1), data[dim0, idx_hi], np.nan)

        diff = np.abs(k_lo - k_hi) / k_lo
        return bool(np.any(diff > tol))

    def fittable_pressures(self) -> list[float]:
        """Identify pressures with enough data points to fit.

        :return: Pressures
        """
        count = np.sum(np.isfinite(self.k_data), axis=0)
        return np.array(self.P)[count >= 3].tolist()

    def unfittable_pressures(self) -> list[float]:
        """Identify pressures without enough data points to fit.

        :return: Pressures
        """
        count = np.sum(np.isfinite(self.k_data), axis=0)
        return np.array(self.P)[count < 3].tolist()

    def is_empty(self) -> bool:
        """Check whether rate is empty (all NaN or no values).

        :return: Boolean
        """
        if np.size(self.k_data) == 0:
            return True

        return np.all(np.isnan(self.k_data)).item()

    def has_pressure_data(self, P: Sequence[float]) -> bool:
        """Check for presence of pressures.

        :param P: Pressures
        :return: Boolean
        """
        if np.size(P) == 0:
            return True

        P0_ = np.expand_dims(np.array(self.P), axis=0)
        P_ = np.expand_dims(np.array(P), axis=1)
        close_arr = np.isclose(P0_, P_)
        if not np.all(np.any(close_arr, axis=1)).item():
            return False
        select = np.any(close_arr, axis=0)
        isfinite = np.isfinite(self.k_data[:, select])
        return np.all(np.any(isfinite, axis=0)).item()

    def clear_pressure_range(self, P0: float, P1: float) -> "Rate":
        """Clear pressure range."""
        P_mid = [p for p in self.P if P0 < p and p < P1]
        return self.clear_pressures([P0, *P_mid, P1])

    def clear_pressures(self, P: Sequence[float]) -> "Rate":
        """Clear pressures.

        :param P: Pressures
        :return: Rate
        """
        P_orig = np.array(self.P)
        P_orig_ = np.expand_dims(P_orig, axis=1)
        P_drop_ = np.expand_dims(P, axis=0)
        drop = np.any(np.isclose(P_orig_, P_drop_), axis=1)
        k_data = self.k_data.copy()
        k_data[:, drop] = np.nan
        return self.__class__(
            order=self.order, T=self.T, P=self.P, k_data=k_data, k_high=self.k_high
        )

    def drop_pressures(self, P: Sequence[float]) -> "Rate":
        """Drop pressures.

        :param P: Pressures
        :return: Rate
        """
        P_orig = np.array(self.P)
        P_orig_ = np.expand_dims(P_orig, axis=1)
        P_drop_ = np.expand_dims(P, axis=0)
        keep = ~np.any(np.isclose(P_orig_, P_drop_), axis=1)
        P_keep = P_orig[keep].tolist()
        k_data = self.k_data[:, keep].copy()
        return self.__class__(
            order=self.order, T=self.T, P=P_keep, k_data=k_data, k_high=self.k_high
        )

    def drop_temperatures(self, T: Sequence[float]) -> "Rate":
        """Drop temperatures.

        :param T: Temperatures to drop
        :return: Rate object
        """
        T_orig = np.array(self.T)
        T_orig_ = np.expand_dims(T_orig, axis=1)
        T_drop_ = np.expand_dims(T, axis=0)
        keep = ~np.any(np.isclose(T_orig_, T_drop_), axis=1)
        T_keep = T_orig[keep].tolist()
        k_data = self.k_data[keep, :].copy()
        k_high = None if self.k_high is None else np.extract(keep, self.k_high).tolist()
        return self.__class__(
            order=self.order, T=T_keep, P=self.P, k_data=k_data, k_high=k_high
        )

    def clear(self) -> "Rate":
        """Return a cleared copy of the rate (all values set to NaN)."""
        k_data = np.full_like(self.k_data, np.nan, dtype=float)
        return self.model_copy(update={"k_data": k_data, "k_high": None})

    def merge_equivalent(self, other: "Rate", *, tol: float = 0.1) -> "Rate":
        """Merge equivalent rates.

        Matching rates are averaged. Mismatched rates are replaced with NaN.
        Matches are determined based on a tolerance threshold.

        :param other: Rate
        :return: Rate
        """
        k_data_avg = (self.k_data + other.k_data) / 2
        k_data_diff = np.abs(self.k_data - other.k_data)
        match = (k_data_diff / k_data_avg) < tol
        k_data = np.where(match, k_data_avg, np.nan)
        k_high = None
        if self.k_high and other.k_high:
            k_high_avg = np.add(self.k_high, other.k_high) / 2
            k_high_diff = np.abs(np.subtract(self.k_high, other.k_high))
            match = (k_high_diff / k_high_avg) < tol
            k_high = np.where(match, k_high_avg, np.nan).tolist()
        return self.model_copy(update={"k_data": k_data, "k_high": k_high})


class BoundedMixin(pydantic.BaseModel):
    """Mixin to define bounded calculator."""

    T_min: float | None = None
    T_max: float | None = None

    def in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> NDArray[np.bool_]:
        """Determine whether temperature(s) are in bounds.

        :param T: Temperature(s)
        :return: Boolean value(s)
        """
        T = np.array(T, dtype=np.float64)  # noqa: N806
        greater_than_min = (
            np.ones_like(T, dtype=bool) if self.T_min is None else self.T_min <= T
        )
        less_than_max = (
            np.ones_like(T, dtype=bool) if self.T_max is None else self.T_max >= T
        )
        return greater_than_min & less_than_max

    def all_in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> bool:
        """Determine whether all temperature(s) are in bounds.

        :param T: Temperature(s)
        :return: `True` if they are
        """
        return np.all(self.in_bounds(T)).item()

    def assert_all_in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> None:
        """Assert that all temperature(s) are in bounds.

        :param T: Temperature(s)
        """
        assert self.all_in_bounds(T), f"{self.T_min} !<= {T} !<= {self.T_max}"


class RateFit(BaseRate, BoundedMixin):
    """Rate fit abstract base classs."""

    efficiencies: dict[str, float] = pydantic.Field(default_factory=dict)

    @property
    def third_body(self) -> str | None:
        """Get third body."""
        eff = self.efficiencies

        if eff.get("M") == 1:
            return "M"

        return next((c for c, e in eff.items() if e == 1.0), None)

    @property
    def is_pressure_dependent(self) -> bool:
        """Determine if the rate is pressure dependent."""
        return True

    @pydantic.field_validator("efficiencies", mode="before")
    @classmethod
    def _sanitize_efficiencies(cls, value: object) -> object:
        if isinstance(value, Mapping):
            return {k: v for k, v in value.items() if v is not None}
        return value

    def is_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was cleared."""
        msg = f"Cleared checking not implemented for {self.__class__.__name__}"
        raise NotImplementedError(msg)

    def is_partially_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was partially cleared."""
        msg = f"Cleared checking not implemented for {self.__class__.__name__}"
        raise NotImplementedError(msg)


class ArrheniusRateFit(RateFit):
    """Arrhenius rate fit."""

    A: float = 1.0
    b: float = 0.0
    E: float = 0.0

    # Private attributes
    type_: ClassVar[str] = "arrhenius"
    _scalers: ClassVar[Scalers] = {"A": np.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "A": D.rate_constant,
        "E": D.energy_per_substance,
    }

    @property
    def is_pressure_dependent(self) -> bool:
        """Determine if the rate is pressure dependent."""
        return False

    @unit_.manage_units([D.temperature, D.pressure], D.rate_constant)
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float128]:
        """Evaluate rate constant."""
        T_, _ = func.normalize_arguments((T, P))  # noqa: N806
        T_ = np.where(self.in_bounds(T_), T_, np.nan)
        R = const.value(C.gas, UNITS)  # noqa: N806
        kTP = self.A * (T_**self.b) * np.exp(-self.E / (R * T_))  # noqa: N806
        return func.normalize_values(kTP, (T, P))

    @classmethod
    @unit_.manage_units([D.temperature, D.rate_constant])
    def fit(
        cls,
        Ts: ArrayLike,  # noqa: N803
        ks: ArrayLike,
        *,
        A_fill: float | None = None,
        bad_fit: Literal["fill"]
        | Literal["warn"]
        | Literal["raise"]
        | Literal["ignore"] = "warn",
        validate: bool = True,
        order: int = 1,
        units: UnitsData | None = None,  # noqa: ARG003
    ) -> "ArrheniusRateFit":
        """Fit data to Arrhenius rate fit.

        :param T: Temperatures
        :param k: Rates
        :param A_fill: Optional dummy parameter for unfittable rates
            (Otherwise, an error will be thrown.)
        :param bad_fit: How to handle bad fits that create floating point errors;
            (Options: "fill", replace with the fill value, or numpy.seterr options)
        :return: Rate fit
        """
        T = np.array(Ts, dtype=np.float64)  # noqa: N806
        _1 = np.ones_like(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        M = np.column_stack([_1, np.log(T), -1 / (R * T)])  # noqa: N806
        v = np.log(ks)

        ok = np.isfinite(v)
        M = M[ok, :]  # noqa: N806
        v = v[ok]

        if len(v) < 3:
            if A_fill is None:
                msg = (
                    f"Cannot fit with fewer than 3 data points: {v}\n"
                    "You can circumvent this by setting the A_fill parameter "
                    "as a placeholder for unfittable rates."
                )
                raise ValueError(msg)
            return cls(order=order, A=A_fill, b=0, E=0)

        (lnA, b, E), *_ = np.linalg.lstsq(M, v, rcond=1e-24)  # noqa: N806
        obj = cls(order=order, A=np.exp(lnA), b=b, E=E)

        if not validate or bad_fit == "fill":
            try:
                with np.errstate(all="raise"):
                    vals = obj(T=T)
                    if not np.all(np.isfinite(vals)):
                        msg = "Fitted rate gives non-finite values over the input temperature range."
                        raise FloatingPointError(msg)
            except FloatingPointError as e:
                if A_fill is not None:
                    return cls(order=order, A=A_fill, b=0, E=0)

                if A_fill is None:
                    msg = f"{e}\nOption `bad_fit='fill'` was given, but no A_fill value was provided."
                else:
                    msg = f"{e}\nNote: You can handle errors using the `bad_fit` keyword argument."

                raise FloatingPointError(msg)
        else:
            with np.errstate(all=bad_fit):
                obj(T=T)

        return obj

    def is_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was cleared."""
        return np.allclose([self.A, self.b, self.E], [A_fill, 0.0, 0.0], atol=0)

    def is_partially_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was partially cleared."""
        return self.is_cleared(A_fill=A_fill)


class FalloffRateFit(RateFit, abc.ABC):  # type: ignore[misc]
    """Falloff rate fit."""

    A_high: float
    b_high: float
    E_high: float
    A_low: float
    b_low: float
    E_low: float
    function: BlendingFunction_
    activated: bool = False

    # Private attributes
    type_: ClassVar[str] = "falloff"
    _scalers: ClassVar[Scalers] = {"A_high": np.multiply, "A_low": np.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "A_high": D.rate_constant,
        "E_high": D.energy_per_substance,
        "A_low": D.rate_constant,
        "E_low": D.energy_per_substance,
    }

    @unit_.manage_units([D.temperature, D.pressure], D.rate_constant)
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float128]:
        """Evaluate rate constant."""
        T_, P_ = func.normalize_arguments((T, P))  # noqa: N806
        P_r = self.effective_reduced_pressure(T_, P_)  # noqa: N806
        if self.activated:
            k_low, _ = self.arrhenius_functions
            kTP = k_low(T_) / (1 + P_r) * self.function(T_, P_r)  # noqa: N806
        else:
            _, k_high = self.arrhenius_functions
            kTP = k_high(T_) * P_r / (1 + P_r) * self.function(T_, P_r)  # noqa: N806
        return func.normalize_values(kTP, (T, P))

    @property
    def arrhenius_functions(
        self,
    ) -> tuple[ArrheniusRateFit, ArrheniusRateFit]:
        """Get low and high temperature arrhenius rate fits."""
        k_low = ArrheniusRateFit(
            A=self.A_high,
            b=self.b_high,
            E=self.E_high,
            order=self.order,
        )
        k_high = ArrheniusRateFit(
            A=self.A_high,
            b=self.b_high,
            E=self.E_high,
            order=self.order,
        )
        return k_low, k_high

    def effective_concentration(
        self,
        T: NDArray[np.float128],  # noqa: N803
        P: NDArray[np.float128],  # noqa: N803
    ) -> NDArray[np.float128]:
        """Get effective concentration(s) from temperature(s) and pressure(s).

        effective [M] = P / R T  (ideal gas law)

        :param T: Temperature(s)
        :param P: Pressure(s)
        :return: Effective concentration(s)
        """
        # Evaluate, using pint to handle units
        R_ = const.quantity(C.gas)  # noqa: N806
        T_ = pint.Quantity(T, UNITS.temperature)  # noqa: N806
        P_ = pint.Quantity(P, UNITS.pressure)  # noqa: N806
        m_ = P_ / (R_ * T_)

        # Return value in concentration units
        return m_.m_as(UNITS.concentration)

    def effective_reduced_pressure(
        self,
        T: NDArray[np.float128],  # noqa: N803
        P: NDArray[np.float128],  # noqa: N803
    ) -> NDArray[np.float128]:
        """Get effective concentration(s) from temperature(s) and pressure(s).

        effective P_r = k_low [M] / k_high  (ideal gas law)

        :param T: Temperature(s)
        :param P: Pressure(s)
        :return: Effective reduced pressure(s)
        """
        m_eff = self.effective_concentration(T, P)
        k_low, k_high = self.arrhenius_functions
        return k_low(T) * m_eff / k_high(T)


class PlogRateFit(RateFit):
    """Plog rate fit."""

    As: list[float]
    bs: list[float]
    Es: list[float]
    Ps: list[float]

    # Private attributes
    type_: ClassVar[str] = "plog"
    _scalers: ClassVar[Scalers] = {"As": np.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "As": D.rate_constant,
        "Es": D.energy_per_substance,
        "Ps": D.pressure,
    }

    @unit_.manage_units([D.temperature, D.pressure], D.rate_constant)
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float128]:
        """Evaluate rate constant for a single pressure."""
        T_, P_ = func.normalize_arguments((T, P))  # noqa: N806
        P0 = self.nearest_pressure(P_, which=0)  # noqa: N806
        P1 = self.nearest_pressure(P_, which=1)  # noqa: N806
        kT0 = self.nearest_arrhenius_values(T_, P_, which=0)  # noqa: N806
        kT1 = self.nearest_arrhenius_values(T_, P_, which=1)  # noqa: N806

        # Evaluate intermediate pressures
        log_P, log_P0, log_P1 = map(np.log, (P_, P0, P1))  # noqa: N806
        kTP = kT0 + (kT1 - kT0) * (log_P - log_P0) / (log_P1 - log_P0)  # noqa: N806

        # Evaluate on-boundary pressures (needed to fill in last pressure value)
        if np.ndim(P_) > 0:
            kTP[..., np.equal(P_, P0)] = kT0[..., np.equal(P_, P0)]
        elif P_ == P0:
            kTP = kT0  # noqa: N806

        return func.normalize_values(kTP, (T, P))

    @property
    def arrhenius_functions(self) -> list[ArrheniusRateFit]:
        """Get arrhenius rate fits in order."""
        return [
            ArrheniusRateFit(A=A, b=b, E=E, order=self.order)
            for A, b, E in zip(self.As, self.bs, self.Es, strict=True)
        ]

    @property
    def pressures(self) -> NDArray[np.float128]:
        """Pressures."""
        return np.array(self.Ps, dtype=np.float128)

    @property
    def pressure_indices(self) -> list[int]:
        """Pressure indices."""
        return list(range(self.pressures.shape[0]))

    def nearest_arrhenius_values(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike,  # noqa: N803
        which: int = 0,
    ) -> NDArray[np.float128]:
        """Get nearest lower or higher pressure.

        :param P: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest function(s)
        """
        iP = self.nearest_index(P, which=which)  # noqa: N806
        kTs = [k(T) for k in self.arrhenius_functions]  # noqa: N806
        return np.where(
            np.isin(iP, self.pressure_indices),
            np.choose(iP, kTs, mode="clip"),
            np.nan,
        )

    def nearest_pressure(
        self,
        P: ArrayLike,  # noqa: N803
        which: int = 0,
    ) -> NDArray[np.float128]:
        """Get nearest lower or higher pressure.

        :param P: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest pressure(s)
        """
        iP = self.nearest_index(P, which=which)  # noqa: N806
        ps = self.pressures
        return np.where(
            np.isin(iP, self.pressure_indices),
            np.take(ps, iP, mode="clip"),
            np.nan,
        )

    def nearest_index(
        self,
        P: ArrayLike,  # noqa: N803
        which: int = 0,
    ) -> NDArray[np.int_]:
        """Get nearest lower or higher index.

        :param P: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest index (indices)
        """
        return np.searchsorted(self.Ps, P, side="right") - 1 + which

    @classmethod
    @unit_.manage_units([D.temperature, D.rate_constant])
    def fit(  # noqa: PLR0913
        cls,
        Ts: Sequence[float],  # noqa: N803
        Ps: Sequence[float],  # noqa: N803
        k_data: ArrayLike,
        *,
        k_high: ArrayLike | None = None,
        A_fill: float | None = None,
        bad_fit: Literal["fill"]
        | Literal["warn"]
        | Literal["raise"]
        | Literal["ignore"] = "warn",
        bad_fit_fill_pressures: Sequence[float] = (),
        validate: bool = True,
        order: int = 1,
        units: UnitsData | None = None,
    ) -> "PlogRateFit":
        """Fit data to Plog rate fit.

        :param T: Temperatures
        :param P: Temperatures
        :param k_data: Finite pressure rates
        :param k_high: High-pressure-limit rates
        :param A_fill: Optional dummy parameter for unfittable rates
            (Otherwise, an error will be thrown.)
        :param bad_fit: How to handle bad fits that create floating point errors;
            (Options: "fill", replace with the fill value, or numpy.seterr options)
        :return: Rate fit
        """
        k_fits = []
        for k_data, P in zip(np.transpose(k_data), Ps, strict=True):
            if validate:
                print(f"Fitting Arrhenius rate for {P = }")

            bad_fit = (
                "fill" if np.any(np.isclose(P, bad_fit_fill_pressures)) else bad_fit
            )
            k_fit = ArrheniusRateFit.fit(
                Ts=Ts,
                ks=k_data,
                A_fill=A_fill,
                bad_fit=bad_fit,
                validate=validate,
                order=order,
                units=units,
            )
            k_fits.append(k_fit)

        k_high_fit = None
        if k_high is not None:
            k_high_fit = ArrheniusRateFit.fit(
                Ts=Ts, ks=k_high, validate=validate, order=order, units=units
            )
            msg = f"Currently not fitting high-pressure limit {k_high_fit}"
            warnings.warn(msg, stacklevel=2)

        return cls(
            order=order,
            As=[f.A for f in k_fits],
            bs=[f.b for f in k_fits],
            Es=[f.E for f in k_fits],
            Ps=Ps,  # pyright: ignore[reportArgumentType]
        )

    def is_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was cleared."""
        return all(k.is_cleared(A_fill=A_fill) for k in self.arrhenius_functions)

    def is_partially_cleared(self, A_fill: float) -> bool:
        """Determine if this rate was partially cleared."""
        return any(k.is_cleared(A_fill=A_fill) for k in self.arrhenius_functions)


class ChebRateFit(RateFit):
    """Chebyshev rate fit."""

    coeffs: NDArray_
    T_range: tuple[float, float]
    P_range: tuple[float, float]

    # Private attributes
    type_: ClassVar[str] = "cheb"
    _scalers: ClassVar[Scalers] = {"coeffs": np.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "coeffs": D.rate_constant,
        "T_range": D.temperature,
        "P_range": D.pressure,
    }

    @unit_.manage_units([D.temperature, D.pressure], D.rate_constant)
    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float128]:
        """Evaluate rate constant for a single pressure."""
        # Skip input processing, since chebgrid2d automatically forms the grid
        T0, T1 = self.T_range  # noqa: N806
        P0, P1 = self.P_range  # noqa: N806

        inv_ = np.reciprocal
        log_ = np.log10

        T_r = (2 * inv_(T) - inv_(T0) - inv_(T1)) / (inv_(T1) - inv_(T0))  # noqa: N806
        P_r = (2 * log_(P) - log_(P0) - log_(P1)) / (log_(P1) - log_(P0))  # noqa: N806

        # AVC: I don't understand why I need to transpose the coefficient matrix here
        kTP = chebyshev.chebgrid2d(T_r, P_r, self.coeffs.T)  # noqa: N806
        return func.normalize_values(kTP, (T, P))


Rate_ = Annotated[
    pydantic.SkipValidation[BaseRate],
    pydantic.BeforeValidator(lambda x: BaseRate.model_validate(x)),
    pydantic.PlainSerializer(lambda x: BaseRate.model_validate(x).model_dump()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(pydantic.BaseModel)),
    ),
]


# Conversions
def chemkin_string(rate_const: RateFit, eq_width: int = 0) -> str:
    """Write Chemkin rate to a string.

    :param rate_const: Rate constant
    :param eq_width: Equation width for alignment
    :return: Chemkin rate string
    """
    # Generate dummy head line
    head_line = chemkin.write_numbers([1.0, 0.0, 0.0])

    # Calculate the total width of the top line for alignment
    head_width = eq_width + len(head_line) + 1

    # Generate auxiliary lines and replace head line, if appropriate
    match rate_const:
        case ArrheniusRateFit():
            head_params = [rate_const.A, rate_const.b, rate_const.E]
            head_line = chemkin.write_numbers(head_params)
            aux_lines = []
        case FalloffRateFit(activated=False):
            high_params = [rate_const.A_high, rate_const.b_high, rate_const.E_high]
            low_params = [rate_const.A_low, rate_const.b_low, rate_const.E_low]
            head_line = chemkin.write_numbers(high_params)
            aux_lines = [chemkin.write_aux("LOW", low_params, head_width=head_width)]
        case FalloffRateFit(activated=True):
            high_params = [rate_const.A_high, rate_const.b_high, rate_const.E_high]
            low_params = [rate_const.A_low, rate_const.b_low, rate_const.E_low]
            head_line = chemkin.write_numbers(low_params)
            aux_lines = [chemkin.write_aux("HIGH", high_params, head_width=head_width)]
        case PlogRateFit():
            plog_params = [rate_const.Ps, rate_const.As, rate_const.bs, rate_const.Es]
            aux_lines = [
                chemkin.write_aux("PLOG", row) for row in zip(*plog_params, strict=True)
            ]
        case ChebRateFit():
            shape = np.shape(rate_const.coeffs)
            aux_lines = [
                chemkin.write_aux("TCHEB", rate_const.T_range, head_width=head_width),
                chemkin.write_aux("PCHEB", rate_const.P_range, head_width=head_width),
                chemkin.write_aux("CHEB", shape, head_width=head_width, as_int=True),
                *(
                    chemkin.write_aux("CHEB", row, head_width=head_width)
                    for row in rate_const.coeffs.tolist()
                ),
            ]
        case _:
            msg = f"Rate constant has unknown type {type(rate_const)}:\n{rate_const}"
            raise ValueError(msg)

    if isinstance(rate_const, FalloffRateFit):
        aux_lines.extend(
            blend.chemkin_aux_lines(rate_const.function, head_width=head_width),
        )

    eff_line = chemkin.write_efficiencies(
        rate_const.efficiencies,
        third_body=rate_const.third_body,
    )
    if eff_line is not None:
        aux_lines.append(eff_line)

    lines = [head_line, *aux_lines]
    return "\n".join(lines)


def from_mess_channel_output(mess_chan_out: str, order: int) -> Rate:
    """Extract rate data from MESS output.

    :param mess_chan_out: MESS output channel string
    :param order: Order
    :return: Rate data
    """
    res = mess.parse_output_channel(mess_chan_out)
    return from_mess_channel_output_parse_results(res, order=order)


# Parse helpers
def from_mess_channel_output_parse_results(
    res: mess.MessOutputChannelParseResults,
    order: int,
) -> Rate:
    """Extract rate data from MESS output parse results.

    :param res: MESS output parse results
    :param order: Order
    :return: Rate data
    """
    return Rate(
        order=order,
        T=res.T,
        P=res.P,
        k_data=np.transpose(res.k_data),
        k_high=res.k_high,
        units={"substance": "molec"},
    )


def from_chemkin_parse_results(
    res: chemkin.ChemkinRateParseResults,
    units: UnitsData | None = None,
) -> RateFit:
    """Extract rate data from Chemkin parse results.

    Chemkin parse results are modified in-place

    :param res: Chemkin rate parse results
    :return: Rate data
    """
    # Extract efficiencies
    efficiencies = res.efficiencies.copy()
    res.efficiencies.clear()

    # Determine reaction order
    order = len(res.reactants) + bool(efficiencies)

    # Determine units
    units = Units() if units is None else Units.model_validate(units)

    if "CHEB" in res.aux_numbers:
        # Read coefficients
        cheb: list[float] = res.aux_numbers.pop("CHEB")
        shape = tuple(map(int, cheb[:2]))
        coeffs = np.reshape(cheb[2:], shape)
        # Read ranges
        t_range = res.aux_numbers.pop("TCHEB")
        p_range = res.aux_numbers.pop("PCHEB")
        return ChebRateFit(
            coeffs=coeffs,
            T_range=t_range,
            P_range=p_range,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "PLOG" in res.aux_numbers:
        ps, As, bs, Es = zip(  # noqa: N806
            *mit.chunked(res.aux_numbers.pop("PLOG"), 4, strict=True),
            strict=True,
        )
        return PlogRateFit(
            As=As,
            bs=bs,
            Es=Es,
            Ps=ps,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "LOW" in res.aux_numbers:
        A_high, b_high, E_high = res.arrhenius  # noqa: N806
        A_low, b_low, E_low = res.aux_numbers.pop("LOW")  # noqa: N806
        function = blend.from_chemkin_parse_results(res)
        return FalloffRateFit(
            A_low=A_low,
            b_low=b_low,
            E_low=E_low,
            A_high=A_high,
            b_high=b_high,
            E_high=E_high,
            function=function,
            activated=False,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "HIGH" in res.aux_numbers:
        A_low, b_low, E_low = res.arrhenius  # noqa: N806
        A_high, b_high, E_high = res.aux_numbers.pop("HIGH")  # noqa: N806
        function = blend.from_chemkin_parse_results(res)
        return FalloffRateFit(
            A_low=A_low,
            b_low=b_low,
            E_low=E_low,
            A_high=A_high,
            b_high=b_high,
            E_high=E_high,
            function=function,
            activated=True,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    A, b, E = res.arrhenius  # noqa: N806
    return ArrheniusRateFit(
        A=A,
        b=b,
        E=E,
        efficiencies=efficiencies,
        order=order,
        units=units,
    )


# Display
def display(  # noqa: PLR0913
    rate_: BaseRate | Sequence[BaseRate],
    *,
    T_range: tuple[float, float] = (400, 1250),  # noqa: N803
    P: float = 1,  # noqa: N803
    units: UnitsData | None = None,
    label: str | Sequence[str] | None = None,
    color: str | Sequence[str] | None = None,
    x_label: str | None = None,  # noqa: RUF001
    y_label: str | None = None,  # noqa: RUF001
    x_unit: str | None = None,  # noqa: RUF001
    y_unit: str | None = None,  # noqa: RUF001
    check_order: bool = True,
    plot_type: Literal["arrh"] | Literal["simple"] = "arrh",
) -> alt.Chart:
    """Display one or more reaction rates on an Arrhenius plot.

    :param rxn_: Reaction rate(s)
    :param T_range: Temperature range, defaults to (400, 1250)
    :param P: Pressure
    :param label_: Label(s), defaults to None
    :param color_: Color(s), defaults to None
    :param x_label: X-axis label
    :param y_label: Y-axis label
    """
    rates = [rate_] if isinstance(rate_, BaseRate) else rate_
    labels = [label] if isinstance(label, str) else label
    colors = [color] if isinstance(color, str) else color
    rate0, *rates_ = rates
    if check_order:
        for other_rate in rates_:
            if not rate0.order == other_rate.order:
                msg = f"Mismatched reaction orders: {rate0} !~ {other_rate}"
                raise ValueError(msg)
    order = rate0.order

    plot_ = plot.arrhenius if plot_type == "arrh" else plot.simple

    def make_chart(
        ixs: Sequence[int],
        rates: Sequence[BaseRate],
        labels: Sequence[str] | None,
        colors: Sequence[str] | None,
        mark: str,
    ) -> alt.Chart:
        rates_ = [rates[i] for i in ixs]
        labels_ = None if labels is None else [labels[i] for i in ixs]
        colors_ = None if colors is None else [colors[i] for i in ixs]
        (T, *Ts), ks = zip(  # noqa: N806
            *(r.plot_data(T=T_range, P=P, units=units) for r in rates_),
            strict=True,
        )
        for T_ in Ts:  # noqa: N806
            assert np.allclose(T, T_), f"{T} !~ {T_}"
        return plot_(
            ks=ks,
            T=T,
            order=order,
            units=units,
            labels=labels_,
            colors=colors_,
            x_label=x_label,
            y_label=y_label,
            x_unit=x_unit,
            y_unit=y_unit,
            mark=mark,
        )

    charts = []
    for mark in (plot.Mark.line, plot.Mark.point):
        ixs = [i for i, r in enumerate(rates) if r.plot_mark == mark]
        if ixs:
            chart = make_chart(
                ixs, rates=rates, labels=labels, colors=colors, mark=mark
            )
            charts.append(chart)

    chart, *others = charts
    return (
        chart if not others else alt.layer(*charts).resolve_scale(color="independent")
    )


def display_p(  # noqa: PLR0913
    rate_: BaseRate | Sequence[BaseRate],
    *,
    T: float = 825,  # noqa: N803
    P_range: tuple[float, float] = (0.1, 100),  # noqa: N803
    units: UnitsData | None = None,
    label: str | Sequence[str] | None = None,
    color: str | Sequence[str] | None = None,
    y_label: str | None = None,  # noqa: RUF001
    y_unit: str | None = None,  # noqa: RUF001
    check_order: bool = True,
) -> alt.Chart:
    """Display one or more reaction rates on an Arrhenius plot.

    :param rxn_: Reaction rate(s)
    :param T_range: Temperature range, defaults to (400, 1250)
    :param P: Pressure
    :param label_: Label(s), defaults to None
    :param color_: Color(s), defaults to None
    :param x_label: X-axis label
    :param y_label: Y-axis label
    """
    rates = [rate_] if isinstance(rate_, BaseRate) else rate_
    labels = [label] if isinstance(label, str) else label
    colors = [color] if isinstance(color, str) else color
    rate0, *rates_ = rates
    if check_order:
        for other_rate in rates_:
            if not rate0.order == other_rate.order:
                msg = f"Mismatched reaction orders: {rate0} !~ {other_rate}"
                raise ValueError(msg)
    order = rate0.order

    units = UNITS if units is None else Units.model_validate(units)
    x_unit = unit_.pretty_string(units.pressure)
    x_label = f"𝑃 ({x_unit})"

    y_label = "𝑘" if y_label is None else y_label
    y_unit = (
        unit_.pretty_string(units.rate_constant(order)) if y_unit is None else y_unit
    )
    if y_unit:
        y_label = f"{y_label} ({y_unit})"

    def make_chart(
        ixs: Sequence[int],
        rates: Sequence[BaseRate],
        labels: Sequence[str] | None,
        colors: Sequence[str] | None,
        mark: str,
    ) -> alt.Chart:
        rates_ = [rates[i] for i in ixs]
        labels_ = None if labels is None else [labels[i] for i in ixs]
        colors_ = None if colors is None else [colors[i] for i in ixs]
        (P, *Ps), ks = zip(  # noqa: N806
            *(r.plot_data(T=T, P=P_range, units=units) for r in rates_),
            strict=True,
        )
        for P_ in Ps:  # noqa: N806
            assert np.allclose(P, P_), f"{P} !~ {P_}"
        return plot.general(
            y_data=ks,
            x_data=P,
            labels=labels_,
            colors=colors_,
            x_label=x_label,
            y_label=y_label,
            x_scale=plot.log_scale(P_range),
            x_axis=plot.log_scale_axis(P_range),
            mark=mark,
        )

    charts = []
    for mark in (plot.Mark.line, plot.Mark.point):
        ixs = [i for i, r in enumerate(rates) if r.plot_mark == mark]
        if ixs:
            chart = make_chart(
                ixs, rates=rates, labels=labels, colors=colors, mark=mark
            )
            charts.append(chart)

    chart, *others = charts
    return (
        chart if not others else alt.layer(*charts).resolve_scale(color="independent")
    )
