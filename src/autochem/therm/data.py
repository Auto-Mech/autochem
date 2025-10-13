"""Thermodynamic data."""

import abc
import functools
import itertools
from collections.abc import Sequence
from typing import Annotated, ClassVar, Literal, Self

import altair as alt
import numpy as np
import pandas as pd
import pint
import pydantic
import xarray
from numpy.typing import ArrayLike, NDArray
from pydantic_core import core_schema

from .. import unit_
from ..unit_ import UNITS, C, D, Dimension, UnitManager, Units, UnitsData, dim
from ..util import FormulaData, chemkin, form, pac99, plot
from ..util.type_ import Formula_, Frozen, Scalable, Scalers, SubclassTyped
from .func import Bounded, Nasa7Calculator, ThermCalculator


class Key:
    """Property key values."""

    # Independent variables
    T = "T"
    # Thermodynamic functions
    Cv = "Cv"
    Cp = "Cp"
    S = "S"
    H = "H"
    dH = "dH"  # noqa: N815
    # Partition function and derivatives
    Z0 = "Z0"
    Z1 = "Z1"
    Z2 = "Z2"


class BaseTherm(ThermCalculator, UnitManager, Frozen, Scalable, SubclassTyped, abc.ABC):
    """Abstract base class for thermodynamic data."""

    formula: Formula_
    charge: int = 0

    def display(  # noqa: PLR0913
        self,
        props: Sequence[Literal["Cv", "Cp", "S", "H", "dH"]] = ("Cp", "S", "H"),
        *,
        others: "Sequence[BaseTherm]" = (),
        others_labels: Sequence[str] = (),
        T_range: tuple[float, float] = (200, 3000),  # noqa: N803
        units: UnitsData | None = None,
        label: str = "This work",
        x_label: str = "𝑇",  # noqa: RUF001
        y_labels: Sequence[str | None] | None = None,
        horizontal: bool = False,
    ) -> alt.Chart:
        """Display as a thermodynamic function plot.

        :param props: Thermodynamic properties to display
        :param others: Other thermodynamic data to compare to
        :param others_labels: Labels for other thermodynamic data
        :param T_range: Temperature range
        :param units: Units
        :param x_label: X-axis label
        :param y_labels: Y-axis labels, by property
        :param horizontal: Whether to display horizontally
        :return: Chart
        """
        y_labels = y_labels or [None] * len(props)
        charts = [
            self._display(
                prop=prop,
                others=others,
                others_labels=others_labels,
                T_range=T_range,
                units=units,
                label=label,
                x_label=x_label,
                y_label=y_label,
            )
            for prop, y_label in zip(props, y_labels, strict=True)
        ]
        concat_ = alt.hconcat if horizontal else alt.vconcat
        return concat_(*charts)

    def _display(  # noqa: PLR0913
        self,
        prop: Literal["Cv", "Cp", "S", "H", "dH"],
        others: "Sequence[BaseTherm]" = (),
        others_labels: Sequence[str] = (),
        T_range: tuple[float, float] = (200, 3000),  # noqa: N803
        units: UnitsData | None = None,
        label: str = "This work",
        x_label: str = "𝑇",  # noqa: RUF001
        y_label: str | None = None,
    ) -> alt.Chart:
        """Display as a thermodynamic function plot.

        :param prop: Thermodynamic properties to display
        :param others: Other thermodynamic data to compare to
        :param others_labels: Labels for other thermodynamic data
        :param T_range: Temperature range
        :param units: Units
        :param x_label: X-axis label
        :param y_labels: Y-axis labels, by property
        :return: Chart
        """
        units = UNITS if units is None else Units.model_validate(units)

        # Property units
        prop_unit_dct = {
            Key.Cv: units.energy_per_substance / units.temperature,
            Key.Cp: units.energy_per_substance / units.temperature,
            Key.S: units.energy_per_substance / units.temperature,
            Key.H: units.energy_per_substance,
            Key.dH: units.energy_per_substance,
        }
        prop_func_dct = {
            Key.Cv: heat_capacity_constant_volume,
            Key.Cp: heat_capacity_constant_pressure,
            Key.S: entropy,
            Key.H: enthalpy,
            Key.dH: delta_enthalpy,
        }
        prop_label_dct = {
            Key.Cv: "𝐶ᵥ",
            Key.Cp: "𝐶ₚ",
            Key.S: "𝑆",  # noqa: RUF001
            Key.H: "𝐻",  # noqa: RUF001
            Key.dH: "Δ𝐻",
        }

        # Process units
        x_unit = unit_.pretty_string(units.temperature)
        y_unit = unit_.pretty_string(prop_unit_dct.get(prop))

        # Add units to labels
        x_label = f"{x_label} ({x_unit})"
        y_label = f"{y_label or prop_label_dct.get(prop)} ({y_unit})"

        # Get property function
        func_ = prop_func_dct[prop]

        # Gather objects and labels
        assert len(others) == len(others_labels), f"{others_labels} !~ {others}"
        all_objs = [self, *others]
        all_labels = [label, *others_labels]
        all_colors = plot.LINE_COLOR_CYCLE[: len(all_labels)]

        # Gather data from functons
        T = np.linspace(*T_range, num=500)  # noqa: N806
        data_dct = {L: func_(o, T) for L, o in zip(all_labels, all_objs, strict=True)}
        data = pd.DataFrame({"x": T, **data_dct})

        # Prepare encoding parameters
        x = alt.X("x", title=x_label)
        y = alt.Y("value:Q", title=y_label)
        color = (
            alt.Color(
                "key:N",
                scale=alt.Scale(domain=all_labels, range=all_colors),
            )
            if others
            else alt.value(all_colors[0])
        )

        # Create chart
        return (
            alt.Chart(data)
            .mark_line()
            .transform_fold(fold=list(data_dct.keys()))
            .encode(x=x, y=y, color=color)
        )


class Therm(BaseTherm):
    """Raw thermodynamic data.

    :param T: Temperatures (K)
    :param Z0: Logs of one-particle partition function per volume, ln(Q_1' [cm^-3])
    :param Z1: First temperature derivatives of Z0, d(ln(Q_1'))/dT [K^-1]
    :param Z2: Second temperature derivatives of Z0, d^2(ln(Q_1'))/dT^2 [K^-2]
    :param Hf: Enthalpy of formation at 298.15 K [cal/mol]
    """

    T: list[float]
    Z0: list[float]
    Z1: list[float]
    Z2: list[float]
    Hf: float | None = None

    # Private attributes
    type_: ClassVar[str] = "data"
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "T": D.temperature,
        "Z0": dim.log(D.volume),
        "Z1": D.temperature**-1,
        "Z2": D.temperature**-2,
    }

    @property
    def T_min(self) -> float:  # noqa: N802
        """Get minimum temperature."""
        return min(self.T)

    @property
    def T_max(self) -> float:  # noqa: N802
        """Get maximum temperature."""
        return max(self.T)

    @pydantic.model_validator(mode="after")
    def sort_by_temperature(self) -> Self:
        """Sort arrays by temperature."""
        frozen = self.model_config["frozen"]
        self.model_config["frozen"] = False
        self.T, self.Z0, self.Z1, self.Z2 = map(
            list,
            zip(
                *sorted(zip(self.T, self.Z0, self.Z1, self.Z2, strict=True)),
                strict=True,
            ),
        )
        self.model_config["frozen"] = frozen
        return self

    # Replace this with a model validator (before or after), allowing extra arguments
    def __init__(
        self,
        formula: dict[str, int],
        Hf: float | None = None,  # noqa: N803
        Tf: float = 298.15,  # noqa: N803
        **kwargs: object,
    ) -> None:
        """Initialize, handling enthalpy of formation.

        :param Hf: Enthalpy of formation at 0 K or 298.15 K
        :param Tf: Reference temperature for enthalpy of formation, 0 K or 298.15 K
        """
        # Initialize the object first, to allow enthalpy calculation if needed
        super().__init__(formula=formula, Hf=Hf, **kwargs)

        if int(Tf) == 0:
            # Access through __dict__ here since model is frozen
            self.__dict__["Hf"] += self.delta_enthalpy(
                T=298,
                method="nearest",
            ) - elemental_delta_enthalpy_room_temperature(self.formula)
        elif int(Tf) != 298:
            msg = f"Invalid reference temperature Tf = {Tf}"
            raise ValueError(msg)

    # Properties
    @unit_.manage_units([], D.energy_per_substance)
    def enthalpy_of_formation(self, units: UnitsData | None = None) -> float:  # noqa: ARG002
        """Get enthalpy of formation at 298.15 K."""
        assert self.Hf is not None, "Enthalpy of formation not set"
        return self.Hf

    @property
    def data_set(self) -> xarray.Dataset:
        """Access data as an xarray Dataset."""
        coord_key = Key.T
        coord_vals = self.T
        data_arrs = {
            Key.Z0: self.Z0,
            Key.Z1: self.Z1,
            Key.Z2: self.Z2,
            Key.Cv: self.heat_capacity_data(const="V"),
            Key.Cp: self.heat_capacity_data(const="P"),
            Key.S: self.entropy_data(),
            Key.dH: self.delta_enthalpy_data(),
        }
        return xarray.Dataset(
            data_vars={k: ([coord_key], v) for k, v in data_arrs.items()},
            coords={coord_key: coord_vals},
        )

    def __call__(self, T: ArrayLike) -> NDArray[np.float64]:
        """Evaluate partition function."""
        Z0: NDArray[np.float64] = (
            self.data_set["Z0"]
            .sel(  # noqa: N806
                {Key.T: T},
                method="ffill",
            )
            .data
        )
        return np.exp(Z0)

    # Thermodynamic function data points
    def temperature_data(self) -> NDArray[np.float64]:
        """Get temperature data."""
        return np.array(self.T, dtype=np.float64)

    @unit_.manage_units([], D.energy_per_substance / D.temperature)
    def heat_capacity_data(
        self,
        const: Literal["P", "V"] = "P",
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Calculate the heat capacity at constant volume or pressure.

        Formula:

            C_v = R (2 T d(ln(Q_1'))/dT + T^2 d^2(ln(Q_1'))/dT^2)
                = R (2 T Z_1 + T^2 Z_2)
            C_p = C_v + R = R (1 + 2 T Z_1 + T^2 Z_2)

        :param const: Whether to hold pressure ("P") or volume ("V") constant
        :param units: Units
        :return: Heat capacity
        """
        # Evaluate
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(self.T, dtype=np.float64)  # noqa: N806
        Z1 = np.array(self.Z1, dtype=np.float64)  # noqa: N806
        Z2 = np.array(self.Z2, dtype=np.float64)  # noqa: N806
        C_ = R * (1 + 2 * T * Z1 + T**2 * Z2)  # noqa: N806
        C_ -= R if const == "V" else 0.0  # noqa: N806
        return C_

    # TODO(avcopan): Fix unit manager to handle pressure keyword input  # noqa: FIX002
    # https://github.com/Auto-Mech/autochem/issues/670
    @unit_.manage_units([], D.energy_per_substance / D.temperature)
    def entropy_data(
        self,
        P: float = 1,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Calculatate entropy.

        Formula:

            S = R (T d(ln(Q_1'))/dT + ln(Q_1') - ln(c) + 1)
              = R (T Z1 + Z0 - ln(c) + 1)

        in terms of the one-particle partition function per unit volume, Q_1', and the
        standard concentration, c = N_A / V, which comes from accounting for volume.
        The latter is determined for a standard state pressure via the ideal gas law:

           c = P / k_B T

        This must be converted to internal units.

        :param P: Pressure
        :param units: Units
        :return: Entropy
        """
        T = np.array(self.T, dtype=np.float64)  # noqa: N806

        # Evaluate standard concentration (molecules* / volume) from pressure  *implicit
        units = Units.model_validate(units) if units is not None else UNITS
        k_B_ = pint.Quantity("boltzmann_constant")  # noqa: N806
        T_ = pint.Quantity(T, UNITS.temperature)  # noqa: N806
        P_ = pint.Quantity(P, units.pressure)  # noqa: N806
        c = (P_ / (k_B_ * T_)).m_as("1/cm**3")

        # Evaluate
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(self.T, dtype=np.float64)  # noqa: N806
        Z0 = np.array(self.Z0, dtype=np.float64)  # noqa: N806
        Z1 = np.array(self.Z1, dtype=np.float64)  # noqa: N806
        return R * (T * Z1 + Z0 - np.log(c) + 1)

    @unit_.manage_units([], D.energy_per_substance)
    def enthalpy_data(self, units: UnitsData | None = None) -> NDArray[np.float64]:  # noqa: ARG002
        """Calculate enthalpy data.

        Formula:

            H(T) = Hf(298.15 K) + (dH(T) -  dH(298.15 K))

        :param units: Units
        :return: Enthalpy
        """
        return self.Hf + self.delta_enthalpy_data() - self.delta_enthalpy(T=298.15)

    @unit_.manage_units([], D.energy_per_substance)
    def delta_enthalpy_data(
        self,
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Calculate delta enthalpy data.

        Formula:

            dH = R (T^2 d(ln(Q_1'))/dT + T)
               = R (T^2 Z_1 + T)

        :param units: Units
        :return: Enthalpy
        """
        # Evaluate
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(self.T, dtype=np.float64)  # noqa: N806
        Z1 = np.array(self.Z1, dtype=np.float64)  # noqa: N806
        return R * (T**2 * Z1 + T)

    # Thermodynamic function calculators
    def heat_capacity(
        self,
        T: ArrayLike,  # noqa: N803
        const: Literal["P", "V"] = "P",
        units: UnitsData | None = None,  # noqa: ARG002
        method: str = "nearest",
    ) -> NDArray[np.float64]:
        """Evaluate heat capacity, C_V(T) or C_P(T).

        :param T: Temperature(s)
        :param const: Whether to hold pressure ("P") or volume ("V") constant
        :param method: Xarray data selection method
        :param units: Unit system
        :return: Function value(s)
        """
        key = Key.Cp if const == "P" else Key.Cv
        return self.data_set[key].sel(T=T, method=method).data

    def entropy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
        method: str = "nearest",
    ) -> NDArray[np.float64]:
        """Evaluate entropy, S(T).

        :param T: Temperature(s)
        :param method: Xarray data selection method
        :param units: Unit system
        :return: Function value(s)
        """
        return self.data_set[Key.S].sel(T=T, method=method).data

    @unit_.manage_units([], D.energy_per_substance)
    def enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
        method: str = "nearest",
    ) -> NDArray[np.float64]:
        """Calculate enthalpy data.

        Formula:

            H(T) = Hf(298.15 K) + (dH(T) -  dH(298.15 K))

        :param T: Temperature for evaluation
        :param method: Xarray data selection method
        :param units: Units
        :return: Enthalpy
        """
        return (
            self.Hf
            + self.delta_enthalpy(T=T, method=method)
            - self.delta_enthalpy(T=298.15, method=method)
        )

    @unit_.manage_units([], D.energy_per_substance)
    def delta_enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
        method: str = "nearest",
    ) -> NDArray[np.float64]:
        """Calculate enthalpy.

        Formula:

            dH = R (T^2 d(ln(Q_1'))/dT + T)
               = R (T^2 Z_1 + T)

        :param T: Temperature for evaluation
        :param method: Xarray data selection method
        :param units: Units
        :return: Enthalpy
        """
        return self.data_set[Key.dH].sel(T=T, method=method).data


class ThermFit(BaseTherm, Bounded, abc.ABC):
    """Fitted thermodynamic data."""


class Nasa7ThermFit(ThermFit):
    """Fitted thermodynamic data."""

    T_mid: float
    coeffs_low: list[float]
    coeffs_high: list[float]

    # Private attributes
    type_: ClassVar[str] = "nasa7"
    _scalers: ClassVar[Scalers] = {
        "coeffs_low": np.multiply,
        "coeffs_high": np.multiply,
    }

    def piecewise_calculators(self) -> list[Nasa7Calculator]:
        """Get calculators for a numpy.piecewise evaluation.

        :return: Pair of calculators
        """
        calc_low = Nasa7Calculator.from_coefficients(
            T_min=self.T_min,
            T_max=self.T_mid,
            coeffs=self.coeffs_low,
        )
        calc_high = Nasa7Calculator.from_coefficients(
            T_min=self.T_mid,
            T_max=self.T_max,
            coeffs=self.coeffs_high,
        )
        return [calc_low, calc_high]

    def piecewise_conditions(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> list[NDArray[np.bool_]]:
        """Get conditions for a numpy.piecewise evaluation.

        :param T: Temperature(s)
        :return: Pair of boolean arrays
        """
        calc_low, calc_high = self.piecewise_calculators()
        return [
            calc_low.in_bounds(T, include_max=False),
            calc_high.in_bounds(T, include_max=True),
        ]

    def heat_capacity(
        self,
        T: ArrayLike,  # noqa: N803
        const: Literal["P", "V"] = "P",
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate heat capacity, C_V(T) or C_P(T).

        :param T: Temperature(s)
        :param const: Whether to hold pressure ("P") or volume ("V") constant
        :param units: Unit system
        :return: Function value(s)
        """
        conds = self.piecewise_conditions(T)
        funcs = [
            functools.partial(calc.heat_capacity, const=const)
            for calc in self.piecewise_calculators()
        ]
        return np.piecewise(T, conds, funcs)

    def entropy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate entropy, S(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        conds = self.piecewise_conditions(T)
        funcs = [calc.entropy for calc in self.piecewise_calculators()]
        return np.piecewise(T, conds, funcs)

    def enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy, H(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        conds = self.piecewise_conditions(T)
        funcs = [calc.enthalpy for calc in self.piecewise_calculators()]
        return np.piecewise(T, conds, funcs)

    def delta_enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy change, dH(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        conds = self.piecewise_conditions(T)
        funcs = [calc.delta_enthalpy for calc in self.piecewise_calculators()]
        return np.piecewise(T, conds, funcs)

    @classmethod
    def fit(  # noqa: PLR0913
        cls,
        T: ArrayLike,  # noqa: N803
        Cp: ArrayLike,  # noqa: N803
        S: ArrayLike,  # noqa: N803
        H: ArrayLike,  # noqa: N803
        formula: FormulaData,
        charge: int = 0,
        T_mid: float = 1000,  # noqa: N803
        T_min: float | None = None,  # noqa: N803
        T_max: float | None = None,  # noqa: N803
    ) -> "Nasa7ThermFit":
        """Fit data to Nasa-7 therm fit object.

        Construct a matrix mapping (unknown) polynomial coefficients to property values.

            Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
            S/R = a1 ln(T) + a2 T + (a3/2) T^2 + (a4/3) T^3 + (a5/4) T^4 + a7
            H/R = a1 T + (a2/2) T^2 + (a3/3) T^3 + (a4/4) T^4 + (a5/5) T^5 + a6


        The matrix is solved using least squares to find the polynomial coefficients:

            [     1,    T1,      T1^2,      T1^3,      T1^4, 0, 0]        [Cp(T1)/R]
            [     1,    T2,      T2^2,      T2^3,      T2^4, 0, 0] [a1]   [Cp(T2)/R]
            [   ...,   ...,       ...,       ...,       ..., 0, 0] [a2]   [   ...  ]
            [ln(T1),    T1, (1/2)T1^2, (1/3)T1^3, (1/4)T1^4, 0, 1] [a3]   [ S(T1)/R]
            [ln(T2),    T2, (1/2)T2^2, (1/3)T2^3, (1/4)T2^4, 0, 1] [a4] = [ S(T2)/R]
            [   ...,   ...,       ...,       ...,       ..., 0, 1] [a5]   [   ...  ]
            [    T1,  T1^2, (1/3)T1^3, (1/4)T1^4, (1/5)T1^5, 1, 0] [a6]   [ H(T1)/R]
            [    T2,  T2^2, (1/3)T2^3, (1/4)T2^4, (1/5)T2^5, 1, 0] [a7]   [ H(T2)/R]
            [   ...,   ...,       ...,       ...,       ..., 1, 0]        [   ...  ]

            Dimension: 3*n_T x 7    (n_T = number of temperature values)

        :param T: Temperatures
        :param Cp: Constant-pressure heat capacities
        :param S: Entropies
        :param H: Enthalpies
        :param T_mid: Middle temperature
        :param T_min: Minimum temperature
        :param T_max: Maximum temperature
        :return: Fitted object
        """
        T = np.array(T, dtype=np.float64)  # noqa: N806
        Cp = np.array(Cp, dtype=np.float64)  # noqa: N806
        S = np.array(S, dtype=np.float64)  # noqa: N806
        H = np.array(H, dtype=np.float64)  # noqa: N806
        T_min = T_min or np.min(T)  # noqa: N806
        T_max = T_max or np.max(T)  # noqa: N806
        low = (T_min <= T) & (T_mid >= T)
        high = (T_mid <= T) & (T_max >= T)
        calc_low = Nasa7Calculator.fit(T=T[low], Cp=Cp[low], S=S[low], H=H[low])
        calc_high = Nasa7Calculator.fit(T=T[high], Cp=Cp[high], S=S[high], H=H[high])
        return cls(
            formula=formula,
            charge=charge,
            T_min=calc_low.T_min,
            T_mid=T_mid,
            T_max=calc_high.T_max,
            coeffs_low=calc_low.coefficients,
            coeffs_high=calc_high.coefficients,
        )


Therm_ = Annotated[
    pydantic.SkipValidation[BaseTherm],
    pydantic.BeforeValidator(lambda x: BaseTherm.model_validate(x)),
    pydantic.PlainSerializer(lambda x: BaseTherm.model_validate(x).model_dump()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(pydantic.BaseModel)),
    ),
]


# Properties
def heat_capacity(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    const: Literal["P", "V"] = "P",
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate heat capacity, Cv(T) or Cp(T).

    :param T: Temperature(s)
    :param const: Whether to hold pressure ("P") or volume ("V") constant
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.heat_capacity(T=T, const=const, units=units)


def heat_capacity_constant_pressure(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate heat capacity at constant pressure, Cp(T).

    :param T: Temperature(s)
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.heat_capacity_constant_pressure(T=T, units=units)


def heat_capacity_constant_volume(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate heat capacity at constant volume, Cv(T).

    :param T: Temperature(s)
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.heat_capacity_constant_volume(T=T, units=units)


def entropy(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate entropy, S(T).

    :param T: Temperature(s)
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.entropy(T=T, units=units)


def enthalpy(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate enthalpy, H(T).

    :param T: Temperature(s)
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.enthalpy(T=T, units=units)


def delta_enthalpy(
    therm: ThermCalculator,
    T: ArrayLike,  # noqa: N803
    units: UnitsData | None = None,
) -> NDArray[np.float64]:
    """Evaluate enthalpy change, dH(T).

    :param T: Temperature(s)
    :param units: Unit system
    :return: Function value(s)
    """
    return therm.delta_enthalpy(T=T, units=units)


# Conversions
def from_messpf_output_string(  # noqa: PLR0913
    pf_str: str,
    formula: FormulaData,
    charge: int = 0,
    Hf: float | None = None,  # noqa: N803
    Tf: float = 0,  # noqa: N803
    units: UnitsData | None = None,
) -> Therm:
    """Extract thermo data from MESS-PF output string.

    :param pf_str: MESS-PF output string
    :param formula: Formula
    :param charge: Molecular charge
    :param Hf: Enthalpy of formation at 0 K or 298 K
    :param Tf: Enthalpy of formation reference temperature, 0 K or 298 K
    :param units: Units of enthalpy of formation
    :return: Thermo data
    """
    units0 = UNITS if units is None else Units.model_validate(units)
    Hf = None if Hf is None else dim.convert(units0, UNITS, D.energy, Hf)  # noqa: N806
    formula = form.normalize_input(formula)
    lines = list(map(str.strip, pf_str.strip().splitlines()))
    lines = list(itertools.dropwhile(lambda s: not s.startswith("Z_0"), lines))
    data = [list(map(float, line.split())) for line in lines[1:]]
    T, Z0, Z1, Z2, *_ = map(list, zip(*data, strict=True))  # noqa: N806
    # Replace 298.2 with 298.15 (needed for PAC99 input to run)
    T = [298.15 if round(T) == 298 else T for T in T]  # noqa: N806
    return Therm(
        T=T,
        Z0=Z0,
        Z1=Z1,
        Z2=Z2,
        Hf=Hf,
        Tf=Tf,
        formula=formula,
        charge=charge,
    )


def from_chemkin_string(
    therm_str: str,
    T_mid: float = 1000,  # noqa: N803
) -> ThermFit:
    """Read species thermo from Chemkin string.

    :param therm_str: Chemkin species therm string
    :param T_mid: Default mid-point temperature
    :return: Species thermo
    """
    # Parse string
    res = chemkin.parse_thermo(therm_str)

    # Extract thermo data
    return from_chemkin_parse_results(res, T_mid=T_mid)


def from_chemkin_parse_results(
    res: chemkin.ChemkinThermoParseResults,
    T_mid: float | None = None,  # noqa: N803
) -> ThermFit:
    """Extract thermo data from Chemkin parse results.

    :param res: Chemkin thermo parse results
    :param T_mid: Default mid-point temperature
    :return: Thermo data fit
    """
    T_mid = T_mid or 1000.0  # noqa: N806

    # Determine charge, if any
    charge = 0
    formula = res.formula.copy()
    if "E" in formula:
        charge = -formula.pop("E")

    # Read in coefficients
    if len(res.coeffs) == 14:
        return Nasa7ThermFit(
            T_min=res.T_min,
            T_max=res.T_max,
            coeffs_low=res.coeffs[7:],
            coeffs_high=res.coeffs[:7],
            T_mid=res.T_mid or T_mid,
            formula=formula,
            charge=charge,
        )

    msg = f"Unable to interpret parse results: {res}"
    raise ValueError(msg)


def from_pac99_output_string(therm_str: str) -> Nasa7ThermFit:
    """Read species thermo from PAC99 output string.

    :param therm_str: PAC99 .c97 output string
    :return: NASA7 thermo fit
    """
    # Parse string
    res = pac99.parse_thermo(therm_str)

    # Extract thermo data
    return from_pac99_output_parse_results(res)


def from_pac99_output_parse_results(
    res: pac99.Pac99ThermoParseResults,
) -> Nasa7ThermFit:
    """Extract thermo data from Chemkin parse results.

    :param res: Chemkin thermo parse results
    :param T_mid: Default mid-point temperature
    :return: Thermo data fit
    """
    # Determine charge, if any
    charge = 0
    formula = res.formula.copy()
    if "E" in formula:
        charge = -formula.pop("E")

    assert len(res.ranges) == len(res.coeffs_lst) == 2, res

    coeffs_low, coeffs_high = res.coeffs_lst
    (T_min, T_mid), (T_mid_, T_max) = res.ranges  # noqa: N806
    assert T_mid == T_mid_, f"{T_mid} != {T_mid_}"

    # Read in coefficients
    return Nasa7ThermFit(
        T_min=T_min,
        T_max=T_max,
        coeffs_low=coeffs_low,
        coeffs_high=coeffs_high,
        T_mid=T_mid,
        formula=formula,
        charge=charge,
    )


# Helpers
KJ_TO_CAL = pint.Quantity(1, "kJ").m_as("cal")
ENTHALPY_CHANGE_0K_TO_298K: dict[str, float] = {
    # Monatomic values (kJ)
    "He": 6.197,
    "Ne": 6.197,
    "Ar": 6.197,
    "Kr": 6.197,
    "Xe": 6.197,
    "S": 4.412,
    "P": 5.36,
    "C": 1.05,
    "Si": 3.217,
    "Ge": 4.636,
    "Sn": 6.323,
    "Pb": 6.87,
    "B": 1.222,
    "Al": 4.54,
    "Zn": 5.647,
    "Cd": 6.247,
    "Hg": 9.342,
    "Cu": 5.004,
    "Ag": 5.745,
    "Ti": 4.824,
    "U": 6.364,
    "Th": 6.35,
    "Be": 1.95,
    "Mg": 4.998,
    "Ca": 5.736,
    "Li": 4.632,
    "Na": 6.46,
    "K": 7.088,
    "Rb": 7.489,
    "Cs": 7.711,
    # Diatomic values (kJ)
    "O": 8.68 / 2,
    "H": 8.468 / 2,
    "F": 8.825 / 2,
    "Cl": 9.181 / 2,
    "Br": 24.52 / 2,
    "N": 8.67 / 2,
}
ENTHALPY_CHANGE_0K_TO_298K = {
    k: v * KJ_TO_CAL for k, v in ENTHALPY_CHANGE_0K_TO_298K.items()
}


def elemental_delta_enthalpy_room_temperature(
    formula: FormulaData,
) -> float:
    """Get change in enthalpy of formation 0 K to room temperature (298 K).

    :param formula: Formula
    """
    formula = form.normalize_input(formula)
    assert all(k in ENTHALPY_CHANGE_0K_TO_298K for k in formula), (
        f"Invalid symbols: {formula}"
    )
    return sum(ENTHALPY_CHANGE_0K_TO_298K[k] * v for k, v in formula.items())
