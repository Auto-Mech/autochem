"""Single-temperature-range thermodynamic fit function models."""

import abc
from typing import ClassVar, Literal

import numpy as np
import pydantic
from numpy.typing import ArrayLike, NDArray

from .. import unit_
from ..unit_ import UNITS, C, D, UnitsData
from ..util.type_ import Frozen


class Bounded(pydantic.BaseModel):
    """Mixin to define bounded calculator."""

    T_min: float
    T_max: float

    def in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
        *,
        include_max: bool = True,
    ) -> NDArray[np.bool_]:
        """Determine whether temperature(s) are in bounds.

        :param T: Temperature(s)
        :return: Boolean value(s)
        """
        T = np.array(T, dtype=np.float64)  # noqa: N806
        greater_than_min = self.T_min <= T
        less_than_max = (self.T_max >= T) if include_max else (self.T_max > T)
        return greater_than_min & less_than_max

    def all_in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> bool:
        """Determine whether all temperature(s) are in bounds.

        :param T: Temperature(s)
        :return: `True` if they are
        """
        return np.all(self.in_bounds(T))

    def assert_all_in_bounds(
        self,
        T: ArrayLike,  # noqa: N803
    ) -> None:
        """Assert that all temperature(s) are in bounds.

        :param T: Temperature(s)
        """
        assert self.all_in_bounds(T), f"{self.T_min} !<= {T} !<= {self.T_max}"


class ThermCalculator(Frozen, abc.ABC):
    """Abstract base class for therm calculators."""

    @abc.abstractmethod
    def heat_capacity(
        self,
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

    def heat_capacity_constant_pressure(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate heat capacity at constant pressure, Cp(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        return self.heat_capacity(T, const="P", units=units)

    def heat_capacity_constant_volume(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate heat capacity at constant volume, Cv(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        return self.heat_capacity(T, const="V", units=units)

    @abc.abstractmethod
    def entropy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate entropy, S(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """

    @abc.abstractmethod
    def enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy, H(T).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """

    @abc.abstractmethod
    def thermal_enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate thermal enthalpy, H(T) - H(0).

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """


class Nasa7Calculator(Bounded, ThermCalculator):
    """Nasa Polynomial calculator (7 coefficients)."""

    a1: float
    a2: float
    a3: float
    a4: float
    a5: float
    a6: float
    a7: float

    # Private attributes
    type_: ClassVar[str] = "nasa7"

    @property
    def coefficients(self) -> list[float]:
        """Get coefficients of the calculator."""
        return [self.a1, self.a2, self.a3, self.a4, self.a5, self.a6, self.a7]

    @classmethod
    def from_coefficients(
        cls,
        T_min: float,  # noqa: N803
        T_max: float,  # noqa: N803
        coeffs: list[float],
    ) -> "Nasa7Calculator":
        """Create a Nasa7Calculator from coefficients."""
        a1, a2, a3, a4, a5, a6, a7 = coeffs
        return cls(
            T_min=T_min,
            T_max=T_max,
            a1=a1,
            a2=a2,
            a3=a3,
            a4=a4,
            a5=a5,
            a6=a6,
            a7=a7,
        )

    @unit_.manage_units([], D.energy_per_substance / D.temperature)
    def heat_capacity(
        self,
        T: ArrayLike,  # noqa: N803
        const: Literal["P", "V"] = "P",
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate heat capacity, Cv(T) or Cp(T).

        Formula:
            Cp(T) = R (a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4)
            Cv(T) = Cp(T) - R

        :param T: Temperature(s)
        :param const: Whether to hold pressure ("P") or volume ("V") constant
        :param units: Unit system
        :return: Function value(s)
        """
        self.assert_all_in_bounds(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(T, dtype=np.float64)  # noqa: N806
        C_ = R * (  # noqa: N806
            self.a1 + self.a2 * T + self.a3 * T**2 + self.a4 * T**3 + self.a5 * T**4
        )
        C_ -= R if const == "V" else 0.0  # noqa: N806
        return C_

    @unit_.manage_units([], D.energy_per_substance / D.temperature)
    def entropy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate entropy, S(T).

        Formula:
            S(T) = R (a1 ln(T) + a2 T + (a3/2) T^2 + (a4/3) T^3 + (a5/4) T^4 + a7)

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        self.assert_all_in_bounds(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        return R * self.a7 + self.thermal_entropy(T)

    @unit_.manage_units([], D.energy_per_substance)
    def enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy, H(T).

        Formula:
            H(T) = R (a1 T + (a2/2) T^2 + (a3/3) T^3 + (a4/4) T^4 + (a5/5) T^5 + a6)

        Coefficient a6 is defined to satisfy H(298.15) = heat of formation at 298.15.

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        self.assert_all_in_bounds(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        return R * self.a6 + self.thermal_enthalpy(T)

    @unit_.manage_units([], D.energy_per_substance / D.temperature)
    def thermal_entropy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate entropy, S(T).

        Formula:
            S(T) = R (a1 ln(T) + a2 T + (a3/2) T^2 + (a4/3) T^3 + (a5/4) T^4 + a7)

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        self.assert_all_in_bounds(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(T, dtype=np.float64)  # noqa: N806
        return R * (
            self.a1 * np.log(T)
            + self.a2 * T
            + (self.a3 / 2) * T**2
            + (self.a4 / 3) * T**3
            + (self.a5 / 4) * T**4
        )

    @unit_.manage_units([], D.energy_per_substance)
    def thermal_enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,  # noqa: ARG002
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy change, dH(T).

        Formula:
            dH(T) = R (a1 T + (a2/2) T^2 + (a3/3) T^3 + (a4/4) T^4 + (a5/5) T^5)

        :param T: Temperature(s)
        :param units: Unit system
        :return: Function value(s)
        """
        self.assert_all_in_bounds(T)

        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        T = np.array(T, dtype=np.float64)  # noqa: N806
        return R * (
            self.a1 * T
            + (self.a2 / 2) * T**2
            + (self.a3 / 3) * T**3
            + (self.a4 / 4) * T**4
            + (self.a5 / 5) * T**5
        )

    @classmethod
    def fit(
        cls,
        T: ArrayLike,  # noqa: N803
        Cp: ArrayLike,  # noqa: N803
        S: ArrayLike,  # noqa: N803
        H: ArrayLike,  # noqa: N803
        T_mid: float = 1000,  # noqa: N803
        T_min: float | None = None,  # noqa: N803
        T_max: float | None = None,  # noqa: N803
    ) -> tuple["Nasa7Calculator", "Nasa7Calculator"]:
        """Fit data to Nasa-7 calculator coefficients.

        Step 1: Fit heat capacity data to get a1 - a5 for low and high
        temperature ranges.

        Step 2: Determine integration constants as follows.

            a6_low = (
                sum_low     (H(T) - Hth_fit_low(T))
                + sum_high  (H(T) - Hth_fit_high(T) - delta_H_low)
            ) / (N * R)

            delta_H_low = Hth_fit_low(T_mid) - Hth_fit_high(T_mid)

        Equations for a6_high, a7_low, a7_high are obtained by replacing "low"
        with "high" and "H" with "S".

        :param T: Temperatures
        :param Cp: Constant-pressure heat capacities
        :param S: Entropies
        :param H: Enthalpies
        :return: Fitted object
        """
        calc_low = cls.heat_capacity_fit(
            T=T, Cp=Cp, T_min=T_min, T_max=T_mid, T_ref=T_mid
        )
        calc_high = cls.heat_capacity_fit(
            T=T, Cp=Cp, T_min=T_mid, T_max=T_max, T_ref=T_mid
        )

        low = (T_min <= T) * (T_mid >= T)
        high = (T_mid <= T) * (T_max >= T)

        Hth_high = calc_high.thermal_enthalpy(T_mid)  # noqa: N806
        Hth_low = calc_low.thermal_enthalpy(T_mid)  # noqa: N806
        Sth_high = calc_high.thermal_entropy(T_mid)  # noqa: N806
        Sth_low = calc_low.thermal_entropy(T_mid)  # noqa: N806
        delta_H_low = Hth_low - Hth_high  # noqa: N806
        delta_H_high = Hth_high - Hth_low  # noqa: N806
        delta_S_low = Sth_low - Sth_high  # noqa: N806
        delta_S_high = Sth_high - Sth_low  # noqa: N806

        nT = np.sum(low) + np.sum(high)  # noqa: N806
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806

        # Enthalpy integration constants
        num_a6_low = np.sum(H[low] - calc_low.thermal_enthalpy(T[low])) + np.sum(
            H[high] - calc_high.thermal_enthalpy(T[high]) - delta_H_low
        )
        num_a6_high = np.sum(
            H[low] - calc_low.thermal_enthalpy(T[low]) - delta_H_high
        ) + np.sum(H[high] - calc_high.thermal_enthalpy(T[high]))

        # Entropy integration constants
        num_a7_low = np.sum(S[low] - calc_low.thermal_entropy(T[low])) + np.sum(
            S[high] - calc_high.thermal_entropy(T[high]) - delta_S_low
        )
        num_a7_high = np.sum(
            S[low] - calc_low.thermal_entropy(T[low]) - delta_S_high
        ) + np.sum(S[high] - calc_high.thermal_entropy(T[high]))

        a6_low = num_a6_low / (nT * R)
        a6_high = num_a6_high / (nT * R)
        a7_low = num_a7_low / (nT * R)
        a7_high = num_a7_high / (nT * R)
        return calc_low.model_copy(
            update={"a6": a6_low, "a7": a7_low}
        ), calc_high.model_copy(update={"a6": a6_high, "a7": a7_high})

    @classmethod
    def heat_capacity_fit(
        cls,
        T: ArrayLike,  # noqa: N803
        Cp: ArrayLike,  # noqa: N803
        T_ref: float = 1000,  # noqa: N803
        T_min: float | None = None,  # noqa: N803
        T_max: float | None = None,  # noqa: N803
    ) -> "Nasa7Calculator":
        """Fit data to Nasa-7 calculator coefficients.

        Parametrization:

            Cp/R = A0 + A1 T + A2 T^2 + A3 T^3 + A4 T^4
            H/R = A0 T + (A1/2) T^2 + (A2/3) T^3 + (A3/4) T^4 + (A4/5) T^5 + A5
            S/R = A0 ln(T) + A1 T + (A2/2) T^2 + (A3/3) T^3 + (A4/4) T^4 + A6

        1. Reduced linear equation for A0 - A4:

            [t1 - 1, t1^2 - 1, t1^3 - 1, t1^4 - 1] [a1]   [(Cp(T1) - Cp(T_ref))/R]
            [t2 - 1, t2^2 - 1, t2^3 - 1, t2^4 - 1] [a2] = [(Cp(T2) - Cp(T_ref))/R]
            [   ...,      ...,      ...,      ...] [a3]   [                ...]
            [tN - 1, tN^2 - 1, tN^3 - 1, tN^4 - 1] [a4]   [(Cp(TN) - Cp(T_ref))/R]

            t = T / T_ref

            A0 = Cp(T_ref) / R - (a1 + a2 + a3 + a4)
            A1 = a1 / T_ref
            A2 = a2 / T_ref^2
            A3 = a3 / T_ref^3
            A4 = a4 / T_ref^4

        :param T: Temperatures
        :param Cp: Constant-pressure heat capacities
        :param S: Entropies
        :param H: Enthalpies
        :return: Fitted object
        """
        T = np.array(T, dtype=np.float64)  # noqa: N806
        select = (T_min <= T) * (T_max >= T)

        Cp_refs = np.compress(np.isclose(T, T_ref), Cp)  # noqa: N806
        if not Cp_refs.size:
            msg = f"No temperature matches {T_ref = }: {T = }"
            raise ValueError(msg)

        (Cp_ref,) = Cp_refs  # noqa: N806

        t = np.compress(select, T / T_ref)
        _1 = np.ones_like(t)

        # Transformation matrix
        M = np.column_stack([t - _1, t**2 - _1, t**3 - _1, t**4 - _1])  # noqa: N806

        # Data vector
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        v = (np.compress(select, Cp) - Cp_ref) / R

        (a1, a2, a3, a4), *_ = np.linalg.lstsq(M, v)

        return cls(
            T_min=T_min,
            T_max=T_max,
            a1=Cp_ref / R - (a1 + a2 + a3 + a4),
            a2=a1 / T_ref,
            a3=a2 / T_ref**2,
            a4=a3 / T_ref**3,
            a5=a4 / T_ref**4,
            a6=0,
            a7=0,
        )
