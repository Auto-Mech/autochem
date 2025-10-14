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
    def delta_enthalpy(
        self,
        T: ArrayLike,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[np.float64]:
        """Evaluate enthalpy change, dH(T).

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
        T = np.array(T, dtype=np.float64)  # noqa: N806
        return R * (
            self.a1 * np.log(T)
            + self.a2 * T
            + (self.a3 / 2) * T**2
            + (self.a4 / 3) * T**3
            + (self.a5 / 4) * T**4
            + self.a7
        )

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
        return R * self.a6 + self.delta_enthalpy(T)

    @unit_.manage_units([], D.energy_per_substance)
    def delta_enthalpy(
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
    ) -> "Nasa7Calculator":
        """Fit data to Nasa-7 calculator coefficients.

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
        :return: Fitted object
        """
        T = np.array(T, dtype=np.float64)  # noqa: N806
        _0 = np.zeros_like(T)
        _1 = np.ones_like(T)

        # Transformation matrix
        M_Cp = np.column_stack([_1, T, T**2, T**3, T**4, _0, _0])  # noqa: N806
        M_S = np.column_stack(  # noqa: N806
            [np.log(T), T, T**2 / 2, T**3 / 3, T**4 / 4, _0, _1],
        )
        M_H = np.column_stack([T, T**2 / 2, T**3 / 3, T**4 / 4, T**5 / 5, _1, _0])  # noqa: N806
        M = np.vstack((M_Cp, M_S, M_H))  # noqa: N806

        # Data vector
        R = unit_.const.value(C.gas, UNITS)  # noqa: N806
        v = np.concatenate((Cp / R, S / R, H / R))

        T_min, T_max = np.min(T), np.max(T)  # noqa: N806
        (a1, a2, a3, a4, a5, a6, a7), *_ = np.linalg.lstsq(M, v, rcond=1e-24)
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
