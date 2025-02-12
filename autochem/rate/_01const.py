"""Rate constant models."""

import abc
from typing import Annotated, ClassVar

import more_itertools as mit
import numpy
import pint
import pydantic
import xarray
from numpy.polynomial import chebyshev
from numpy.typing import ArrayLike, NDArray
from pydantic_core import core_schema

from .. import unit_
from ..unit_ import UNITS, Dimension, UnitManager, Units, UnitsData
from ..util import chemkin
from ..util.type_ import Frozen, NDArray_, Scalable, Scalers, SubclassTyped
from ._00func import (
    BlendingFunction_,
    extract_blending_function_from_chemkin_parse_results,
)
from ._00func import (
    chemkin_aux_lines as blending_function_chemkin_aux_lines,
)


class RateConstant(UnitManager, Frozen, Scalable, SubclassTyped, abc.ABC):
    """Abstract base class for rate constants."""

    order: int = 1

    @property
    def unit(self) -> pint.Unit:
        return UNITS.rate_constant(self.order)

    def __init__(self, units: UnitsData | None = None, **kwargs):
        super().__init__(units=units, **kwargs)

    @abc.abstractmethod
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant, k(t, p, x).

        Uses:
        - If temperature and pressure are both numbers, the rate constant will be
          returned as a number, k(t, p).
        - If either temperature or pressure are both lists, the rate constant will be
          returned as a 1D array, [k(t1, p), k(t2, p), ...] or [k(t, p1), k(t, p2), ...]
        - If temperature and pressure are lists, the rate constant will be returned as a
          2D array, [[k(t1, p1), k(t1, p2), ...], [k(t2, p1), k(t2, p2), ...]]

        :param t: Temperature(s)
        :param p: Pressure(s)
        :param units: Input / desired output units
        :return: Value(s)
        """
        pass

    def process_input(
        self, t: ArrayLike, p: ArrayLike
    ) -> tuple[NDArray[numpy.float64], NDArray[numpy.float64]]:
        """Normalize rate constant input.

        :param t: Temperature(s)
        :param p: Pressure(s)
        :return: Temperature(s) and pressure(s)
        """
        # Convert to numpy arrays
        t = numpy.array(numpy.squeeze(t), dtype=numpy.float64)
        p = numpy.array(numpy.squeeze(p), dtype=numpy.float64)

        # Handle dimensions if 2D
        is_2d = numpy.ndim(t) > 0 and numpy.ndim(p) > 0
        t = t[:, numpy.newaxis] if is_2d else t
        p = p[numpy.newaxis, :] if is_2d else p

        return t, p


class RawRateConstant(RateConstant):
    ts: list[float]
    ps: list[float]
    k_array: NDArray_

    # Private attributes
    _t_key = "t"
    _p_key = "p"
    _type: str = "raw"
    _scalers: ClassVar[Scalers] = {"k_array": numpy.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "ts": Dimension.temperature,
        "ps": Dimension.pressure,
        "k_array": Dimension.rate_constant,
    }

    @property
    def ktp(self):
        return xarray.DataArray(
            data=self.k_array, coords={self._t_key: self.ts, self._p_key: self.ps}
        )

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant."""
        t, p = self.process_input(t, p)
        ktp: NDArray[numpy.float64] = self.ktp.sel(t=t, p=p, method="ffill").data
        return ktp


class ParamRateConstant(RateConstant):
    """Abstract base class for parametrized rate constants."""

    efficiencies: dict[str, float] = pydantic.Field(default_factory=dict)

    @property
    def third_body(self) -> str | None:
        """Get third body."""
        eff = self.efficiencies

        if eff.get("M") == 1:
            return "M"

        return next((c for c, e in eff.items() if e == 1.0), None)


class ArrheniusRateConstant(ParamRateConstant):
    A: float = 1.0
    b: float = 0.0
    E: float = 0.0

    # Private attributes
    _type: str = "arrhenius"
    _scalers: ClassVar[Scalers] = {"A": numpy.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "A": Dimension.rate_constant,
        "E": Dimension.energy_per_substance,
    }

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant."""
        t, p = self.process_input(t, p)
        r = unit_.system.gas_constant_value(UNITS)
        ktp = self.A * (t**self.b) * numpy.exp(-self.E / (r * t))
        return ktp


class BlendedRateConstant(ParamRateConstant, abc.ABC):  # type: ignore[misc]
    A_high: float
    b_high: float
    E_high: float
    A_low: float
    b_low: float
    E_low: float
    function: BlendingFunction_

    # Private attributes
    _scalers: ClassVar[Scalers] = {"A_high": numpy.multiply, "A_low": numpy.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "A_high": Dimension.rate_constant,
        "E_high": Dimension.energy_per_substance,
        "A_low": Dimension.rate_constant,
        "E_low": Dimension.energy_per_substance,
    }

    @property
    def arrhenius_functions(
        self,
    ) -> tuple[ArrheniusRateConstant, ArrheniusRateConstant]:
        k_low = ArrheniusRateConstant(
            A=self.A_high, b=self.b_high, E=self.E_high, order=self.order
        )
        k_high = ArrheniusRateConstant(
            A=self.A_high, b=self.b_high, E=self.E_high, order=self.order
        )
        return k_low, k_high

    def effective_concentration(
        self, t: NDArray[numpy.float64], p: NDArray[numpy.float64]
    ) -> NDArray[numpy.float64]:
        """Get effective concentration(s) from temperature(s) and pressure(s).

        effective [M] = P / R T  (ideal gas law)

        :param t: Temperature(s)
        :param p: Pressure(s)
        :return: Effective concentration(s)
        """
        # Evaluate, using pint to handle units
        r_ = pint.Quantity("molar_gas_constant")
        t_ = pint.Quantity(t, UNITS.temperature)
        p_ = pint.Quantity(p, UNITS.pressure)
        m_ = p_ / (r_ * t_)

        # Return value in concentration units
        return m_.m_as(UNITS.concentration)

    def effective_reduced_pressure(
        self, t: NDArray[numpy.float64], p: NDArray[numpy.float64]
    ) -> NDArray[numpy.float64]:
        """Get effective concentration(s) from temperature(s) and pressure(s).

        effective P_r = k_low [M] / k_high  (ideal gas law)

        :param t: Temperature(s)
        :param p: Pressure(s)
        :return: Effective reduced pressure(s)
        """
        m_eff = self.effective_concentration(t, p)
        k_low, k_high = self.arrhenius_functions
        return k_low(t) * m_eff / k_high(t)


class FalloffRateConstant(BlendedRateConstant):

    # Private attributes
    _type: str = "falloff"

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant."""
        t, p = self.process_input(t, p)
        _, k_high = self.arrhenius_functions
        p_r = self.effective_reduced_pressure(t, p)
        ktp = k_high(t) * p_r / (1 + p_r) * self.function(t, p_r)
        return ktp


class ActivatedRateConstant(BlendedRateConstant):

    # Private attributes
    _type: str = "activated"

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant."""
        t, p = self.process_input(t, p)
        k_low, _ = self.arrhenius_functions
        p_r = self.effective_reduced_pressure(t, p)
        ktp = k_low(t) / (1 + p_r) * self.function(t, p_r)
        return ktp


class PlogRateConstant(ParamRateConstant):
    As: list[float]
    bs: list[float]
    Es: list[float]
    ps: list[float]

    # Private attributes
    _type: str = "plog"
    _scalers: ClassVar[Scalers] = {"As": numpy.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "As": Dimension.rate_constant,
        "Es": Dimension.energy_per_substance,
        "ps": Dimension.pressure,
    }

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant for a single pressure."""
        t, p = self.process_input(t, p)
        p0 = self.nearest_pressure(p, which=0)
        p1 = self.nearest_pressure(p, which=1)
        kt0 = self.nearest_arrhenius_values(p, t, which=0)
        kt1 = self.nearest_arrhenius_values(p, t, which=1)

        # Evaluate intermediate pressures
        log_p, log_p0, log_p1 = map(numpy.log, (p, p0, p1))
        ktp = kt0 + (kt1 - kt0) * (log_p - log_p0) / (log_p1 - log_p0)

        # Evaluate on-boundary pressures (needed to fill in last pressure value)
        if numpy.ndim(p) > 0:
            ktp[..., numpy.equal(p, p0)] = kt0[..., numpy.equal(p, p0)]
        elif p == p0:
            ktp = kt0

        return ktp

    @property
    def arrhenius_functions(self) -> list[ArrheniusRateConstant]:
        return [
            ArrheniusRateConstant(A=A, b=b, E=E, order=self.order)
            for A, b, E in zip(self.As, self.bs, self.Es, strict=True)
        ]

    def arrhenius_values(self, t: ArrayLike) -> NDArray[numpy.float64]:
        """Functions."""
        return numpy.array(
            [k(t) for k in self.arrhenius_functions], dtype=numpy.float64
        ).T

    @property
    def pressures(self) -> NDArray[numpy.float64]:
        """Pressures."""
        return numpy.array(self.ps, dtype=numpy.float64)

    @property
    def pressure_indices(self) -> list[int]:
        """Pressure indices."""
        return list(range(self.pressures.shape[0]))

    def nearest_arrhenius_values(
        self, p: ArrayLike, t: ArrayLike, which: int = 0
    ) -> NDArray[numpy.float64]:
        """Get nearest lower or higher pressure.

        :param p: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest function(s)
        """
        ix = self.nearest_index(p, which=which)
        kts = self.arrhenius_values(t)
        return numpy.where(
            numpy.isin(ix, self.pressure_indices),
            numpy.take(kts, ix, mode="clip", axis=-1),
            numpy.nan,
        )

    def nearest_pressure(self, p: ArrayLike, which: int = 0) -> NDArray[numpy.float64]:
        """Get nearest lower or higher pressure.

        :param p: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest pressure(s)
        """
        ix = self.nearest_index(p, which=which)
        ps = self.pressures
        return numpy.where(
            numpy.isin(ix, self.pressure_indices),
            numpy.take(ps, ix, mode="clip"),
            numpy.nan,
        )

    def nearest_index(self, p: ArrayLike, which: int = 0) -> NDArray[numpy.int_]:
        """Get nearest lower or higher index.

        :param p: Pressure(s)
        :param which: 0=lower, 1=higher
        :return: Nearest index (indices)
        """
        val = numpy.searchsorted(self.ps, p, side="right") - 1 + which
        return val


class ChebRateConstant(ParamRateConstant):
    coeffs: NDArray_
    t_range: tuple[float, float]
    p_range: tuple[float, float]

    # Private attributes
    _type: str = "cheb"
    _scalers: ClassVar[Scalers] = {"coeffs": numpy.multiply}
    _dimensions: ClassVar[dict[str, Dimension]] = {
        "coeffs": Dimension.rate_constant,
        "t_range": Dimension.temperature,
        "p_range": Dimension.pressure,
    }

    @unit_.manage_units(
        [Dimension.temperature, Dimension.pressure], Dimension.rate_constant
    )
    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant for a single pressure."""
        t, p = self.process_input(t, p)
        t0, t1 = self.t_range
        p0, p1 = self.p_range

        inv_ = numpy.reciprocal
        log_ = numpy.log10

        t_ = (2 * inv_(t) - inv_(t0) - inv_(t1)) / (inv_(t1) - inv_(t0))
        p_ = (2 * log_(p) - log_(p0) - log_(p1)) / (log_(p1) - log_(p0))

        # AVC: I don't understand why I need to transpose the coefficient matrix here
        ktp = chebyshev.chebgrid2d(t_, p_, self.coeffs.T)
        return ktp


RateConstant_ = Annotated[
    pydantic.SkipValidation[RateConstant],
    pydantic.BeforeValidator(lambda x: RateConstant.model_validate(x)),
    pydantic.PlainSerializer(lambda x: RateConstant.model_validate(x).model_dump()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(pydantic.BaseModel))
    ),
]


# Conversions
def chemkin_string(rate_const: ParamRateConstant, eq_width: int = 0) -> str:
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
        case ArrheniusRateConstant():
            head_params = [rate_const.A, rate_const.b, rate_const.E]
            head_line = chemkin.write_numbers(head_params)
            aux_lines = []
        case ActivatedRateConstant():
            high_params = [rate_const.A_high, rate_const.b_high, rate_const.E_high]
            low_params = [rate_const.A_low, rate_const.b_low, rate_const.E_low]
            head_line = chemkin.write_numbers(low_params)
            aux_lines = [chemkin.write_aux("HIGH", high_params, head_width=head_width)]
        case FalloffRateConstant():
            high_params = [rate_const.A_high, rate_const.b_high, rate_const.E_high]
            low_params = [rate_const.A_low, rate_const.b_low, rate_const.E_low]
            head_line = chemkin.write_numbers(high_params)
            aux_lines = [chemkin.write_aux("LOW", low_params, head_width=head_width)]
        case PlogRateConstant():
            plog_params = [rate_const.ps, rate_const.As, rate_const.bs, rate_const.Es]
            aux_lines = [
                chemkin.write_aux("PLOG", row) for row in zip(*plog_params, strict=True)
            ]
        case ChebRateConstant():
            shape = numpy.shape(rate_const.coeffs)
            aux_lines = [
                chemkin.write_aux("TCHEB", rate_const.t_range, head_width=head_width),
                chemkin.write_aux("PCHEB", rate_const.p_range, head_width=head_width),
                chemkin.write_aux("CHEB", shape, head_width=head_width, as_int=True),
                *(
                    chemkin.write_aux("CHEB", row, head_width=head_width)
                    for row in rate_const.coeffs.tolist()
                ),
            ]
        case _:
            raise ValueError(
                f"Rate constant has unknown type {type(rate_const)}:\n{rate_const}"
            )

    if isinstance(rate_const, BlendedRateConstant):
        aux_lines.extend(
            blending_function_chemkin_aux_lines(
                rate_const.function, head_width=head_width
            )
        )

    eff_line = chemkin.write_efficiencies(
        rate_const.efficiencies, third_body=rate_const.third_body
    )
    if eff_line is not None:
        aux_lines.append(eff_line)

    lines = [head_line, *aux_lines]
    return "\n".join(lines)


# Parse helpers
def extract_rate_constant_from_chemkin_parse_results(
    res: chemkin.ChemkinRateParseResults, units: UnitsData | None = None
) -> ParamRateConstant:
    """Extract blending function from Chemkin parse results.

    Chemkin parse results are modified in-place

    :param res: Chemkin parse results
    :return: Blending function
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
        coeffs = numpy.reshape(cheb[2:], shape)
        # Read ranges
        t_range = res.aux_numbers.pop("TCHEB")
        p_range = res.aux_numbers.pop("PCHEB")
        return ChebRateConstant(
            coeffs=coeffs,
            t_range=t_range,
            p_range=p_range,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "PLOG" in res.aux_numbers:
        ps, As, bs, Es = zip(
            *mit.chunked(res.aux_numbers.pop("PLOG"), 4, strict=True), strict=True
        )
        return PlogRateConstant(
            As=As,
            bs=bs,
            Es=Es,
            ps=ps,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "LOW" in res.aux_numbers:
        A_high, b_high, E_high = res.arrhenius
        A_low, b_low, E_low = res.aux_numbers.pop("LOW")
        function = extract_blending_function_from_chemkin_parse_results(res)
        return FalloffRateConstant(
            A_low=A_low,
            b_low=b_low,
            E_low=E_low,
            A_high=A_high,
            b_high=b_high,
            E_high=E_high,
            function=function,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    if "HIGH" in res.aux_numbers:
        A_low, b_low, E_low = res.arrhenius
        A_high, b_high, E_high = res.aux_numbers.pop("HIGH")
        function = extract_blending_function_from_chemkin_parse_results(res)
        return ActivatedRateConstant(
            A_low=A_low,
            b_low=b_low,
            E_low=E_low,
            A_high=A_high,
            b_high=b_high,
            E_high=E_high,
            function=function,
            efficiencies=efficiencies,
            order=order,
            units=units,
        )

    A, b, E = res.arrhenius
    return ArrheniusRateConstant(
        A=A, b=b, E=E, efficiencies=efficiencies, order=order, units=units
    )
