"""Blending function models."""

import abc
from typing import Annotated, TypeVar

import numpy
import pydantic
from numpy.typing import ArrayLike
from pydantic_core import core_schema

from ..util import chemkin
from ..util.type_ import Frozen, Scalable, SubclassTyped

F = TypeVar("F", bound="BlendingFunction")


class BlendingFunction(Frozen, Scalable, SubclassTyped, abc.ABC):
    """Abstract base class for blending functions."""

    @abc.abstractmethod
    def __call__(self, t: ArrayLike, p_r: ArrayLike) -> numpy.ndarray:
        """Evaluate blending function, f(t, p_r).

        :param t: Temperature(s)
        :param p_r: Reduced pressure(s)
        :return: Blending function value(s)
        """
        pass


class LindemannBlendingFunction(BlendingFunction):

    # Private attributes
    _type: str = "lindemann"

    def __call__(self, t: ArrayLike, p_r: ArrayLike) -> numpy.ndarray:
        """Evaluate blending function, f(t, p_r)."""
        return numpy.array(1.0)


class TroeBlendingFunction(BlendingFunction):
    a: float
    t3: float
    t1: float
    t2: float | None = None

    # Private attributes
    _type: str = "troe"

    def __call__(self, t: ArrayLike, p_r: ArrayLike) -> numpy.ndarray:
        """Evaluate blending function, f(t, p_r)."""
        log_f = numpy.log10(self.f_cent(t)) / (1 + self.f1(t, p_r) ** 2)
        return numpy.power(10, log_f)

    def f_cent(self, t: ArrayLike) -> numpy.ndarray:
        """Evaluate center broadening factor."""
        a, t3, t1, t2 = (self.a, self.t3, self.t1, self.t2)
        f_cent = (1 - a) * numpy.exp(-numpy.divide(t, t3))
        f_cent += a * numpy.exp(-numpy.divide(t, t1))
        f_cent += 0.0 if t2 is None else numpy.exp(-numpy.divide(t2, t))
        return f_cent

    def n(self, t: ArrayLike) -> numpy.ndarray:
        """Evaluate N."""
        return 0.75 - 1.27 * numpy.log10(self.f_cent(t))

    def c(self, t: ArrayLike) -> numpy.ndarray:
        """Evaluate C."""
        return -0.4 - 0.67 * numpy.log10(self.f_cent(t))

    def f1(self, t: ArrayLike, p_r: ArrayLike) -> numpy.ndarray:
        """Evaluate f1."""
        n = self.n(t)
        c = self.c(t)
        log_p_r_plus_c = numpy.log10(p_r) + c
        return log_p_r_plus_c / (n - 0.14 * log_p_r_plus_c)


class SriBlendingFunction(BlendingFunction):
    a: float
    b: float
    c: float
    d: float = 1.0
    e: float = 0.0

    # Private attributes
    _type: str = "sri"

    def __call__(self, t: ArrayLike, p_r: ArrayLike) -> numpy.ndarray:
        """Evaluate blending function, f(t, p_r)."""
        a, b, c, d, e = (self.a, self.b, self.c, self.d, self.e)
        return (
            d
            * (a * numpy.exp(-numpy.divide(b, t)) + numpy.exp(-numpy.divide(t, c)))
            ** (1 / (1 + numpy.log10(p_r) ** 2))
            * numpy.power(t, e)
        )


BlendingFunction_ = Annotated[
    pydantic.SkipValidation[BlendingFunction],
    pydantic.BeforeValidator(lambda x: BlendingFunction.model_validate(x)),
    pydantic.PlainSerializer(lambda x: BlendingFunction.model_validate(x).model_dump()),
    pydantic.GetPydanticSchema(
        lambda _, handler: core_schema.with_default_schema(handler(pydantic.BaseModel))
    ),
]


# Properties
def chemkin_aux_lines(func: BlendingFunction, head_width: int = 55) -> list[str]:
    """Determine auxiliary Chemkin lines for blending function.

    :param func: Blending function
    :return: Auxiliary lines
    """
    match func:
        case TroeBlendingFunction():
            troe_params = [func.a, func.t3, func.t1]
            troe_params += [] if func.t2 is None else [func.t2]
            return [chemkin.write_aux("TROE", troe_params, head_width=head_width)]
        case SriBlendingFunction():
            sri_params = [func.a, func.b, func.c, func.d, func.e]
            return [chemkin.write_aux("SRI", sri_params, head_width=head_width)]
        case _:
            assert isinstance(func, LindemannBlendingFunction)
            return []


# Parse helpers
def extract_blending_function_from_chemkin_parse_results(
    res: chemkin.ChemkinRateParseResults,
) -> BlendingFunction:
    """Extract blending function from Chemkin parse results.

    Chemkin parse results are modified in-place

    :param res: Chemkin parse results
    :return: Blending function
    """
    if "SRI" in res.aux_numbers:
        numbers = res.aux_numbers.pop("SRI")
        keys = ["a", "b", "c", "d", "e"]
        data = dict(zip(keys, numbers, strict=False))
        return SriBlendingFunction.model_validate(data)

    if "TROE" in res.aux_numbers:
        numbers = res.aux_numbers.pop("TROE")
        keys = ["a", "t3", "t1", "t2"]
        data = dict(zip(keys, numbers, strict=False))
        return TroeBlendingFunction.model_validate(data)

    return LindemannBlendingFunction()
