"""Rates."""

import itertools
import math
from collections.abc import Mapping, Sequence
from typing import ClassVar

import numpy
import pint
import pydantic
from numpy.typing import ArrayLike, NDArray

from ..unit_ import UnitsData
from ..util import chemkin
from ..util.type_ import Scalable, Scalers
from ._01const import (
    ArrheniusRateConstant,
    ParamRateConstant,
    RateConstant_,
    extract_rate_constant_from_chemkin_parse_results,
)
from ._01const import (
    chemkin_string as rate_constant_chemkin_string,
)


class Rate(Scalable):
    """Rate class."""

    reactants: list[str]
    products: list[str]
    reversible: bool = True
    rate_constant: RateConstant_ = pydantic.Field(
        default_factory=lambda data: ArrheniusRateConstant(
            A=1, b=0, E=0, order=len(data["reactants"])
        )
    )

    # Private attributes
    _scalers: ClassVar[Scalers] = {"rate_constant": (lambda c, x: c * x)}

    @property
    def unit(self) -> pint.Unit:
        """Unit."""
        return self.rate_constant.unit

    @property
    def third_body(self) -> str | None:
        """Third body."""
        if isinstance(self.rate_constant, ParamRateConstant):
            return self.rate_constant.third_body

        return None

    @property
    def is_pressure_dependent(self) -> bool:
        """Whether the rate is pressure dependent."""
        return isinstance(self.rate_constant, ArrheniusRateConstant)

    def __call__(
        self, t: ArrayLike, p: ArrayLike = 1, units: UnitsData | None = None
    ) -> NDArray[numpy.float64]:
        """Evaluate rate constant.

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
        return self.rate_constant(t=t, p=p, units=units)


# Constructors
def from_chemkin_string(
    rate_str: str, units: UnitsData | None = None, strict: bool = True
) -> Rate:
    """Read rate from Chemkin string.

    :param rate_str: Chemkin string
    :param units: Units
    :param strict: Whether to fail if there are unused aux keys
    :return: Rate
    """
    # Parse string
    res = chemkin.parse_rate(rate_str)

    # Extract rate constant
    rate_constant = extract_rate_constant_from_chemkin_parse_results(res, units=units)

    # Check that all information was used
    if strict:
        assert not res.aux_numbers, f"Unused auxiliary values: {res.aux_numbers}"
        assert not res.aux_misc, f"Unused auxiliary values: {res.aux_misc}"
        assert not res.efficiencies, f"Unused efficiencies: {res.efficiencies}"

    # Instantiate object
    return Rate(
        reactants=res.reactants,
        products=res.products,
        reversible=res.reversible,
        rate_constant=rate_constant,
    )


# Transformations
def expand_lumped(rate: Rate, exp_dct: Mapping[str, Sequence[str]]) -> list[Rate]:
    """Expand a lumped reaction rates into its components.

    Assumes an even ratio among unlumped coefficients, in the absence of information.

        unlumped rate coefficient
        = lumped rate coefficient x
            nexp ^ stoich / multiset(nexp, stoich) for each lumped reactant
            1             / multiset(nexp, stoich) for each lumped product

    Here, nexp is the number of components in the lump and stoich is its stoichiometry
    in the reaction.

    There is no physical meaning to the individual rates, and some of the reactions may
    be unphysical. This only serves to reproduce the same net rate while distinguishing
    lump components.

    :param rate: Rate
    :param exp_dct: Mapping of lumped species to lump components
    :return: Component rates
    """

    def _expand(name: str, rev: bool) -> tuple[float, list[dict[int, str]]]:
        """Determine reaction expansion and scale factor for one lumped species.

        Reaction expansion given as list of index -> name mappings representing
        different combinations of lump components.
        """
        # Get species expansion
        exp = exp_dct.get(name)
        if exp is None:
            return 1.0, [{}]

        # Get name combinations
        name_pool = rate.products if rev else rate.reactants
        stoich = name_pool.count(name)
        name_combs = list(itertools.combinations_with_replacement(exp, stoich))
        # Determine factor
        factor = 1.0 if rev else len(exp) ** stoich
        factor /= len(name_combs)
        # Determine reaction expansion dictionaries
        name_idxs = [i for i, n in enumerate(name_pool) if n == name]
        exp_dcts = [
            dict(zip(name_idxs, name_comb, strict=True)) for name_comb in name_combs
        ]
        return factor, exp_dcts

    rexps = [_expand(n, rev=False) for n in set(rate.reactants) if n in exp_dct]
    pexps = [_expand(n, rev=True) for n in set(rate.products) if n in exp_dct]

    rfactors, rexp_dcts = zip(*rexps, strict=True) if rexps else ((), ())
    pfactors, pexp_dcts = zip(*pexps, strict=True) if pexps else ((), ())

    # Scale rate by calculated factor
    rate0 = rate.model_copy()
    rate0 *= math.prod(rfactors + pfactors)

    # Expand reactions
    rexp_combs = [
        {k: v for d in ds for k, v in d.items()} for ds in itertools.product(*rexp_dcts)
    ]
    pexp_combs = [
        {k: v for d in ds for k, v in d.items()} for ds in itertools.product(*pexp_dcts)
    ]
    rates = []
    for rexp_comb, pexp_comb in itertools.product(rexp_combs, pexp_combs):
        rate_ = rate0.model_copy()
        rate_.reactants = [
            rexp_comb.get(i) if i in rexp_comb else s
            for i, s in enumerate(rate0.reactants)
        ]
        rate_.products = [
            pexp_comb.get(i) if i in pexp_comb else s
            for i, s in enumerate(rate0.products)
        ]
        rates.append(rate_)
    return rates


# Conversions
def chemkin_equation(rate: Rate) -> str:
    """Get Chemkin equation string.

    :param rate: Rate
    :return: Chemkin equation string
    """
    return chemkin.write_equation(
        reactants=rate.reactants,
        products=rate.products,
        reversible=rate.reversible,
        third_body=rate.third_body,
        pressure_dependent=rate.is_pressure_dependent,
    )


def chemkin_string(rate: Rate, eq_width: int = 55, dup: bool = False) -> str:
    """Get Chemkin rate string.

    :param rate: Rate
    :param eq_width: Width for equation
    :param duplicate: Whether this is a duplicate reaction
    :return: Chemkin rate string
    """
    eq = chemkin_equation(rate)
    rate_str = rate_constant_chemkin_string(rate.rate_constant, eq_width=eq_width)
    reac_str = f"{eq:<{eq_width}} {rate_str}"
    return chemkin.write_with_dup(reac_str, dup=dup)
