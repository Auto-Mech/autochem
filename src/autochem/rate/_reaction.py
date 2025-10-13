"""Rates."""

import itertools
import math
from collections.abc import Mapping, Sequence
from typing import ClassVar

import altair as alt
import numpy as np
import pydantic

from ..unit_ import UnitsData
from ..util import chemkin, mess, plot
from ..util.type_ import Scalable, Scalers
from . import data
from .data import ArrheniusRateFit, PlogRateFit, Rate, Rate_, RateFit


class Reaction(Scalable):
    """Rate class."""

    reactants: list[str]
    products: list[str]
    reversible: bool = True
    rate: Rate_ = pydantic.Field(
        default_factory=lambda data: ArrheniusRateFit(
            A=1, b=0, E=0, order=len(data["reactants"])
        ),
    )

    # Private attributes
    _scalers: ClassVar[Scalers] = {"rate": (lambda c, x: c * x)}

    def __add__(self, other: "Reaction") -> "Reaction":
        res = self.model_copy()
        res.rate = self.rate + other.rate
        return res


# Constructors
def from_chemkin_string(
    rxn_str: str,
    *,
    units: UnitsData | None = None,
    strict: bool = True,
) -> Reaction:
    """Read rate from Chemkin string.

    :param rate_str: Chemkin string
    :param units: Units
    :param strict: Whether to fail if there are unused aux keys
    :return: Rate
    """
    # Parse string
    res = chemkin.parse_rate(rxn_str)

    # Extract rate constant
    rate_fit = data.from_chemkin_parse_results(res, units=units)

    # Check that all information was used
    if strict:
        assert not res.aux_numbers, f"Unused auxiliary values: {res.aux_numbers}"
        assert not res.aux_misc, f"Unused auxiliary values: {res.aux_misc}"
        assert not res.efficiencies, f"Unused efficiencies: {res.efficiencies}"

    # Instantiate object
    return Reaction(
        reactants=res.reactants,
        products=res.products,
        reversible=res.reversible,
        rate=rate_fit,
    )


def from_mess_channel_output(
    mess_chan_out: str, *, reversible: bool = True
) -> Reaction:
    """Extract rate data from MESS output.

    :param mess_chan_out: MESS output channel string
    :param order: Order
    :return: Rate data
    """
    res = mess.parse_output_channel(mess_chan_out)
    assert res.id1 is not None, res
    assert res.id2 is not None, res

    order = 1 if res.id1.startswith("W") else 2
    reactants = [res.id1]
    products = [res.id2]
    rate = data.from_mess_channel_output_parse_results(res, order=order)

    # Instantiate object
    return Reaction(
        reactants=reactants,
        products=products,
        reversible=reversible,
        rate=rate,
    )


# Transformations
def expand_lumped(
    rxn: Reaction, exp_dct: Mapping[str, Sequence[str]]
) -> list[Reaction]:
    """Expand a lumped reaction rates into its components.

    Assumes an even ratio among unlumped coefficients, in the absence of information.

        unlumped rate
        = lumped rate x
            nexp ^ stoich / multiset(nexp, stoich) for each lumped reactant
            1             / multiset(nexp, stoich) for each lumped product

    Here, nexp is the number of components in the lump and stoich is its stoichiometry
    in the reaction.

    There is no physical meaning to the individual rates, and some of the reactions may
    be unphysical. This only serves to reproduce the same net rate while distinguishing
    lump components.

    :param rxn: Reaction
    :param exp_dct: Mapping of lumped species to lump components
    :return: Component reactions
    """

    def _expand(name: str, *, rev: bool = True) -> tuple[float, list[dict[int, str]]]:
        """Determine reaction expansion and scale factor for one lumped species.

        Reaction expansion given as list of index -> name mappings representing
        different combinations of lump components.
        """
        # Get species expansion
        exp = exp_dct.get(name)
        if exp is None:
            return 1.0, [{}]

        # Get name combinations
        name_pool = rxn.products if rev else rxn.reactants
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

    rexps = [_expand(n, rev=False) for n in set(rxn.reactants) if n in exp_dct]
    pexps = [_expand(n, rev=True) for n in set(rxn.products) if n in exp_dct]

    rfactors, rexp_dcts = zip(*rexps, strict=True) if rexps else ((), ())
    pfactors, pexp_dcts = zip(*pexps, strict=True) if pexps else ((), ())

    # Scale rate by calculated factor
    rxn0 = rxn.model_copy()
    rxn0 *= math.prod(rfactors + pfactors)

    # Expand reactions
    rexp_combs = [
        {k: v for d in ds for k, v in d.items()} for ds in itertools.product(*rexp_dcts)
    ]
    pexp_combs = [
        {k: v for d in ds for k, v in d.items()} for ds in itertools.product(*pexp_dcts)
    ]
    rxns = []
    for rexp_comb, pexp_comb in itertools.product(rexp_combs, pexp_combs):
        rxn_ = rxn0.model_copy()
        rxn_.reactants = [
            rexp_comb.get(i) if i in rexp_comb else s
            for i, s in enumerate(rxn0.reactants)
        ]
        rxn_.products = [
            pexp_comb.get(i) if i in pexp_comb else s
            for i, s in enumerate(rxn0.products)
        ]
        rxns.append(rxn_)
    return rxns


# Conversions
def chemkin_equation(rxn: Reaction) -> str:
    """Get Chemkin equation string.

    :param rxn: Reaction
    :return: Chemkin equation string
    """
    assert isinstance(rxn.rate, RateFit), "Rate must be a RateFit"

    return chemkin.write_equation(
        reactants=rxn.reactants,
        products=rxn.products,
        reversible=rxn.reversible,
        third_body=rxn.rate.third_body,
        pressure_dependent=rxn.rate.is_pressure_dependent,
    )


def chemkin_string(rxn: Reaction, *, eq_width: int = 55, dup: bool = False) -> str:
    """Get Chemkin rate string.

    :param rxn: Reaction
    :param eq_width: Width for equation
    :param duplicate: Whether this is a duplicate reaction
    :return: Chemkin rate string
    """
    eq = chemkin_equation(rxn)
    rate_str = data.chemkin_string(rxn.rate, eq_width=eq_width)
    rxn_str = f"{eq:<{eq_width}} {rate_str}"
    return chemkin.write_with_dup(rxn_str, dup=dup)


# Fitting
def fit_high(rxn: Reaction) -> Reaction:
    """Fit high-pressure limit rate data to single Arrhenius.

    :param rxn: Reaction with rate data
    :return: Reaction with rate fit
    """
    rxn = rxn.model_copy()
    rate = rxn.rate
    assert isinstance(rate, Rate), rate
    rxn.rate = ArrheniusRateFit.fit(T=rate.T, k=rate.k_high, order=rate.order)
    return rxn


def fit_plog(rxn: Reaction, *, sanitize: bool = False) -> Reaction:
    """Fit rate data to Plog.

    :param rxn: Reaction with rate data
    :return: Reaction with rate fit
    """
    rxn = rxn.model_copy()
    rate = rxn.rate
    assert isinstance(rate, Rate), rate
    if sanitize:
        rate = rate.without_nan()
    rxn.rate = PlogRateFit.fit(
        T=rate.T,
        P=rate.P,
        k_data=rate.k_data,
        k_high=rate.k_high,
        order=rate.order,
    )
    return rxn


# Display
def display(  # noqa: PLR0913
    rxn: Reaction | Sequence[Reaction],
    *,
    T_range: tuple[float, float] = (400, 1250),  # noqa: N803
    P: float = 1,  # noqa: N803
    units: UnitsData | None = None,
    label: str | Sequence[str] | None = None,
    color: str | Sequence[str] | None = None,
    x_label: str = "1000/𝑇",  # noqa: RUF001
    y_label: str = "𝑘",  # noqa: RUF001
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
    rxns = [rxn] if isinstance(rxn, Reaction) else rxn
    labels = [label] if isinstance(label, str) else label
    colors = [color] if isinstance(color, str) else color

    rates = [r.rate for r in rxns]
    rate0, *rates_ = rates
    for rate_ in rates_:
        assert rate0.order == rate_.order, f"{rate0} !~ {rate_}"
    order = rate0.order

    nr = len(rates)
    labels = labels or ([f"k{i + 1}" for i in range(nr)] if nr > 1 else None)

    def make_chart(
        ixs: Sequence[int],
        rates: Sequence[data.BaseRate],
        labels: Sequence[str] | None,
        colors: Sequence[str] | None,
        mark: str,
    ) -> alt.Chart:
        rates_ = [rates[i] for i in ixs]
        labels_ = None if labels is None else [labels[i] for i in ixs]
        colors_ = None if colors is None else [colors[i] for i in ixs]
        (T, *Ts), ks = zip(  # noqa: N806
            *(r.plot_data(T_range=T_range, P=P, units=units) for r in rates_),
            strict=True,
        )
        for T_ in Ts:  # noqa: N806
            assert np.allclose(T, T_), f"{T} !~ {T_}"
        return plot.arrhenius(
            ks=ks,
            T=T,
            order=order,
            units=units,
            labels=labels_,
            colors=colors_,
            x_label=x_label,
            y_label=y_label,
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
