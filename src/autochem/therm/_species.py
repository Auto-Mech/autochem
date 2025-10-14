"""Thermodynamic data."""

import datetime
from collections.abc import Sequence
from typing import ClassVar, Literal

import altair as alt
import numpy as np
import pyparsing as pp

from ..unit_ import UnitsData
from ..util import FormulaData, chemkin, form, pac99
from ..util.type_ import Scalable, Scalers
from . import data
from .data import Nasa7ThermFit, Therm, Therm_


# TODO(avcopan): Include subtype: `Species[Therm]` / `Species[ThermFit]`  # noqa: FIX002
# https://github.com/Auto-Mech/autochem/issues/671
class Species(Scalable):
    """A species with thermodynamic data."""

    name: str
    therm: Therm_

    # Private attributes
    _scalers: ClassVar[Scalers] = {"therm": (lambda c, x: c * x)}


# Properties
def temperature_minimum(spc: Species) -> float:
    """Get the minimum temperature of the species.

    :param spc: Species
    :return: Minimum temperature
    """
    if isinstance(spc.therm, Nasa7ThermFit):
        return spc.therm.T_min
    msg = f"Species {spc.name} has no minimum temperature for {type(spc.therm)}"
    raise NotImplementedError(msg)


def temperature_middle(spc: Species) -> float:
    """Get the middle temperature of the species.

    :param spc: Species
    :return: Middle temperature
    """
    if isinstance(spc.therm, Nasa7ThermFit):
        return spc.therm.T_mid
    msg = f"Species {spc.name} has no middle temperature for {type(spc.therm)}"
    raise NotImplementedError(msg)


def temperature_maximum(spc: Species) -> float:
    """Get the maximum temperature of the species.

    :param spc: Species
    :return: Maximum temperature
    """
    if isinstance(spc.therm, Nasa7ThermFit):
        return spc.therm.T_max
    msg = f"Species {spc.name} has no maximum temperature for {type(spc.therm)}"
    raise NotImplementedError(msg)


# Conversions
def from_chemkin_string(
    spc_str: str,
    T_mid: float | None = None,  # noqa: N803
) -> Species:
    """Read species thermo from Chemkin string.

    :param spc_str: Chemkin species therm string
    :param T_min: Default minimum temperature
    :param T_mid: Default middle temperature
    :param T_max: Default maximum temperature
    :return: Species thermo
    """
    # Parse string
    res = chemkin.parse_thermo(spc_str)

    # Extract thermo data
    therm_fit = data.from_chemkin_parse_results(res, T_mid=T_mid)

    return Species(name=res.name, therm=therm_fit)


def from_messpf_output_string(  # noqa: PLR0913
    pf_str: str,
    formula: FormulaData,
    name: str | None = None,
    charge: int = 0,
    Hf: float | None = None,  # noqa: N803
    Tf: float = 0,  # noqa: N803
    units: UnitsData | None = None,
) -> Species:
    """Read species thermo from MESS-PF output string.

    :param pf_str: MESS-PF output string
    :param formula: Species formula, as string or dict
    :param name: Species name
    :param charge: Molecular charge
    :return: Species thermo
    """
    if name is None:
        prefix_expr = pp.Literal("T,") + pp.Literal("K")
        name_expr = pp.Word(pp.printables)
        expr = pp.SkipTo(prefix_expr) + prefix_expr + name_expr("name")
        name = expr.parse_string(pf_str).get("name")

    therm = data.from_messpf_output_string(
        pf_str,
        formula=formula,
        charge=charge,
        Hf=Hf,
        Tf=Tf,
        units=units,
    )
    return Species(name=name, therm=therm)


def from_pac99_output_string(spc_str: str) -> Species:
    """Read species thermo from PAC99 output string.

    :param spc_str: PAC99 .c97 output string
    :return: Species thermo
    """
    # Parse string
    res = pac99.parse_thermo(spc_str)

    # Extract thermo data
    therm_fit = data.from_pac99_output_parse_results(res)

    return Species(name=res.name, therm=therm_fit)


def chemkin_string(spc: Species) -> str:
    """Get Chemkin thermo string.

    :param therm: Species thermo
    :return: Chemkin thermo string
    """
    therm = spc.therm
    T_min = T_max = T_mid = None  # noqa: N806
    match therm:
        case Nasa7ThermFit():
            T_min = therm.T_min  # noqa: N806
            T_max = therm.T_max  # noqa: N806
            T_mid = therm.T_mid  # noqa: N806
            coeffs = therm.coeffs_high + therm.coeffs_low
        case _:
            msg = f"Thermodynamic data type {type(therm)} not implemented"
            raise NotImplementedError(msg)

    line1 = chemkin.write_therm_entry_header(
        name=spc.name,
        form_dct=therm.formula,
        T_min=T_min,
        T_max=T_max,
        T_mid=T_mid,
        charge=therm.charge,
    )
    lines = chemkin.write_therm_entry_coefficient_lines(coeffs)

    lines = [line1, *lines]
    return "\n".join(f"{L: <78}{i + 1:>2d}" for i, L in enumerate(lines))


def pac99_input_string(
    spc: Species,
    T_min: float = 200,  # noqa: N803
    T_mid: float = 1000,  # noqa: N803
    T_max: float = 3000,  # noqa: N803
) -> str:
    """Generate a PAC99 input string for fitting to a NASA-7 polynomial.

    :param spc: Species thermo with thermo data (not a fit)
    :return: PAC99 input string
    """
    assert isinstance(spc.therm, Therm)
    fml_str = form.string(spc.therm.formula, ones=True)
    Hf298 = spc.therm.enthalpy_of_formation(units={"energy": "J"})  # noqa: N806
    Ts = spc.therm.T  # noqa: N806
    # Enthalpy units are set to kJ by "KJOULE" keyword below
    dH298 = spc.therm.delta_enthalpy(T=298, method="nearest", units={"energy": "kJ"})  # noqa: N806
    dHs = spc.therm.delta_enthalpy_data(units={"energy": "kJ"})  # noqa: N806
    # Entropy and heat capacity are always in Joules
    # (see https://ntrs.nasa.gov/citations/19930003779, p. 33)
    Ss = spc.therm.entropy_data(P=1, units={"pressure": "bar", "energy": "J"})  # noqa: N806
    Cs = spc.therm.heat_capacity_data(const="P", units={"energy": "J"})  # noqa: N806
    date = format(datetime.date.today(), r"%Y%m%d")  # noqa: DTZ011
    return "\n".join(
        [
            pac99_input_line("NAME", spc.name),
            pac99_input_line(fml_str, None, None, "HF298", Hf298, "JOULES", decimals=0),
            pac99_input_line("DATE", date),
            pac99_input_line("REFN", "ME"),
            pac99_input_line(
                "LSTS",
                "OLD",
                None,
                "T",
                T_min,
                "T",
                T_min,
                "T",
                T_mid,
                decimals=0,
            ),
            pac99_input_line("LSTS", "T", T_max, decimals=0),
            pac99_input_line("OUTP", "LSQS"),
            pac99_input_line("METH", "READIN", None, "KJOULE", None, "BAR"),
            pac99_input_line(None, "T", 0, "CP", 0, "S", 0, "H-H2", -dH298),
            *(
                pac99_input_line(None, "T", T, "CP", C, "S", S, "H-H0", dH)
                for T, C, S, dH in zip(Ts, Cs, Ss, dHs, strict=True)
            ),
            pac99_input_line("FINISH"),
        ],
    )


def fit(
    spc: Species,
    T_min: float | None = None,  # noqa: N803
    T_mid: float = 1000,  # noqa: N803
    T_max: float | None = None,  # noqa: N803
    type_: Literal["nasa7"] = "nasa7",  # noqa: ARG001
) -> Species:
    """Fit data to therm fit object.

    :param spc: Species thermo with thermo data (not a fit)
    :return: Species thermo with fitted data
    """
    therm: Therm = spc.therm
    T = therm.temperature_data()  # noqa: N806
    Cp = therm.heat_capacity_data(const="P")  # noqa: N806
    S = therm.entropy_data(P=1, units={"pressure": "bar"})  # noqa: N806
    H = therm.enthalpy_data()  # noqa: N806

    T_min = T_min or np.min(T)  # noqa: N806
    T_max = T_max or np.max(T)  # noqa: N806

    spc_fit = spc.model_copy()
    spc_fit.therm = Nasa7ThermFit.fit(
        T=T,
        Cp=Cp,
        S=S,
        H=H,
        formula=therm.formula,
        charge=therm.charge,
        T_min=T_min,
        T_mid=T_mid,
        T_max=T_max,
    )
    return spc_fit


# Display
def display(  # noqa: PLR0913
    spc: Species,
    *,
    props: Sequence[Literal["Cv", "Cp", "S", "H", "dH"]] = ("Cp", "S", "H"),
    others: Sequence[Species] = (),
    others_labels: Sequence[str] = (),
    T_range: tuple[float, float] = (200, 3000),  # noqa: N803
    units: UnitsData | None = None,
    label: str = "This work",
    x_label: str = "𝑇",  # noqa: RUF001
    y_labels: Sequence[str | None] | None = None,
    horizontal: bool = False,
) -> alt.Chart:
    """Display as an Arrhenius plot, optionally comparing to other rates.

    :param spc: Species thermo
    :param props: Thermodynamic properties to display
    :param others: Other reaction rates for comparison
    :param others_labels: Labels for other reaction rates
    :param t_range: Temperature range
    :param p: Pressure
    :param units: Units
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :return: Chart
    """
    return spc.therm.display(
        props=props,
        others=[o.therm for o in others],
        others_labels=others_labels,
        T_range=T_range,
        units=units,
        label=label,
        x_label=x_label,
        y_labels=y_labels,
        horizontal=horizontal,
    )


# Helpers
def pac99_input_line(  # noqa: PLR0913
    key: str | None = None,
    label1: str | None = None,
    num1: float | None = None,
    label2: str | None = None,
    num2: float | None = None,
    label3: str | None = None,
    num3: float | None = None,
    label4: str | None = None,
    num4: float | None = None,
    decimals: int = 3,
) -> str:
    """Format a line of PAC99 input."""
    widths = [6] + [6, 12] * 4
    vals = [
        key,
        label1,
        pac99_number(num1, decimals=decimals),
        label2,
        pac99_number(num2, decimals=decimals),
        label3,
        pac99_number(num3, decimals=decimals),
        label4,
        pac99_number(num4, decimals=decimals),
    ]
    line = ""
    for idx, (width, val) in enumerate(zip(widths, vals, strict=True)):
        if val is not None:
            indent = sum(widths[:idx])
            line = f"{line: <{indent}}{val: <{width}}"
    return line


def pac99_number(num: float | None = None, decimals: int = 3) -> str | None:
    """Format number for PAC99 input."""
    if num is None:
        return None

    num_str = f"{num:.{decimals}f}"
    num_str = num_str if "." in num_str else f"{num_str}."
    return f"{num_str: >10}  "
