"""Utility functions for reading and writing PAC99 data."""

import re

import more_itertools as mit
import pydantic
import pyparsing as pp
from pyparsing import common as ppc


# Thermo
class Pac99ThermoParseResults(pydantic.BaseModel):
    """PAC99 thermo parse results."""

    name: str
    formula: dict[str, int]
    ranges: list[tuple[float, float]]
    coeffs_lst: list[list[float]]
    identifier: str
    comment: str
    phase: int
    weight: float
    Hf_298: float


def parse_thermo(therm_str: str) -> Pac99ThermoParseResults:
    """Extract all thermo information from a PAC99 output string.

    :param therm_str: PAC99 output string
    :return: Parse results including name, formula, coefficients, etc.
    """
    line1, line2, *lines = therm_str.strip().splitlines()

    # Parse line 1
    name = line1[:24].strip()
    comment = line1[24:].strip()

    # Parse line 2
    count = int(line2[:2])
    identifier = line2[3:9].strip()
    formula = dict(FORM_ENTRIES.parse_string(line2[9:50]).as_list())
    phase = int(line2[51])
    weight = float(line2[52:65])
    Hf_298 = float(line2[65:80])  # noqa: N806

    # Parse remaining lines, only storing those that match the NASA-7 form
    coeff_head_ = (
        ppc.number(Key.min)
        + ppc.number(Key.max)
        + ppc.integer(Key.num_heat_coeffs)
        + pp.Group(ppc.number * 8)(Key.heat_exps)
        + ppc.number(Key.H298)
    )
    heat_exps_nasa = (0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0)
    coeff_ = pp.Combine(ppc.number + pp.Literal("D") + ppc.signed_integer)
    coeff_.set_parse_action(lambda toks: float(re.sub("[Dd]", "e", toks[0])))
    coeffs_ = coeff_ * 5

    ranges = []
    coeffs_lst = []
    for coeff_head, coeff_line1, coeff_line2 in mit.chunked(lines, 3, strict=True):
        # Parse the header
        res = coeff_head_.parse_string(coeff_head)
        num_heat_coeffs = res.get(Key.num_heat_coeffs)
        heat_exps = tuple(res.get(Key.heat_exps))

        # If it has the NASA-7 form, read in data
        if num_heat_coeffs == 5:
            assert heat_exps == heat_exps_nasa, f"{heat_exps} != {heat_exps_nasa}"

            # Get temperature range
            T_min = res.get(Key.min)  # noqa: N806
            T_max = res.get(Key.max)  # noqa: N806
            ranges.append((T_min, T_max))

            # Parse coefficients
            coeffs1 = coeffs_.parse_string(coeff_line1).as_list()
            coeffs2 = coeffs_.parse_string(coeff_line2).as_list()
            assert not any(coeffs2[:-2])
            coeffs_lst.append(coeffs1 + coeffs2[-2:])

    assert len(ranges) == len(coeffs_lst) <= count

    return Pac99ThermoParseResults(
        name=name,
        formula=formula,
        ranges=ranges,
        coeffs_lst=coeffs_lst,
        identifier=identifier,
        comment=comment,
        phase=phase,
        weight=weight,
        Hf_298=Hf_298,
    )


# Helpers
class Key:
    """Pyparsing token keys."""

    min = "min"
    max = "max"
    num_heat_coeffs = "num_heat_coeffs"
    heat_exps = "heat_exps"
    a5 = "a5"
    a6 = "a6"
    H298 = "H298"


#  - Pyparsing expressions
#   - Thermo entry formula
FORM_KEY = pp.OneOrMore(pp.Char(pp.alphas))
FORM_VAL = ppc.number.copy()
FORM_VAL.set_parse_action(lambda toks: int(toks[0]))
FORM_ENTRY = pp.Group(FORM_KEY + FORM_VAL)
FORM_ENTRIES = pp.OneOrMore(FORM_ENTRY)
#   - Thermo temperature line
