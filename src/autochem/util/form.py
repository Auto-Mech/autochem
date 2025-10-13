"""Formula utilities."""

from collections import defaultdict
from typing import TypeAlias

import pyparsing as pp
from pyparsing import pyparsing_common as ppc

SYMBOL = pp.Word(pp.alphas.upper(), pp.alphas.lower())
COUNT = pp.Opt(ppc.integer).setParseAction(lambda x: x if x else [1])
STOICH = pp.Group(SYMBOL + COUNT)
FORMULA = pp.OneOrMore(STOICH)
Formula: TypeAlias = dict[str, int]
FormulaData: TypeAlias = Formula | str


def normalize_input(fml: FormulaData) -> Formula:
    """Normalize formula input, which could be given as a string.

    :param fml_inp: Formula input
    :return: Formula
    """
    if isinstance(fml, str):
        fml_lst: list[tuple[str, int]] = FORMULA.parse_string(fml).as_list()
        fml = defaultdict(int)
        for sym, count in fml_lst:
            fml[sym.title()] += int(count)

    return {k.title(): int(v) for k, v in fml.items() if v}


def string(fml_inp: FormulaData, *, ones: bool = False, upper: bool = False) -> str:
    """Convert formula input to string.

    :param fml_inp: Formula input
    :param ones: Whether to write 1 for single atoms
    :param upper: Whether to write all-uppercase symbols
    :return: Formula string
    """
    fml = normalize_input(fml_inp)
    fml = {k.upper(): v for k, v in fml.items()} if upper else fml
    return "".join([f"{k}{v}" if v > 1 or ones else k for k, v in fml.items()])
