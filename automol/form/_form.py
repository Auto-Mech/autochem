"""Library for dealing with molecular formula,
represented as dict[atom symbol: atom number].
"""

import collections
import functools
import itertools

import pyparsing as pp
from _collections_abc import Sequence
from pyparsing import pyparsing_common as ppc

from phydat import ptab

from ..util import dict_

SYMBOL = pp.Word(pp.alphas.upper(), pp.alphas.lower())
STOICH = pp.Opt(pp.Literal("*") | ppc.integer).setParseAction(lambda x: x if x else [1])
ATOM_COUNT = pp.Group(SYMBOL + STOICH)
FORMULA = pp.OneOrMore(ATOM_COUNT)
Formula = dict[str, int]


def electron_count(fml: Formula) -> int:
    """Count the number of electrons for the atoms in a molecular formula.

    :param fml: Stochiometric chemical formula
    :return: Number of electrons.
    """
    assert _is_standard(fml)

    elec_count = 0
    for key in fml:
        value = fml[key]
        elec_count += value * ptab.to_number(key)

    return elec_count


def atom_count(fml: Formula) -> int:
    """Count the number of atoms in this molecular formula.

    :param fml: Stochiometric chemical formula
    :return: Number of atoms.
    """
    assert _is_standard(fml)

    return sum(fml.values())


def heavy_atom_count(fml: Formula) -> int:
    """Count the number of heavy atoms in this molecular formula.

    :param fml: Stochiometric chemical formula
    :return: Number of heavy atoms
    """
    assert _is_standard(fml)
    fml = without(fml, ["H"])
    return sum(fml.values())


def element_count(fml: Formula, symb: str) -> int:
    """Count the number of a given element in this molecular formula.

    :param fml: Stochiometric chemical formula
    :param symb: Atomic symbol of element to be counted
    :return: Number of a certain element in the formula
    """
    assert _is_standard(fml)

    return fml[symb] if symb in fml else 0


def without(fml: dict[str, int], symbs: tuple = ()) -> dict[str, int]:
    """Return a formula without hydrogen.

    :param fml: Chemical formula, without hydrogen
    :symbs: Chemical symbols
    :return: Dictionary with new formula, without hydrogen
    """
    return {k: v for k, v in fml.items() if k not in symbs}


def match(fml1: dict[str, int], fml2: dict[str, int]) -> bool:
    """Check for a match between two formulas, allowing wildcard values.

    A stoichiometry of -1 indicates a wildcard value

    :param fml1: A chemical formula
    :param fml2: Another chemical formula
    :return: `True` if so, `False` if not
    """
    excl_symbs1 = dict_.keys_by_value(fml1, lambda x: x < 0)
    excl_symbs2 = dict_.keys_by_value(fml2, lambda x: x < 0)
    excl_symbs = excl_symbs1 | excl_symbs2

    fml1 = without(fml1, excl_symbs)
    fml2 = without(fml2, excl_symbs)
    return fml1 == fml2


def add_element(fml: Formula, symb: str, num: int = 1) -> Formula:
    """Add or subtract (if num < 0) this element from the molecular formula.

    :param fml: Stochiometric chemical formula
    :param symb: Atomic symbol of element to be added
    :param num: Number of the element to add to the formula
    :return: Formula with added element
    """
    assert ptab.to_number(symb)
    assert _is_standard(fml)

    symb = ptab.to_symbol(symb)
    fml = fml.copy()
    if symb in fml:
        fml[symb] += num
    else:
        fml[symb] = num

    assert fml[symb] > 0

    return fml


def join(fml1: Formula, fml2: Formula) -> int:
    """Join two formulas together.

    :param fml1: Stochiometric chemical formula 1
    :param fml2: Stochiometric chemical formula 2
    :return: Formula with the sum of both formulas
    """
    fml = dict(fml1)
    for symb, num in fml2.items():
        fml = add_element(fml, symb, num=num)

    return fml


def join_sequence(fmls: Formula) -> int:
    """Join a sequence of formulas together.

    :param fml: Stochiometric chemical formula
    :return: Sum of the formulas
    """
    return functools.reduce(join, fmls)


def sorted_symbols_in_sequence(fmls: list[dict]) -> list[dict]:
    """Sort a sequence of formulas based on Hill-sorting.

    :param fmls: A sequence of formulas
    :return: The same sequence, but sorted
    """
    return sorted_symbols(join_sequence(fmls).keys())


# Str<->Dict Converters
def string(fml: Formula, hyd: bool = True) -> str:
    """Convert formula dictionary to formula string in the Hill convention.
    Resultant string is identical to InChI formula string.

    :param fml: stochiometric chemical formula
    :param hyd: include hydrogens?
    :return: True if formula includes hydrogen, False if no hydrogen
    """
    fml_lst = [
        (symb, fml[symb]) for symb in sorted_symbols(fml.keys()) if symb != "H" or hyd
    ]

    fml_str = "".join(
        map(str, itertools.filterfalse(lambda x: x == 1, itertools.chain(*fml_lst)))
    )

    return fml_str


def string2(fml: Formula) -> str:
    """Convert formula dictionary to formula string that includes 1s in when
    there is only one atom.

    :param fml: stochiometric chemical formula
    :return: Formula which includes 1s when there is only one atom
    """
    fml = collections.OrderedDict(sorted(fml.items()))

    fml_str = "".join(map(str, itertools.chain.from_iterable(fml.items())))

    return fml_str


def from_string(fml_str: str) -> Formula:
    """Convert formula string to formula dictionary.

    Wildcard values can be specified as, for example, 'O2H*', which will be interpreted
    as {'O': 2, 'H': -1} where the -1 indicates a wildcard stoichiometry.

    :param fml_str: stochiometric chemical formula string
    :return: The formula
    """
    fml = dict(FORMULA.parseString(fml_str).asList())
    fml = dict_.transform_values(fml, lambda x: -1 if x == "*" else x)
    return fml


def sorted_symbols(
    seq: Sequence,
    symbs_first: Sequence[str] = ("C", "H"),
    symbs_last: Sequence[str] = (),
) -> tuple[str, ...]:
    """Produce a sorted list of atomic symbols; some elements given priority.
    By default, C placed first, then H, then others in alphabetical order.

    :param seq: Formula or sequence of atomic symbols
    :param symbs_first: Atomic symbols to place first
    :param symbs_last: Atomic symbols to place last
    :return: Sorted list
    """

    def _sort_key(char):
        if char in symbs_first:
            val = symbs_first.index(char)
        elif char in symbs_last:
            val = len(symbs_first) + 1 + symbs_last.index(char)
        else:
            val = len(symbs_first)
        return (val, char)

    return tuple(sorted(seq, key=_sort_key))


def argsort_symbols(
    seq: Sequence,
    symbs_first: Sequence[str] = ("C", "H"),
    symbs_last: Sequence[str] = (),
    idx: int | None = None,
) -> tuple[str, ...]:
    """Determine the sort order for a sequence of atomic symbols.

    :param seq: Formula or sequence of atomic symbols
    :param symbs_first: Atomic symbols to place first
    :param symbs_last: Atomic symbols to place last
    :param idx: Index of symbol for sorting
    :return: Sorted syymboles
    """

    def _sort_key(entry):
        if idx is not None:
            entry = tuple(entry[0]) + entry[1:]
            start = entry[:idx]
            char = entry[idx]
            rest = entry[(idx + 1) :]
        else:
            start = ()
            char = entry[0]
            rest = entry[1:]

        if char in symbs_first:
            val = symbs_first.index(char)
        elif char in symbs_last:
            val = len(symbs_first) + 1 + symbs_last.index(char)
        else:
            val = len(symbs_first)
        return (start, val, char, rest)

    return tuple(
        idx
        for (val, idx) in sorted(((v, i) for (i, v) in enumerate(seq)), key=_sort_key)
    )


def sort_vector(
    fml: Sequence[str], symbs: Sequence[str] | None = None
) -> tuple[int, ...]:
    """Generate a sort vector for sorting various formulas against each other.

    :param fml_str: stochiometric chemical formula string
    :param symbs: atomic symbols in the desired sort order (optional)
    :return: Sorted vector
    """
    if symbs is None:
        symbs = sorted_symbols(fml.keys())

    vec = tuple(fml[s] if s in fml else 0 for s in symbs)
    return vec


def _is_standard(fml: Formula) -> bool:
    """Assess if the formula conforms to the standard form.

    :param fml: stochiometric chemical formula
    :return: True if formula is in standard form, false if not
    """
    symbs = list(fml.keys())

    return symbs == list(filter(ptab.to_number, map(ptab.to_symbol, symbs)))
