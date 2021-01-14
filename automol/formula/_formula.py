""" Library for dealing with molecular formula,
    represented as dict[atom symbol: atom number]
"""

import functools
import itertools
import collections
from phydat import ptab


def electron_count(fml):
    """ Count the number of electrons for the atoms in a molecular formula.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: int
    """

    assert _is_standard(fml)

    elec_count = 0
    for key in fml:
        value = fml[key]
        elec_count += value*ptab.to_number(key)

    return elec_count


def atom_count(fml):
    """ Count the number of atoms in this molecular formula.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: int
    """

    assert _is_standard(fml)

    return sum(fml.values())


def element_count(fml, symb):
    """ Count the number of a given element in this molecular formula.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :param symb: atomic symbol of element to be counted
        :type symb: str
        :rtype: int
    """

    assert _is_standard(fml)

    return fml[symb] if symb in fml else 0


def add_element(fml, symb, num=1):
    """ add or subtract (if num < 0) this element from the molecular formula

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :param symb: atomic symbol of element to be added
        :type symb: str
        :param num: number of the element to add to the formula
        :type num: int
        :rtype: dict[str:int]
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


def join(fml1, fml2):
    """ Join two formulas together.

        :param fml1: stochiometric chemical formula 1
        :type fml1: dict[str:int]
        :param fml2: stochiometric chemical formula 2
        :type fml2: dict[str:int]
        :rtype: int
    """

    fml = dict(fml1)
    for symb, num in fml2.items():
        fml = add_element(fml, symb, num=num)

    return fml


def join_sequence(fmls):
    """ Join a sequence of formulas together:

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: int
    """

    return functools.reduce(join, fmls)


def string(fml):
    """ Convert formula dictionary to formula string in the Hill convention. 
        Resultant string is identical to InChI formula string.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: str
    """
    
    fml_lst = [(symb, fml[symb]) for symb in sorted_symbols(fml.keys())]

    fml_str = ''.join(map(
        str,
        itertools.filterfalse(lambda x: x == 1, itertools.chain(*fml_lst))))

    return fml_str


def string2(fml):
    """ Convert formula dictionary to formula string that includes 1s in when there
        is only one atom.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: str
    """

    fml = collections.OrderedDict(sorted(fml.items()))

    fml_str = ''.join(map(str, itertools.chain.from_iterable(fml.items())))

    return fml_str


def sorted_symbols(seq, symbs=('C', 'H')):
    """ Produce a sorted list of atomic symbols; some elements given priority.
        Be default, C placed first, then H, others follow in alphabetical order.

        :param seq: formula or sequence of atomic symbols
        :type seq: dict, list, or tuple
        :param symbs: atomic symbols to place first
        :type symbs: sequence of strings
        :rtyp: tuple(str)
    """

    def _sort_key(char):
        val = symbs.index(char) if char in symbs else len(symbs)
        return (val, char)

    return tuple(sorted(seq, key=_sort_key))


def argsort_symbols(seq, symbs=('C', 'H')):
    """ Get the sort order for a sequence of atomic symbols.

        :param seq: formula or sequence of atomic symbols
        :type seq: dict, list, or tuple
        :param symbs: atomic symbols to place first
        :type symbs: sequence of strings
        :rtyp: tuple(str)
    """

    def _sort_key(entry):
        char = entry[0]
        rest = entry[1:]
        val = symbs.index(char) if char in symbs else len(symbs)
        return (val, char, rest)

    return tuple(idx for (val, idx) in
                 sorted(((v, i) for (i, v) in enumerate(seq)), key=_sort_key))


def _is_standard(fml):
    """ Assess if the formula conforms to the standard form.

        :param fml: stochiometric chemical formula
        :type fml: dict[str:int]
        :rtype: bool
    """

    symbs = list(fml.keys())

    return symbs == list(filter(ptab.to_number, map(ptab.to_symbol, symbs)))
