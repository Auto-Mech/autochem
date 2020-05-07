""" molecular formula
"""
import itertools
import collections
from qcelemental import periodictable as pt


def electron_count(fml):
    """ the number of atoms in this molecular formula
    """
    assert _is_standard(fml)
    elec_count = 0
    for key in fml:
        value = fml[key]
        elec_count += value*pt.to_Z(key)
    return elec_count


def atom_count(fml):
    """ the number of atoms in this molecular formula
    """
    assert _is_standard(fml)
    return sum(fml.values())


def element_count(fml, sym):
    """ the number of a given element in this molecular formula
    """
    assert _is_standard(fml)
    return fml[sym] if sym in fml else 0


def hydrogen_count(fml):
    """ the number of hydrogens in this molecular formula
    """
    return element_count(fml, 'H')


def add_element(fml, sym, num=1):
    """ add or subtract (if num < 0) this element from the molecular formula
    """
    assert pt.to_Z(sym)
    assert _is_standard(fml)
    sym = pt.to_E(sym)
    fml = fml.copy()
    if sym in fml:
        fml[sym] += num
    else:
        fml[sym] = num
    assert fml[sym] > 0
    return fml


def add_hydrogen(fml, num=1):
    """ add hydrogen to this molecular formula
    """
    return add_element(fml, 'H', num)


def join(fml1, fml2):
    """ join two formulas together
    """
    fml = dict(fml1)
    for sym, num in fml2.items():
        fml = add_element(fml, sym, num=num)
    return fml


def string(fml):
    """ convert formula dictionary to formula string in the Hill convention

    (should come out identical to an InChI formula string)
    """
    ncar = element_count(fml, 'C')
    nhyd = element_count(fml, 'H')

    fml.pop('C')
    fml.pop('H')

    fml_lst = list(sorted(fml.items()))
    if nhyd:
        fml_lst.insert(0, ('H', nhyd))

    if ncar:
        fml_lst.insert(0, ('C', ncar))

    fml_str = ''.join(map(str,
        itertools.filterfalse(lambda x: x==1, itertools.chain(*fml_lst))))

    return fml_str


def _is_standard(fml):
    syms = list(fml.keys())
    return syms == list(filter(pt.to_Z, map(pt.to_E, syms)))
