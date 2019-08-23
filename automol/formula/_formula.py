""" molecular formula
"""
from qcelemental import periodictable as pt


def electron_count(fml):
    """ the number of atoms in this molecular formula
    """
    assert _is_standard(fml)
    electron_count = 0
    for key in fml:
        value = fml[key]
        electron_count += value*pt.to_Z(key)
    return electron_count


def atom_count(fml):
    """ the number of atoms in this molecular formula
    """
    assert _is_standard(fml)
    return sum(fml.values())


def hydrogen_count(fml):
    """ the number of hydrogens in this molecular formula
    """
    assert _is_standard(fml)
    return fml['H']


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


def _is_standard(fml):
    syms = list(fml.keys())
    return syms == list(filter(pt.to_Z, map(pt.to_E, syms)))
