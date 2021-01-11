"""
   Periodic Table operations
"""

from qcelemental import periodictable


def to_symbol(atom):
    return periodictable.to_E(atom)


def to_number(atom):
    return periodictable.to_Z(atom)


def to_mass(atom):
    return periodictable.to_mass(atom)
