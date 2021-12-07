"""
   Periodic Table operations
"""

from qcelemental import periodictable
from qcelemental import vdwradii


def to_symbol(atom):
    """ Obtain the atomic symbol for a given atom.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: str
    """
    return periodictable.to_E(atom)


def to_number(atom):
    """ Obtain the atomic number for a given atom.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: int
    """
    return periodictable.to_Z(atom)


def to_mass_number(atom):
    """ Obtain the mass number for a given atom.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: int
    """
    return periodictable.to_A(atom)


def to_mass(atom):
    """ Obtain the atomic mass for a given atom (in amu).

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: float
    """
    return periodictable.to_mass(atom)


def to_group(atom):
    """ Obtain the number of the group the atom belongs to.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: int
    """
    return periodictable.to_group(atom)


def van_der_waals_radius(symb):
    """ Obtain the van der Waals radius for an atom (in Bohr).

        :param symb: atomic symbol
        :type symb: str
        :rtype: float
    """
    return vdwradii.get(symb, units='bohr')
