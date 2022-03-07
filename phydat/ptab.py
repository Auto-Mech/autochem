"""
   Periodic Table operations
"""

from qcelemental import periodictable
from qcelemental import vdwradii

GROUP_2_VALENCE = {
    None: 0,
    1: 1,   # H
    2: 2,   # Be
    13: 3,  # B
    14: 4,  # C
    15: 3,  # N
    16: 2,  # O
    17: 1,  # F
    18: 0,  # He
}

GROUP_2_LONE_PAIR_COUNT = {
    None: 0,
    1: 0,   # H
    2: 0,   # Be
    13: 0,  # B
    14: 0,  # C
    15: 1,  # N
    16: 2,  # O
    17: 3,  # F
    18: 4,  # He
}


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


def valence(atom):
    """ Obtain the number of bonds typically formed by the atom.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: int
    """
    grp = periodictable.to_group(atom)
    val = GROUP_2_VALENCE[grp] if grp in GROUP_2_VALENCE else None
    return val


def lone_pair_count(atom):
    """ Obtain the number of bonds typically formed by the atom.

        :param atom: atom representation (symbol, number, mass)
        :type atom: str/int
        :rtype: int
    """
    grp = periodictable.to_group(atom)
    lpc = (GROUP_2_LONE_PAIR_COUNT[grp]
           if grp in GROUP_2_LONE_PAIR_COUNT else None)
    return lpc


def van_der_waals_radius(symb):
    """ Obtain the van der Waals radius for an atom (in Bohr).

        :param symb: atomic symbol
        :type symb: str
        :rtype: float
    """
    return vdwradii.get(symb, units='bohr')
