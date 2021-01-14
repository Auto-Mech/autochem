"""
  Various potential interaction forms
"""

import numpy
from phydat import phycon
from automol.geom import symbols, distance


# DCTS OF POTENTIAL PARAMETERS
# eps[whatever], sig[ang] params
LJ_DCT = {
    ('H', 'H'): (0.25, 1.0),
    ('H', 'C'): (0.25, 1.0),
    ('H', 'O'): (0.25, 1.0),
    ('C', 'C'): (0.25, 1.0),
    ('C', 'O'): (0.25, 1.0),
    ('O', 'O'): (0.25, 1.0)
}

# A, B, C params E[kcal] R[Ang]; R cutoff
EXP6_DCT = {
    ('H', 'H'): (2.442e3, 3.74, 48.8, 1.0),
    ('H', 'C'): (6.45e3, 3.67, 116.0, 1.0),
    ('H', 'N'): (6.45e3, 3.67, 116.0, 1.0),
    ('H', 'O'): (6.45e3, 3.67, 116.0, 1.0),
    ('C', 'C'): (7.69e4, 3.6, 460.0, 0.8),
    ('C', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('O', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('Cl', 'Cl'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'C'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'O'): (0.0, 3.6, 460.0, 0.8),
    ('Cl', 'H'): (0.0, 3.6, 460.0, 0.8),
    ('N', 'N'): (7.69e4, 3.6, 460.0, 0.8),
    ('N', 'O'): (7.69e4, 3.6, 460.0, 0.8),
    ('N', 'C'): (7.69e4, 3.6, 460.0, 0.8)
}


def _read_params(dct, symb1, symb2):
    """ Read the parameters from one of the potential
        parameter dictionaries

        :param dct: potential parameter dct
        :type dct: dict[(symbs):(params)]
        :param symb1: atomic symbol of atom1 involved in the interaction
        :param symb1: str
        :param symb2: atomic symbol of atom2 involved in the interaction
        :param symb2: str
        :rtype: tuple(float)
    """

    params = dct.get((symb1, symb2), None)
    if params is None:
        params = dct.get((symb2, symb1), None)

    return params


# POTENTIAL FORMS
def lj_potential(rdist, eps, sig):
    """ Calculate potential energy value of two interacting bodies
        assuming a 12-6 Lennard-Jones potential.

        :param rdist: distance between two interacting bodies (Bohr)
        :type rdist: float
        :param eps: Lennard-Jones epsilon parameter for interaction (_)
        :type eps: float
        :param eps: Lennard-Jones sigma parameter for interaction (_)
        :type eps: float
        :rtpe: float
    """
    return (4.0 * eps) * ((sig / rdist)**12 - (sig / rdist)**6)


def exp6_potential(rdist, apar, bpar, cpar, rcut):
    """ Calculate potential energy value of two interacting bodies
        assuming a modified Buckingham potential.

        :param rdist: distance between two interacting bodies (Bohr)
        :type rdist: float
        :param apar: potential parameter A
        :type apar: float
        :param bpar: potential parameter B
        :type bpar: float
        :param cpar: potential parameter C
        :type cpar: float
        :param rcut: threshhold where interaction potential becomes constant
        :type rcut: float
        :rtpe: float
    """
    if rdist < rcut:
        pot_val = apar * numpy.exp(-1.0*bpar*rcut) - (cpar / rcut**6)
    else:
        pot_val = apar * numpy.exp(-1.0*bpar*rdist) - (cpar / rdist**6)
    return pot_val


# INTERACTION MATRIX WITH ABOVE POTENTIALS
def pairwise_potential_matrix(geo, potential='exp6'):
    """ Build a matrix of values describing the interaction potential
        of all atoms in a geometry.

        :param geo: automol geometry object
        :type geo: tuple(tuple(float))
        :param potential: potential model of the atomic interactions
        :type potential: str
        :rtype: nd.array
    """

    # Initialize matrix
    natoms = len(geo)
    pot_mat = numpy.zeros((natoms, natoms))

    for i in range(natoms):
        for j in range(natoms):
            pot_mat[i, j] = _pairwise_potentials(
                geo, (i, j), potential=potential)

    return pot_mat


def _pairwise_potentials(geo, idx_pair, potential='exp6'):
    """ Calculate the sum of the pairwise potential for a
        given pair of atoms in a geometry.

        :param geo: automol geometry object
        :type geo: tuple(tuple(float))
        :param idx_pair: indices of atoms for which to calculate interaction
        :type idx_pair: tuple(int, int)
        :param potential: potential model of the atomic interactions
        :type potential: str
        :rtype: nd.array
    """

    # Get the indexes and symbols
    idx1, idx2 = idx_pair
    if idx1 != idx2:

        # Get the symbols of the atoms
        symbs = symbols(geo)
        symb1, symb2 = symbs[idx1], symbs[idx2]

        # Calculate interatomic distance
        rdist = (
            distance(geo, idx1, idx2) * phycon.BOHR2ANG
        )

        # Calculate the potential
        if potential == 'exp6':
            params = _read_params(EXP6_DCT, symb1, symb2)
            pot_val = exp6_potential(rdist, *params)
        elif potential == 'lj_12_6':
            params = _read_params(LJ_DCT, symb1, symb2)
            pot_val = lj_potential(rdist, *params)
        else:
            pot_val = None

    else:
        pot_val = 1.0e10

    return pot_val
