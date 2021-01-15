"""
  Various potential interaction forms
"""

import numpy
from phydat import phycon
from automol.geom import symbols, distance


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


# change to a dictionary
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


# Check the repulsion
def low_repulsion_struct(geo_ref, geo_samp,
                         potential='exp6', pairs='offdiag', thresh=40.0,):
    """ Check if the long-range interaction energy for the sample structure
        exceeds that for the reference structure by more than given threshold.

        :param geo_ref: reference structure against which repulsion is assessed
        :type geo_ref: automol geometry data structure
        :param geo_samp: test geometry to assess repulsion
        :type geo_samp: automol geometry data structure
        :param pairs: pair choosing algorithm
        :type pairs: str
        :param thresh: threshold for determining level of repulsion (kcal/mol)
        :type thesh: float
        :rtype: bool
    """

    # Calculate the pairwise potentials
    pot_mat = pairwise_potential_matrix(geo_ref, potential=potential)
    pot_mat_samp = pairwise_potential_matrix(geo_samp, potential=potential)

    # Generate the pairs for the potentials
    pairs = _generate_pairs(geo_ref, pairs=pairs)

    # Calculate sum of potentials
    sum_ref, sum_samp = 0.0, 0.0
    for (idx1, idx2) in pairs:
        sum_ref += pot_mat[idx1, idx2]
        sum_samp += pot_mat_samp[idx1, idx2]

    # Check if the potentials are within threshold
    low_repulsion = bool((sum_samp - sum_ref) <= thresh)

    return low_repulsion


def _generate_pairs(geo):
    """ Determine all of the pairs of atoms to calculate
        interatomic potentials for. Only generate off-diagonal pairs for now.

        :param geo: automol geometry object
        :type geo: tuple(tuple(float))
        :param pairs: pair choosing algorithm
        :type pairs: str
    """

    natoms = automol.geom.count(geo) 

    pairs = tuple()
    for i in range(natoms):
        for j in range(natoms):
            if i != j:
                pairs += ((i, j),)

    return pairs
