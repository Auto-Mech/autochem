"""
  Various potential interaction forms
"""

import numpy
from qcelemental import constants as qcc
from automol.geom import symbols, distance


# DCTS OF POTENTIAL PARAMETERS
# eps[whatever], sig[ang] params
LJ_DCT = {
    ('H', 'H'): [0.25, 1.0],
    ('H', 'C'): [0.25, 1.0],
    ('H', 'O'): [0.25, 1.0],
    ('C', 'C'): [0.25, 1.0],
    ('C', 'O'): [0.25, 1.0],
    ('O', 'O'): [0.25, 1.0],
}

# A, B, C params E[kcal] R[Ang]; R cutoff
EXP6_DCT = {
    ('H', 'H'): [2.442e3, 3.74, 48.8, 1.0],
    ('H', 'C'): [6.45e3, 3.67, 116.0, 1.0],
    ('H', 'O'): [6.45e3, 3.67, 116.0, 1.0],
    ('C', 'C'): [7.69e4, 3.6, 460.0, 0.8],
    ('C', 'O'): [7.69e4, 3.6, 460.0, 0.8],
    ('O', 'O'): [7.69e4, 3.6, 460.0, 0.8],
    ('Cl', 'Cl'): [0.0, 3.6, 460.0, 0.8],
    ('Cl', 'C'): [0.0, 3.6, 460.0, 0.8],
    ('Cl', 'O'): [0.0, 3.6, 460.0, 0.8],
    ('Cl', 'H'): [0.0, 3.6, 460.0, 0.8]
}


def _read_params(dct, symb1, symb2):
    """ Calculate pot
    """

    params = dct.get((symb1, symb2), None)
    if params is None:
        params = dct.get((symb2, symb1), None)

    return params


# POTENTIAL FORMS
def lj_potential(rdist, eps, sig):
    """ Calculate Lennard-Jones Potential
    """
    return (4.0 * eps) * ((sig / rdist)**12 - (sig / rdist)**6)


def exp6_potential(rdist, apar, bpar, cpar, rcut):
    """ Calculate modified Buckhingham potential
    """
    if rdist < rcut:
        pot_val = apar * numpy.exp(-1.0*bpar*rcut) - (cpar / rcut**6)
    else:
        pot_val = apar * numpy.exp(-1.0*bpar*rdist) - (cpar / rdist**6)
    return pot_val


# INTERACTION MATRIX WITH ABOVE POTENTIALS
def pairwise_potential_matrix(geo, potential='exp6'):
    """ Generate a matrix of pairwise potentials
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
        given set of atom pairs
    """

    # Get the indexes and symbols
    idx1, idx2 = idx_pair
    if idx1 != idx2:

        # Get the symbols of the atoms
        symbs = symbols(geo)
        symb1, symb2 = symbs[idx1], symbs[idx2]

        # Calculate interatomic distance
        rdist = (
            distance(geo, idx1, idx2) *
            qcc.conversion_factor('bohr', 'angstrom')
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
        pot_val = 1e10

    return pot_val


def _generate_pairs(geo, pairs='offdiag'):
    """ Generate a list of pairs to calculate potentials
    """

    assert pairs == 'offdiag', (
        'Can only generate list of off-diagonal pairs'
    )

    # OFF-DIAG PAIRS
    pairs = tuple()
    for i in range(len(geo)):
        for j in range(len(geo)):
            if i != j:
                pairs += ((i, j),)

    # ALTERNATE IDX CODE
    # Grab the indices of the heavy atoms for the zmas
    # heavy_idxs = automol.geom.atom_indices(geo, 'H', match=False)

    # Get the h atom idxs as list of list for each heavy atom
    # gra = automol.geom.graph(geo)
    # neigh_dct = automol.graph.atom_neighbor_keys(gra)
    # h_idxs = ()
    # for idx in heavy_idxs:
    #     neighs = neigh_dct[idx]
    #     h_idxs += (tuple(x for x in neighs if x not in heavy_idxs),)

    # Get heavy atom pairs
    # heavy_pairs = tuple(itertools.combinations(heavy_idxs, 2))
    # h_pairs = tuple()
    # for comb in itertools.combinations(h_idxs, 2):
    #     h_pairs += tuple(itertools.product(*comb))

    return pairs
