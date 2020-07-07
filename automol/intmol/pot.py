"""
  Various potential interaction forms
"""

import numpy


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
    ('O', 'O'): [7.69e4, 3.6, 460.0, 0.8]
}


def pairwise_potential_matrix(geo):
    """ Generate a matrix of pairwise potentials
    """

    # Initialize matrix
    natoms = len(geo)
    pot_mat = numpy.zeros((natoms, natoms))

    for i in range(natoms):
        for j in range(natoms):
            pot_mat[i, j] = _pairwise_potentials(geo, (i, j))

    return pot_mat


def _pairwise_potentials(geo, idx_pair, potential='exp6'):
    """ Calculate the sum of the pairwise potential for a
        given set of atom pairs
    """

    # Get the indexes and symbols
    idx1, idx2 = idx_pair
    if idx1 != idx2:

        # Get the symbols of the atoms
        symbols = automol.geom.symbols(geo)
        symb1, symb2 = symbols[idx1], symbols[idx2]

        # Calculate interatomic distance
        rdist = automol.geom.distance(geo, idx1, idx2) * phycon.BOHR2ANG

        # Calculate the interaction potential value
        if potential == 'exp6':
            pot_val = _pairwise_exp6_potential(rdist, symb1, symb2)
        elif potential == 'lj_12_6':
            pot_val = _pairwise_lj_potential(rdist, symb1, symb2)
        else:
            pot_val = None

    else:
        pot_val = 1e10

    return pot_val


def _pairwise_exp6_potential(rdist, symb1, symb2):
    """ Calculate pot
    """

    exp6_params = EXP6_DCT.get((symb1, symb2), None)
    if exp6_params is None:
        exp6_params = EXP6_DCT.get((symb2, symb1), None)

    pot_val = _exp6_potential(rdist, *exp6_params)

    return pot_val


def _exp6_potential(rdist, apar, bpar, cpar, rcut):
    """ Calculate modified Buckhingham potential
    """
    if rdist < rcut:
        pot_val = apar * numpy.exp(-1.0*bpar*rcut) - (cpar / rcut**6)
    else:
        pot_val = apar * numpy.exp(-1.0*bpar*rdist) - (cpar / rdist**6)
    return pot_val


def _pairwise_lj_potential(rdist, symb1, symb2):
    """ Calculate pot
    """

    ljparams = LJ_DCT.get((symb1, symb2), None)
    if ljparams is None:
        ljparams = LJ_DCT.get((symb2, symb1), None)

    pot_val = _lj_potential(rdist, *ljparams)

    return pot_val


def _lj_potential(rdist, eps, sig):
    """ Calculate Lennard-Jones Potential
    """
    return (4.0 * eps) * ((sig / rdist)**12 - (sig / rdist)**6)


def _generate_pairs(geo, pairs='offdiag'):
    """ Generate a list of pairs to calculate potentials
    """

    assert pairs == 'offdiag', (
        'Can only generate list of off-diagonal pairs'
    )

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

    # Loop over neighbor dict to build pairs
    # print('geo', geo)
    # print('heavy idxs', heavy_idxs)
    # print('heacy pairs', heavy_pairs)
    # print('h idxs', h_idxs)
    # print('h pairs', h_pairs)
    # print(neigh_dct)
    #

    pairs = tuple()
    for i in range(len(geo)):
        for j in range(len(geo)):
            if i != j:
                pairs += ((i, j),)

    return pairs
