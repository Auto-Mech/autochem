"""
  Repulsion calculations
"""

from automol.zmat import geometry
from automol.intmol._pot import pairwise_potential_matrix


def low_repulsion_struct(zma_ref, zma_samp, pairs='offdiag', thresh=40.0):
    """ Check if the long-range energy for the sample structure
    exceeds that for the reference structure by more than the thresh
    """

    # # Convert to geoms
    geo_ref = geometry(zma_ref, remove_dummy_atoms=True)
    geo_samp = geometry(zma_samp, remove_dummy_atoms=True)

    # Calculate the pairwise potentials
    pot_mat = pairwise_potential_matrix(geo_ref)
    pot_mat_samp = pairwise_potential_matrix(geo_samp)

    # Generate the pairs for the potentials
    pairs = _generate_pairs(geo_ref, pairs=pairs)

    # Calculate sum of potentials
    sum_ref, sum_samp = 0.0, 0.0
    for (idx1, idx2) in pairs:
        sum_ref += pot_mat[idx1, idx2]
        sum_samp += pot_mat_samp[idx1, idx2]

    print('long_range_pots {:.2f} {:.2f} {:.2f}'.format(
        sum_ref, sum_samp, sum_samp-sum_ref))

    # # Check if the potentials are within threshold
    low_repulsion = bool((sum_samp - sum_ref) <= thresh)

    return low_repulsion


def _generate_pairs(geo, pairs='offdiag'):
    """ Generate a list of pairs to calculate potentials
    """

    assert pairs == 'offdiag', (
        'Can only generate list of off-diagonal pairs'
    )

    # Off-Diagonal pairs
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
