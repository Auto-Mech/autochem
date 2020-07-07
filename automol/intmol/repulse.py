"""
  Repulsion calculations
"""

def _low_repulsion_struct(zma_ref, zma_samp, thresh=40.0):
    """ Check if the long-range energy for the sample structure
    exceeds that for the reference structure by more than the thresh
    """

    # # Convert to geoms
    geo_ref = automol.zmatrix.geometry(zma_ref)
    geo_samp = automol.zmatrix.geometry(zma_samp)

    # Calculate the pairwise potentials
    pot_mat = _pairwise_potential_matrix(geo_ref)
    pot_mat_samp = _pairwise_potential_matrix(geo_samp)

    # Generate the pairs for the potentials
    pairs = _generate_pairs(geo_ref)

    # Calculate sum of potentials
    sum_ref, sum_samp = 0.0, 0.0
    for (idx1, idx2) in pairs:
        sum_ref += pot_mat[idx1, idx2]
        sum_samp += pot_mat_samp[idx1, idx2]

    print('long_range_pots {:.2f} {:.2f} {:.2f}'.format(
        sum_ref, sum_samp, sum_samp-sum_ref))

    # # Check if the potentials are within threshold
    low_repulsion = bool((sum_samp - sum_ref) <= thresh)
    # low_repulsion = True

    return low_repulsion




