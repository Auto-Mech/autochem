"""
  Various helper functions for dealing with graphs
"""


# Handle the formatting of various lists
def atom_idx_to_symb(idxs, idx_symb_dct):
    """ Convert a list of atom idxs (a1, a2, ..., an)
        to atom symbols
    """
    return tuple(idx_symb_dct[idx] for idx in idxs)


def bond_idx_to_symb(idxs, idx_symb_dct):
    """ Convert a list of bond idxs ((a1, b1), (a2, b2), ..., (an, bn))
        to pairs of atom symbols
    """
    return tuple(
        (idx_symb_dct[idx1], idx_symb_dct[idx2]) for (idx1, idx2) in idxs
    )


def filter_idxs(idxs_lst, filterlst=()):
    """ Filter out a tuple
    """

    filtered_lst = tuple()

    for idxs in idxs_lst:
        if not any(set(idxs) <= set(fidxs) for fidxs in filterlst):
            filtered_lst += (idxs,)

    return filtered_lst


def ring_idxs(gra_rings):
    """ Get idxs for rings
    """

    _ring_idxs = tuple()

    for ring in gra_rings:
        idxs = tuple(ring[0].keys())
        if idxs:
            _ring_idxs += (idxs,)

    return _ring_idxs
