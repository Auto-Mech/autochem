""" drivers for coordinate scans
"""

import itertools


# Build the grirds ultimately used for building potentials
def coords(grids):
    """Determine the dimensions of the grid

    m = (m1, m2, m3)
    n = (n1, n2)
    p = mxn = ((m1, n1), (m1, n2), (m2, n1), (m2, n2), (m3, n1), (m3, n2))

    """
    assert len(grids) in (1, 2, 3, 4), "Rotor must be 1-4 dimensions"
    grid_vals = tuple(itertools.product(*grids))
    return grid_vals


# Manipulate potentials
def scale(pot, scale_factor):
    """Scale the potential by scaling factor

    :param pot: potential along a coordinate
    :type pot: dict[tuple(float)] = float
    :param scale_coeff: initial scaling coeffcient
    :type scael_coeff: float
    :param num_tors: number of torsions used in scaling
    :type num_tors: int
    :rtype:
    """

    new_pot = {}
    val = 0
    pos_der = 0
    for i, (idx, val) in enumerate(pot.items()):
        new_pot[idx] = val * scale_factor
        if i == 1:
            pos_der = val
    neg_der = val
    if (
        neg_der < 0.001
        or pos_der < 0.001
        or abs((neg_der - pos_der) / min(neg_der, pos_der)) > 0.4
    ):
        print(
            "WARNING: start and end derivatives of torsional potential don't match",
            neg_der,
            pos_der,
        )
    return new_pot


def relax_scale(pot):
    """Scale the potential by scaling factor

    :param pot: potential along a coordinate
    :type pot: dict[tuple(float)] = float
    :param scale_coeff: initial scaling coeffcient
    :type scael_coeff: float
    :param num_tors: number of torsions used in scaling
    :type num_tors: int
    :rtype:
    """

    new_pot = {}
    for idx, val in pot.items():
        # scale_factor = 1.0
        # scale_factor = 1.0/(1+0.05*val)
        # scale_factor = 1.0/(1+0.01*val**1.5)
        # scale_factor = 1.0/(1+0.09*val)
        scale_factor = 1.0 / (1 + 0.07 * val)
        # original - may be overscaling at low E
        # scale_factor = 1.0/(1+0.11*val)
        new_pot[idx] = val * scale_factor
        # tanh_scale_factor = 0.05
        # new_pot[idx] = numpy.tanh(val * tanh_scale_factor)/tanh_scale_factor

    return new_pot


def remove_empty_terms(pot):
    """Remove terms from the potential that do not have
    a value associated with them
    """
    return {k: v for k, v in pot.items() if v is not None}


def by_index(pot):
    """Build a new potential where the keys of the potential dictionary
    correspond to the indices along values of n-dimensional grids,
    rather than, possibly, the coordinate values of the grids themselves.

    Key Transformation:
    ((grid_val_i, grid_val_j, ...)_i,) -> ((i, j, ...)_i,)

    :param pot: potential along a coordinate
    :type pot: dict[tuple(float)] = float
    :rtype: dict[tuple(int)] = float
    """

    pot_keys = list(pot.keys())
    dim = dimension(pot)

    remap_dcts = []
    for i in range(dim):
        _coords = sorted(list(set(lst[i] for lst in pot_keys)))
        _idxs = list(range(len(_coords)))
        remap_dcts.append(dict(zip(_coords, _idxs)))

    new_dct = {}
    for keys in pot_keys:
        new_tup = ()
        for i, val in enumerate(keys):
            new_tup += (remap_dcts[i][val],)

        new_dct[new_tup] = pot[keys]

    return new_dct


# checks
def is_nonempty(pot):
    """Determine if the potential has any values"""
    return any(val is not None for val in pot.values())


def dimension(pot):
    """Find the dimension of the potential"""
    return len(list(pot.keys())[0])
