""" drivers for coordinate scans
"""

import itertools
import copy
import numpy
import automol.zmat
import automol.util


# Build the grirds ultimately used for building potentials
def grid(zma, coord_name, span, symmetry, increment, from_equilibrium=False):
    """ scan grids
    """

    # Set the space
    interval = (span / symmetry) - increment

    npoints = int(round(interval / increment, 0)) + 1
    _grid = numpy.linspace(0.0, interval, npoints)

    # Displace from the coordinates equilibrium value if desired
    if from_equilibrium:
        val_dct = automol.zmat.value_dictionary(zma)
        ini_val = val_dct[coord_name]
        _grid = automol.util.numpy_to_float(_grid)
        grid_from_equil = tuple(val + ini_val for val in _grid)

    return grid_from_equil


def points(grids):
    """ Determine the dimensions of the grid

        m = (m1, m2, m3)
        n = (n1, n2)
        p = mxn = ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
    """

    grid_points = ((i for i in range(len(grid))) for grid in grids)
    grid_points = tuple(itertools.product(*grid_points))

    return grid_points


def coords(grids):
    """ Determine the dimensions of the grid

        m = (m1, m2, m3)
        n = (n1, n2)
        p = mxn = ((m1, n1), (m1, n2), (m2, n1), (m2, n2), (m3, n1), (m3, n2))

    """

    assert len(grids) in (1, 2, 3, 4), 'Rotor must be 1-4 dimensions'

    # grid_vals = ((x for x in grid) for grid in grids)
    # grid_vals = tuple(itertools.product(*grid_vals))
    grid_vals = tuple(itertools.product(*grids))

    return grid_vals


# Manipulate potentials
def scale(pot, scale_factor):
    """ Scale the potential by scaling factor

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
        new_pot[idx] = val * scale_factor

    return new_pot


def relax_scale(pot):
    """ Scale the potential by scaling factor

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
        scale_factor = 1.0/(1+0.07*val)
        # original - may be overscaling at low E
        # scale_factor = 1.0/(1+0.11*val)
        new_pot[idx] = val * scale_factor
        # tanh_scale_factor = 0.05
        # new_pot[idx] = numpy.tanh(val * tanh_scale_factor)/tanh_scale_factor

    return new_pot


def truncate(pot, sym_num):
    """ Take a potential and limit it's terms by the symmetry number
    """

    if sym_num == 1:
        potr = copy.deepcopy(pot)
    else:
        potr = {}
        lpot = int(len(pot) / sym_num)
        for i, key in enumerate(pot.keys()):
            if i < lpot:
                potr[key] = pot[key]

    return potr


def remove_empty_terms(pot):
    """ Remove terms from the potential that do not have
        a value associated with them
    """
    return {k: v for k, v in pot.items() if v is not None}


def by_index(pot):
    """ Build a new potential where the keys of the potential dictionary
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
def valid(pot):
    """ Check if the potential is valid
    """

    is_valid = True

    dim = dimension(pot)
    for key, val in pot.items():
        if (not isinstance(key, tuple) or
                not len(key) == dim or
                not isinstance(val, float)):
            is_valid = False

    return is_valid


def is_nonempty(pot):
    """ Determine if the potential has any values
    """
    return any(val is not None for val in pot.values())


def dimension(pot):
    """ Find the dimension of the potential
    """
    return len(list(pot.keys())[0])


# I/O
def string(pot):
    """ Write a string for the potential
        maybe add a way to print coords and values?
        right now its just values
        call util.vec.string
    """

    pot_str = ''
    for val in pot.values():
        pot_str += ' {0:.6f}'.format(val)

    return pot_str
