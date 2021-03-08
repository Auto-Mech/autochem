""" drivers for coordinate scans
"""

import itertools
import copy
import numpy
import automol


# Build the grirds ultimately used for building potentials
def grid(zma, coord_name, span, symmetry, increment, from_equilibrium=False):
    """ scan grids
    """

    # Set the space
    interval = (span / symmetry) - increment

    npoints = int(round(interval / increment, 0)) + 1
    print('npoints', npoints)
    _grid = numpy.linspace(0.0, interval, npoints)
    print('grid', _grid)

    # Displace from the coordinates equilibrium value if desired
    if from_equilibrium:
        val_dct = automol.zmat.value_dictionary(zma)
        ini_val = val_dct[coord_name]
        print('init_val', ini_val)
        grid_from_equil = tuple(val.item() + ini_val for val in _grid)

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


# Manipulate potentiasls
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


def by_index(pot):
    """ Build a new potential where coordinates change by index

        :param pot: potential along a coordinate
        :type pot: dict[tuple(float)] = float
        :param coords: coordinates of potential

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


def dimension(pot):
    """ Find the dimension of the potential
    """
    return len(list(pot.keys())[0])


def check_hr_pot(tors_pots, tors_zmas, tors_paths, emax=-0.5, emin=-10.0):
    """ Check hr pot to see if a new mimnimum is needed
    """

    new_min_zma = None

    print('\nAssessing the HR potential...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pots = tors_pots[name].values()
        zmas = tors_zmas[name].values()
        paths = tors_paths[name].values()
        for pot, zma, path in zip(pots, zmas, paths):
            if emin < pot < emax:
                new_min_zma = zma
                emin = pot
                print(' - New minimmum energy ZMA found for torsion')
                print(' - Ene = {}'.format(pot))
                print(' - Found at path: {}'.format(path))
                print(automol.zmat.string(zma))

    return new_min_zma


# I/O
def string(pot):
    """ Write a string for the potential
        maybe add a way to print coords and values?
        right now its just values
        call util.vec.string
    """

    pot_str = ''
    for val in pot.values():
        pot_str += ' {0:.6f}'.format(pot)

    return pot_str
