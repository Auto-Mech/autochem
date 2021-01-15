""" drivers for coordinate scans
"""

import itertools
import copy
import numpy
from scipy.interpolate import interp1d
import automol


# Build the potential objects
def points(grids):
    """ Determine the dimensions of the grid

        m = (m1, m2, m3)
        n = (n1, n2)
        p = mxn = ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
    """
    assert len(grids) in (1, 2, 3, 4), 'Rotor must be 1-4 dimensions'

    grid_points = ((i for i in range(len(grid)))
                   for grid in grids)
    grid_points = tuple(itertools.product(*grid_points))

    return grid_points


def coords(grids):
    """ Determine the dimensions of the grid

        m = (m1, m2, m3)
        n = (n1, n2)
        p = mxn = ((m1, n1), (m1, n2), (m2, n1), (m2, n2), (m3, n1), (m3, n2))

    """
    assert len(grids) in (1, 2, 3, 4), 'Rotor must be 1-4 dimensions'

    grid_vals = ((x for x in grid)
                 for grid in grids)
    grid_vals = tuple(itertools.product(*grid_vals))

    return grid_vals


# Manipulate potentiasls
def scale(pot, scale_coeff, num_tors):
    """ Scale the potential by scaling factor

        :param pot: potential along a coordinate
        :type pot: dict[tuple(float)] = float
        :param scale_coeff: initial scaling coeffcient
        :type scael_coeff: float
        :param num_tors: number of torsions used in scaling
        :type num_tors: int
        :rtype:
    """

    scale_factor = scale_coeff**(2.0/num_tors)

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
        _idxs = [i for i in range(len(_coords))]
        remap_dcts.append(dict(zip(_coords, _idxs)))

    new_dct = {}
    for keys in pot_keys:
        new_tup = ()
        for i, val in enumerate(keys):
            new_tup += (remap_dcts[i][val],)

        new_dct[new_tup] = pot[keys]

    return new_dct


def hrpot_spline_fitter(pot_dct, min_thresh=-0.0001, max_thresh=50.0):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    pot = list(pot_dct.values())

    # Initialize a variable for the size of the potential
    lpot = len(pot)+1
    pot.append(0.0)

    # Print warning messages
    print_pot = False
    if any(val > max_thresh for val in pot):
        print_pot = True
        max_pot = max(pot)
        print('Warning: Found pot val of {0:.2f}'.format(max_pot),
              ' which is larger than',
              'the typical maximum for a torsional potential')
    # reset any negative values for the first grid point to 0.
    if pot[0] < 0.:
        print('ERROR: The first potential value should be 0.')
        pot[0] = 0.
    if any(val < min_thresh for val in pot):
        print_pot = True
        min_pot = min(pot)
        print('Warning: Found pot val of {0:.2f}'.format(min_pot),
              ' which is below',
              '{0} kcal. Refit w/ positives'.format(min_thresh))

    if print_pot:
        print('Potential before spline:', pot)

    # Build a potential list from only successful calculations
    # First replace high potential values with max_thresh
    # Then replace any negative potential values cubic spline fit values
    idx_success = []
    pot_success = []
    for idx in range(lpot):
        if pot[idx] < 600. and pot[idx] > min_thresh:
            idx_success.append(idx)
            if pot[idx] < max_thresh:
                pot_success.append(pot[idx])
            else:
                pot_success.append(max_thresh)

    if len(pot_success) > 3:
    # Build a new potential list using a spline fit of the HR potential
        pot_spl = interp1d(
            numpy.array(idx_success), numpy.array(pot_success), kind='cubic')
        for idx in range(lpot):
            pot[idx] = float(pot_spl(idx))

    # Do second spline fit of only positive values if any negative values found
    if any(val < min_thresh for val in pot):
        print('Still found negative potential values after first spline')
        print('Potential after spline:', pot)
        if len(pot_success) > 3:
            x_pos = numpy.array([i for i in range(lpot)
                                 if pot[i] >= min_thresh])
            y_pos = numpy.array([pot[i] for i in range(lpot)
                                 if pot[i] >= min_thresh])
            pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
            pot_pos_fit = []
            for idx in range(lpot):
                pot_pos_fit.append(pos_pot_spl(idx))
        else:
            pot_pos_fit = []
            for idx in range(lpot):
                pot_pos_fit.append(pot[idx])

        print('Potential after spline:', pot_pos_fit)
        # Perform second check to see if negative potentials have been fixed
        if any(val < min_thresh for val in pot_pos_fit):
            print('Still found negative potential values after second spline')
            print('Replace with linear interpolation of positive values')
            neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < min_thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot_pos_fit[j] >= min_thresh:
                            idx_1 = j
                            break
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot_pos_fit[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
                        pot_pos_fit[idx_1] * ((i-idx_0)/(idx_1-idx_0))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot_pos_fit[i])
            final_potential = clean_pot.copy()

        else:
            final_potential = pot_pos_fit.copy()

    else:
        final_potential = pot.copy()

    final_potential = final_potential[:-1]

    fin_dct = {}
    for i, val in enumerate(final_potential):
        val_fin = min(val, max_thresh)
        fin_dct[(i,)] = val_fin

    return fin_dct


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


# I/O
def string(tors_pots):
    """ Check hr pot to see if a new mimnimum is needed
    """

    print('\nHR potentials...')
    for name in tors_pots:

        print('- Rotor {}'.format(name))
        pot_str = ''
        for pot in tors_pots[name].values():
            pot_str += ' {0:.2f}'.format(pot)

        print('- Pot:{}'.format(pot_str))


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
                print(automol.zmatrix.string(zma))

    return new_min_zma
