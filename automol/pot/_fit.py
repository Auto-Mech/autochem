"""
 Handle fits to potentials
"""


import numpy


def setup_1d_potential(pot_dct, min_thresh=-0.0001, max_thresh=50.0):
    """ Get a physical hindered rotor potential via a series of spline fits
        but just drop negative potentials instead of fitting them
    """
    # unpack potential dictionary
    _, pot, lpot = _initialize_pot(pot_dct)
    pot = _round_near_zero(pot, max_thresh, min_thresh)

    x_idxs, y_pots = _cap_max_thresh_drop_neg(
        pot, lpot, min_thresh, max_thresh)

    new_pot = {}
    for idx, x_val in enumerate(x_idxs[:-1]):
        new_pot[(x_val,)] = y_pots[idx]

    return new_pot


def _initialize_pot(pot_dct):
    grid_vals, pot = pot_dct.keys(), pot_dct.values()
    pot = list(pot)
    # Check if all values are bad
    if all(numpy.isclose(-10.0, val) for val in pot):
        print('ERROR: All potential values missing')

    # Initialize a variable for the size of the potential
    lpot = len(pot)+1
    pot.append(0.0)

    return grid_vals, pot, lpot


def _round_near_zero(pot, max_thresh, min_thresh):
    """ Caps the potentials over a given threshold
    """
    # Print warning messages
    # Build a potential list from only successful calculations
    # First replace high potential values with max_thresh
    print_pot = False
    if any(val > max_thresh for val in pot):
        print_pot = True
        print(f'Warning: Found pot val of {max(pot):.2f}',
              ' which is larger than',
              'the typical maximum for a torsional potential')

    # reset any negative values for the first grid point to 0.
    if pot[0] < 0.:
        print('ERROR: The first potential value should be 0.')
        pot[0] = 0.0
    elif pot[0] < .01:
        pot[0] = 0.0
    if any(val < min_thresh for val in pot):
        print_pot = True
        print(f'Warning: Found pot val of {min(pot):.2f}',
              ' which is below',
              f'{min_thresh} kcal. Refit w/ positives')

    if print_pot:
        print('Potential before spline:', pot)
    return pot


def _cap_max_thresh_drop_neg(pot, lpot, min_thresh, max_thresh,
                             zero_thresh=-5):
    x_idxs = []
    y_pots = []
    for idx in range(lpot):
        if pot[idx] < 600. and pot[idx] > min_thresh:
            x_idxs.append(idx)
            if pot[idx] < max_thresh:
                y_pots.append(pot[idx])
            else:
                y_pots.append(max_thresh)
        elif pot[idx] < 600. and pot[idx] > zero_thresh:
            x_idxs.append(idx)
            y_pots.append(0)
    return numpy.array(x_idxs), numpy.array(y_pots)
