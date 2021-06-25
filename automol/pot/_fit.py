"""
 Handle fits to potentials
"""


import numpy
from scipy.interpolate import interp1d


def fit_1d_potential(pot_dct, min_thresh=-0.0001, max_thresh=50.0):
    """ Get a physical hindered rotor potential via a series of spline fits
    """

    pot = list(pot_dct.values())

    # Check if all values are bad
    if all(numpy.isclose(-10.0, val) for val in pot):
        print('ERROR: All potential values missing')

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


# def spline_fitter(xarr, yarr):
#     """
#     """
#     x_pos = numpy.array([i for i in range(lpot)
#                          if pot[i] >= min_thresh])
#     y_pos = numpy.array([pot[i] for i in range(lpot)
#                          if pot[i] >= min_thresh])
#     pos_pot_spl = interp1d(x_pos, y_pos, kind='cubic')
#
#
# def linear_fitter(pot):
#     """ Do a one-dimensional linear fitter
#     """
#
#     neg_idxs = [i for i in range(lpot) if pot_pos_fit[i] < min_thresh]
#     clean_pot = []
#     if i in neg_idxs:
#         # Find the indices for positive vals around negative value
#         idx_0 = i - 1
#         while idx_0 in neg_idxs:
#             idx_0 = idx_0 - 1
#         for j in range(i, lpot):
#             if pot_pos_fit[j] >= min_thresh:
#                 idx_1 = j
#                 break
#         pot = _linear_fitter(pot)
#
#
# def _linear_fit(pot, idx_0, idx_1):
#     """ Linear fitter
#     """
#     interp_val = (
#         pot[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
#         pot[idx_1] * ((i-idx_0)/(idx_1-idx_0))
#     )
#
#
# def _re_set_high_values(pot, max_thresh=600., min_thresh=-0.0001):
#     """ Rebuild the potential
#     """
#
#     idx_success = []
#     pot_success = []
#     for idx in range(lpot):
#         if pot[idx] < 600. and pot[idx] > min_thresh:
#             idx_success.append(idx)
#             if pot[idx] < max_thresh:
#                 pot_success.append(pot[idx])
#             else:
#                 pot_success.append(max_thresh)
#
#     return pot
