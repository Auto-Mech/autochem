"""
 Handle fits to potentials
"""


import numpy
# from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline


def _initialize_pot(pot_dct):
    grid_vals, pot = zip(*pot_dct.items())
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


def _cap_max_thresh_drop_neg(pot, lpot, min_thresh, max_thresh):
    x_idxs = []
    y_pots = []
    for idx in range(lpot):
        if pot[idx] < 600. and pot[idx] > min_thresh:
            x_idxs.append(idx)
            if pot[idx] < max_thresh:
                y_pots.append(pot[idx])
            else:
                y_pots.append(max_thresh)
    return numpy.array(x_idxs), numpy.array(y_pots)


def _periodic_spline(x_vals, y_vals, orig_len):
    """ perform a cubic periodic spline on an x and y array
    """
    spl_fit = []
    if len(y_vals) > 3:
        spl_fun = CubicSpline(x_vals, y_vals, bc_type='periodic')
        for idx in range(orig_len):
            spl_fit.append(float(spl_fun(idx)))
    return spl_fit


def _spline_cap_max_thresh(pot, lpot, max_thresh, min_thresh):
    """replace any negative potential values cubic spline fit values
    """
    x_idxs, y_pots = _cap_max_thresh_drop_neg(
        pot, lpot, min_thresh, max_thresh)
    spl_fit = _periodic_spline(x_idxs, y_pots, lpot)
    return pot if spl_fit else spl_fit


def _spline_no_alteration(pot, lpot, min_thresh):
    """Do a spline fit for negatives without changing anything else
    """
    x_pos = numpy.array([i for i in range(lpot)
                         if pot[i] >= min_thresh])
    y_pos = numpy.array([pot[i] for i in range(lpot)
                         if pot[i] >= min_thresh])
    spl_fit = _periodic_spline(x_pos, y_pos, lpot)
    return pot if spl_fit else spl_fit


def fit_1d_potential(pot_dct, min_thresh=-0.0001, max_thresh=50.0):
    """ Get a physical hindered rotor potential via a series of spline fits
    """
    # unpack potential dictionary
    _, pot, lpot = _initialize_pot(pot_dct)
    pot = _round_near_zero(pot, max_thresh, min_thresh)
    # Build a new potential list using a spline fit of the HR potential
    pot = _spline_cap_max_thresh(pot, lpot, max_thresh, min_thresh)

    # Do second spline fit of only positive values if any negative values found
    if any(val < min_thresh for val in pot):
        print('Still found negative potential values after first spline')
        print('Potential after first spline:', pot)
        pot = _spline_no_alteration(pot, lpot, min_thresh)
        print('Potential after second spline:', pot)
        # Perform second check to see if negative potentials have been fixed
        if any(val < min_thresh for val in pot):
            print('Still found negative potential values after second spline')
            print('Replace with linear interpolation of positive values')
            neg_idxs = [i for i in range(lpot) if pot[i] < min_thresh]
            clean_pot = []
            for i in range(lpot):
                if i in neg_idxs:
                    # Find the indices for positive vals around negative value
                    idx_0 = i - 1
                    idx_1 = None
                    while idx_0 in neg_idxs:
                        idx_0 = idx_0 - 1
                    for j in range(i, lpot):
                        if pot[j] >= min_thresh:
                            idx_1 = j
                            break
                    if idx_1 is None:
                        for j in range(0, idx_0):
                            idx_1 = j + lpot
                    # Get a new value for this point on the potential by
                    # doing a linear interp of positives
                    interp_val = (
                        pot[idx_0] * (1.0-((i-idx_0)/(idx_1-idx_0))) +
                        (pot[idx_1 if idx_1 < lpot else idx_1 - lpot]
                         * ((i-idx_0)/(idx_1-idx_0)))
                    )
                    clean_pot.append(interp_val)
                else:
                    clean_pot.append(pot[i])
            final_potential = clean_pot.copy()
            print('Potential after linear interp:', final_potential)

        else:
            final_potential = pot.copy()

        # check for quadratic behavior around minimum
    else:
        final_potential = pot.copy()
    # if lpot > 7:
    #     pos = final_potential[1]
    #     neg = final_potential[-2]
    #     if abs(pos - neg)/(pos + neg) > 0.3:
    #         print('Refitting potentials near the minimum to a quadratic')
    #         pos_1 = final_potential[2]
    #         neg_1 = final_potential[-3]
    #         quad_pots = [neg_1, neg, 0, pos, pos_1]
    #         quad_idxs = [0, 1, 2, 3, 4]
    #         quad_pot_coeffs = numpy.polyfit(quad_idxs, quad_pots, 2)
    #         print(quad_pots)
    #         print(quad_pot_coeffs)
    #         final_potential[-2] = (quad_pot_coeffs[2] +
    #                                quad_pot_coeffs[1] * 1 +
    #                                quad_pot_coeffs[0] * 1**2)
    #         final_potential[1] = (quad_pot_coeffs[2] +
    #                               quad_pot_coeffs[1] * 3 +
    #                               quad_pot_coeffs[0] * 3**2)
    #         print('Potential to fit quadratic', final_potential)

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
