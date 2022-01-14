""" Level 2 functions for sampling and related actions
"""

import itertools
import numpy
from automol.vmat import coordinates
from automol.zmat.base._core import value_dictionary
from automol.zmat.base._core import set_values_by_name


def samples(zma, nsamp, range_dct):
    """ randomly sample over torsional dihedrals
    """
    _names = tuple(range_dct.keys())
    ranges = tuple(range_dct.values())
    vals_lst = _sample_over_ranges(ranges, nsamp)

    zmas = tuple(set_values_by_name(zma, dict(zip(_names, vals)), degree=False)
                 for vals in vals_lst)

    return zmas


def torsional_sampling_ranges(tors_names):
    """ Generate the min and max for a range of values that can be used
        to sample the torsional angle, relative to the equibrium structure.

        Function originally restricted the range by the torsional symmetry
        number; however, it seems most effective to sample the full space
        of values from  0 to 2*pi.
    """
    sym_nums = 1.
    return tuple((0, 2*numpy.pi/sym_nums) for tors_name in tors_names)


def constraint_dct(zma, const_names, var_names=()):
    """ Build a dictionary of constraints

        has to round the values for the filesystem
    """

    # Get the list names sorted for dictionary
    rnames = (name for name in const_names if 'R' in name)
    anames = (name for name in const_names if 'A' in name)
    dnames = (name for name in const_names if 'D' in name)
    rnames = tuple(sorted(rnames, key=lambda x: int(x.split('R')[1])))
    anames = tuple(sorted(anames, key=lambda x: int(x.split('A')[1])))
    dnames = tuple(sorted(dnames, key=lambda x: int(x.split('D')[1])))
    constraint_names = rnames + anames + dnames

    # Remove the scan coordinates so they are not placed in the dict
    constraint_names = tuple(name for name in constraint_names
                             if name not in var_names)

    # Build dictionary
    if constraint_names:
        zma_vals = value_dictionary(zma)
        zma_coords = coordinates(zma)
        assert set(constraint_names) <= set(zma_coords.keys()), (
            'Attempting to constrain coordinates not in zma:'
            f'\n{constraint_names}\n{zma_coords}')
        _dct = dict(zip(
            constraint_names,
            (round(zma_vals[name], 2) for name in constraint_names)
        ))
    else:
        _dct = None

    return _dct


def set_constraint_names(zma, tors_names, tors_model):
    """ Determine the names of constraints along a torsion scan
    """

    constraint_models = ('1dhrf', '1dhrfa', 'tau-1dhrf', 'tau-1dhrfa')

    const_names = tuple()
    if tors_names and tors_model in constraint_models:
        if tors_model in ('1dhrf', 'tau-1dhrf'):
            const_names = tuple(
                itertools.chain(*tors_names))
        elif tors_model in ('1dhrfa', 'tau-1dhrfa'):
            coords = list(coordinates(zma))
            const_names = tuple(coord for coord in coords)

    return const_names


def coord_idxs(zma, key):
    """give a bond length key, return the indices of involved bonded atoms
    """
    coords = coordinates(zma)
    idxs = coords.get(key, [None])
    return idxs[0]


def bond_key_from_idxs(zma, idxs):
    """given indices of involved bonded atoms, return bond name
    """
    idxs = list(idxs)
    idxs.sort(reverse=True)
    idxs = tuple(idxs)
    bond_key = None
    coords = coordinates(zma)
    for key in coords:
        for coord in coords.get(key, [None]):
            if idxs == coord:
                bond_key = key
    return bond_key


# helpers
def _sample_over_ranges(rngs, nsamp):
    """ randomly sample over several ranges
    """
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))
