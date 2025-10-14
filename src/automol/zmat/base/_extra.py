""" Level 2 functions for sampling and related actions
"""

import itertools

import numpy

from ._core import set_values_by_name, value_dictionary
from ...vmat import coordinates

def samples(zma, nsamp, range_dct):
    """randomly sample over torsional dihedrals"""
    _names = tuple(range_dct.keys())
    ranges = tuple(range_dct.values())
    vals_lst = _sample_over_ranges(ranges, nsamp)

    zmas = tuple(
        set_values_by_name(zma, dict(zip(_names, vals)), degree=False)
        for vals in vals_lst
    )

    return zmas


def constraint_dict(zma, const_names, var_names=()):
    """Build a dictionary of constraints

    has to round the values for the filesystem
    """

    # Get the list names sorted for dictionary
    rnames = (name for name in const_names if "R" in name)
    anames = (name for name in const_names if "A" in name)
    dnames = (name for name in const_names if "D" in name)
    rnames = tuple(sorted(rnames, key=lambda x: int(x.split("R")[1])))
    anames = tuple(sorted(anames, key=lambda x: int(x.split("A")[1])))
    dnames = tuple(sorted(dnames, key=lambda x: int(x.split("D")[1])))
    constraint_names = rnames + anames + dnames

    # Remove the scan coordinates so they are not placed in the dict
    constraint_names = tuple(name for name in constraint_names if name not in var_names)

    # Build dictionary
    if constraint_names:
        zma_vals = value_dictionary(zma)
        zma_coords = coordinates(zma)
        assert set(constraint_names) <= set(zma_coords.keys()), (
            "Attempting to constrain coordinates not in zma:"
            f"\n{constraint_names}\n{zma_coords}"
        )
        _dct = dict(
            zip(
                constraint_names,
                (round(zma_vals[name], 2) for name in constraint_names),
            )
        )
    else:
        _dct = None

    return _dct


def set_constraint_names(zma, tors_names, tors_model):
    """Determine the names of constraints along a torsion scan"""

    constraint_models = ("1dhrf", "1dhrfa", "tau-1dhrf", "tau-1dhrfa")

    const_names = tuple()
    if tors_names and tors_model in constraint_models:
        if tors_model in ("1dhrf", "tau-1dhrf"):
            const_names = tuple(itertools.chain(*tors_names))
        elif tors_model in ("1dhrfa", "tau-1dhrfa"):
            coords = list(coordinates(zma))
            const_names = tuple(coord for coord in coords)

    return const_names


# helpers
def _sample_over_ranges(rngs, nsamp):
    """randomly sample over several ranges"""
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))
