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
    """ sampling ranges for torsional dihedrals
    """
    # sym_nums = torsional_symmetry_numbers(zma, tors_names,
    # frm_bnd_key=None, brk_bnd_key=None)
    # return tuple((0, 2*numpy.pi/sym_num) for sym_num in sym_nums)
    # originally restricted range by sym_num.
    # But after all it appears that using the
    # full range is best.
    # after all it appears that using the full sampling range is most effective
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
            'Attempting to constrain coordinates not in zma:\n{}\n{}'.format(
                constraint_names, zma_coords)
        )
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

    const_names = tuple()
    if tors_names and tors_model in ('1dhrf', '1dhrfa'):
        if tors_model == '1dhrf':
            const_names = tuple(
                itertools.chain(*tors_names))
        elif tors_model == '1dhrfa':
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


def get_babs1(zma, dist_name):
    """ get name of torsional coordinate associated with babs1 pre-reformatting
    """
    idxs = coord_idxs(zma, dist_name)
    idx = max(idxs)
    babs1 = 'D{:g}'.format(idx)
    return babs1


def get_babs2(zma, dist_name):
    """ get name of torsional coordinate associated with babs2 pre-reformatting
    """
    idxs = coord_idxs(zma, dist_name)
    idx = max(idxs)
    babs2 = 'D{:g}'.format(idx+1)
    return babs2


# # helpers
def _sample_over_ranges(rngs, nsamp):
    """ randomly sample over several ranges
    """
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))
