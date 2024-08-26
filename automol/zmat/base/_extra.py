""" Level 2 functions for sampling and related actions
"""

import itertools

import numpy

from ._core import set_values_by_name, value_dictionary
from ._core import value, set_key_matrix
from ...geom import dihedral_angle
from ...zmat import coordinates as zcoords
from ...vmat import coordinates, key_matrix

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


def samples_avg_dih(zma, geo, tors_dcts, average_dih,ring_tors_dct,dih_remover):
    """randomly sample over torsional dihedrals"""
    vals_lst = []
    iter_combos = []
    _names, avg_dih = [], []

    all_sampled_torsions = []
    for _,samp_range_dct in tors_dcts:
        all_sampled_torsions.extend(list(samp_range_dct.keys()))

    for key_dct,samp_range_dct in tors_dcts:

        repeats = len(samp_range_dct.keys())
        ring_atoms = [int(idx)-1 for idx in key_dct.split('-')]

        if repeats != len(ring_atoms):
            keymat = [list(item) for item in(key_matrix(zma))]
            # Find all DHs of ring atoms from 4th atom
            # is DH in tors-dcts?
            changed_dh = []

            for at in ring_atoms[3:]:
                dih=f"D{at}"
                if dih not in [b for a,b in dih_remover if a != key_dct]:
                    if dih not in all_sampled_torsions:
                        changed_dh.append(dih)
                        idx = ring_atoms.index(at)
                        keymat[at] = [ring_atoms[idx-1], \
                                      ring_atoms[idx-2], ring_atoms[idx-3]]

                        ring_tors_dct[key_dct].update({dih:(0.,0.)})

            zma = set_key_matrix(zma, keymat)

            # rebuild the zmatrix with new value of that dih
            key_coord_dct = zcoords(zma)
            new_key_dct = {}
            for name,coos in key_coord_dct.items():
                atm_idxs = coos[0]
                if len(atm_idxs) == 2:
                    new_key_dct[name] = value(zma, name, angstrom=True)
                elif len(atm_idxs) == 3:
                    new_key_dct[name] = value(zma, name, degree=True)
                elif len(atm_idxs) == 4:
                    new_key_dct[name] = value(zma, name, degree=True)
                    if name in changed_dh: # compute DH with previous three ring atoms
                        new_key_dct[name] = dihedral_angle(geo, *atm_idxs,degree=True)
            zma = set_values_by_name(zma, new_key_dct)


    for key_dct,samp_range_dct in tors_dcts:

        repeats = len(samp_range_dct.keys())

        _names.extend(list(samp_range_dct.keys()))
        iter_combos.append(itertools.product([1,0,-1],repeat=repeats))
        avg_dih.append(average_dih[key_dct])

    samp_mat = itertools.product(*iter_combos)

    for samp in samp_mat:
        vals = []
        for i,dih_value in enumerate(avg_dih):
            vals.extend( [val * dih_value for val in samp[i]] )
        vals_lst.append(tuple(vals))

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
