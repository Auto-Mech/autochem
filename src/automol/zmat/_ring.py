""" Level 4 Z-Matrix functions for generating ring information
"""

import math
import itertools

from ..graph import base as graph_base
from ._conv import (
    distance,
    graph,
)
from .base import (
    coordinates,
    dihedral_angle_names,
    value_dictionary,
    set_values_by_name,
    value,
    set_key_matrix,
)

from ..geom import dihedral_angle
from ..vmat import key_matrix


# Get information for all rings at once
def all_rings_atoms(zma, tsg=None):
    """Get ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    rings_atoms = graph_base.rings_atom_keys(graph(zma, dummy=True))
    if tsg is not None:
        for ring in graph_base.ts.forming_rings_bond_keys(tsg):
            # Determine number of atoms in the ring
            all_atoms = set()
            for ring_bnd in ring:
                all_atoms = all_atoms | ring_bnd
            natoms = len(all_atoms)

            # intialize list with indices from first bond
            atma, atmb = list(ring)[0]
            ring_atoms = [atma, atmb]

            # Iteratively add to ring idx list by finding with bnd has
            # the idx at end of current list to maintain connectivity
            while len(ring_atoms) != natoms:
                for ring_bnd in ring:
                    atma, atmb = ring_bnd
                    if atma == ring_atoms[-1] and atmb not in ring_atoms:
                        ring_atoms.append(atmb)
                    elif atmb == ring_atoms[-1] and atma not in ring_atoms:
                        ring_atoms.append(atma)

            # Check that ring is not already present in rings_atoms
            do_not_add_ring = 0
            for rng in rings_atoms:
                if not set(rng).difference(set(ring_atoms)):
                    do_not_add_ring = 1
            if do_not_add_ring:
                continue

            # Add to overall list
            rings_atoms = list(rings_atoms)
            # Added sort as I expect that connectivity is defined in "usual" way
            ring_atoms = sorted(ring_atoms)
            # Connectivity should still be preserved,
            # If it is not COME BACK HERE AND FIX
            # I think it is not ALWAYS preserved (see fused rings!!)
            #  so maybe I should remove the sort

            rings_atoms.append(tuple(ring_atoms))
            rings_atoms = frozenset(rings_atoms)

    return rings_atoms


def all_rings_distances(zma, rings_atoms):
    """For every ring present in the system. determine the
    distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    return tuple(ring_distances(zma, ring_atoms) for ring_atoms in rings_atoms)


def all_rings_distances_reasonable(zma, rings_atoms):
    """For every ring present in the system. determine the
    distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    # HOW CAN IT EVER BE FALSE IF I CREATE VALUE DCT AND USE IT WITH SAME ZMA
    # Currently only used in tests
    condition = True
    for ring_atoms in rings_atoms:
        dist_val_dct = ring_distances(zma, ring_atoms)
        condition = ring_distances_reasonable(zma, ring_atoms, dist_val_dct)

    return condition


def all_rings_dihedrals(zma, rings_atoms):
    """Get ring dihedral names and their angle values

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """
    #return tuple(ring_dihedrals(zma, ring_atoms) for ring_atoms in rings_atoms)
    return tuple(complete_ring_dihedrals(zma, ring_atoms) for ring_atoms in rings_atoms)


def all_rings_dct(zma, rings_atoms):
    """Build a dictionary which relates the indices of the atoms
    of various rings to their dihedrals and sampling ranges.

    {rng_idx1-rng_idx2-rng_idx3: {Dn: [min, max], Dn2: [min, max]}}
    """

    ring_dct = {}
    for ring_atoms in rings_atoms:
        dct_label = "-".join(str(atm + 1) for atm in ring_atoms)
        ring_dct[dct_label] = ring_samp_ranges(zma, ring_atoms)

    return ring_dct


# Functions for a single ring
def ring_distances(zma, rng_atoms):
    """Return the distances between each pair of ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    dist_value_dct = {}
    for i, _ in enumerate(rng_atoms):
        dist_value_dct[i] = distance(zma, rng_atoms[i - 1], rng_atoms[i])

    return dist_value_dct


def ring_distances_reasonable(zma, rng_atoms, dist_value_dct, thresh=0.3):
    """Are the distances between ring atoms reasonable?

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    condition = True
    for i, rng_atom in enumerate(rng_atoms):
        chk_dist = dist_value_dct[i] - distance(zma, rng_atoms[i - 1], rng_atom)
        if abs(chk_dist) > thresh:
            condition = False

    return condition


def ring_dihedrals(zma, rng_atoms):
    """Get N-3 ring dihedral names and their angle values

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    coos = coordinates(zma)
    da_names = dihedral_angle_names(zma)
    val_dct = value_dictionary(zma)

    ring_value_dct = {}
    for da_name in da_names:
        da_idxs = list(coos[da_name])[0]
        if len(list(set(da_idxs) & set(rng_atoms))) == 4:
            ring_value_dct[da_name] = val_dct[da_name]

    return ring_value_dct

def complete_ring_dihedrals(zma, rng_atoms):
    """Get ring dihedral names and their angle values
    of all the atoms involved in a ring, not only N-3 dihs

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    coos = coordinates(zma)
    da_names = dihedral_angle_names(zma)
    val_dct = value_dictionary(zma)

    ring_value_dct = {}
    for da_name in da_names:
        da_idxs = list(coos[da_name])[0]
        if da_idxs[0] in rng_atoms:
            ring_value_dct[da_name] = val_dct[da_name]

    return ring_value_dct


def ring_samp_ranges(zma, rng_atoms):
    """Set sampling range for ring dihedrals.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    samp_range_dct = {}
    ring_value_dct = ring_dihedrals(zma, rng_atoms)
    for key, val in ring_value_dct.items():
        samp_range_dct[key] = [val - math.pi / 4, val + math.pi / 4]

    return samp_range_dct


def samples_avg_dih(zma, geo, tors_dcts, average_dih,ring_tors_dct,dih_remover):
    """Samples N-3 dihs of a ring giving either the average value of dihs
    of the ring, 0 or the opposite of the average.
    """
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
                if dih in dihedral_angle_names(zma):
                    if dih not in [b for a,b in dih_remover if a != key_dct]:
                        if dih not in all_sampled_torsions:
                            changed_dh.append(dih)
                            idx = ring_atoms.index(at)
                            keymat[at] = [ring_atoms[idx-1], \
                                            ring_atoms[idx-2], ring_atoms[idx-3]]

                            ring_tors_dct[key_dct].update({dih:(0.,0.)})

            zma = set_key_matrix(zma, keymat)

            # rebuild the zmatrix with new value of that dih
            key_coord_dct = coordinates(zma) #zcoords(zma)
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
