""" Level 4 Z-Matrix functions for generating ring information
"""

import math

from automol.graph import base as graph_base
from automol.zmat._conv import (
    distance,
    graph,
)
from automol.zmat.base import (
    coordinates,
    dihedral_angle_names,
    value_dictionary,
)


# Get information for all rings at once
def all_rings_atoms(zma, tsg=None):
    """Get ring atoms.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    if tsg is None:
        rings_atoms = graph_base.rings_atom_keys(graph(zma))
    else:
        rings_atoms = []
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

            # Add to overall list
            rings_atoms.append(ring_atoms)

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
    return tuple(ring_dihedrals(zma, ring_atoms) for ring_atoms in rings_atoms)


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


def ring_distances_reasonable(zma, rng_atoms, dist_value_dct):
    """Are the distances between ring atoms reasonable?

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    condition = True
    for i, rng_atom in enumerate(rng_atoms):
        chk_dist = dist_value_dct[i] - distance(zma, rng_atoms[i - 1], rng_atom)
        if abs(chk_dist) > 0.3:
            condition = False

    return condition


def ring_dihedrals(zma, rng_atoms):
    """Get ring dihedral names and their angle values

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


def ring_samp_ranges(zma, rng_atoms):
    """Set sampling range for ring dihedrals.

    :param zma: Z-Matrix
    :type zma: automol.zmat object
    :param rng_atoms: idxs for atoms inside rings
    :type rng_atoms: list
    """

    samp_range_dct = {}
    ring_value_dct = ring_dihedrals(zma, rng_atoms)
    for key, value in ring_value_dct.items():
        samp_range_dct[key] = [value - math.pi / 4, value + math.pi / 4]

    return samp_range_dct
