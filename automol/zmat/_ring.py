""" Level 4 Z-Matrix functions for generating ring information
"""

import math
import automol.graph


# Get information for all rings at once
def all_rings_atoms(zma, zrxn=None):
    """ Get ring atoms.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    if zrxn is None:
        rings_atoms = automol.graph.rings_atom_keys(automol.zmat.graph(zma))
    else:
        rings = automol.reac.forming_rings_bond_keys(zrxn)

        rings_atoms = []
        for ring_bnds in rings:
            ring_atoms = []
            for ring_bnd in ring_bnds:
                atma, _ = ring_bnd
                if atma not in ring_atoms:
                    ring_atoms.append(atma)
            rings_atoms.append(ring_atoms)

    return rings_atoms


def all_rings_distances(zma, rings_atoms):
    """ For every ring present in the system. determine the
        distances between each pair of ring atoms.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """
    return tuple(ring_distances(zma, ring_atoms)
                 for ring_atoms in rings_atoms)


def all_rings_distances_reasonable(zma, rings_atoms):
    """ For every ring present in the system. determine the
        distances between each pair of ring atoms.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    condition = True
    for ring_atoms in rings_atoms:
        dist_val_dct = ring_distances(
            zma, ring_atoms)
        condition = ring_distances_reasonable(
            zma, ring_atoms, dist_val_dct)

    return condition


def all_rings_dihedrals(zma, rings_atoms):
    """ Get ring dihedral names and their angle values

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """
    return tuple(ring_dihedrals(zma, ring_atoms)
                 for ring_atoms in rings_atoms)


def all_rings_dct(zma, rings_atoms):
    """ Build a dictionary which relates the indices of the atoms
        of various rings to their dihedrals and sampling ranges.

        {rng_idx1-rng_idx2-rng_idx3: {Dn: [min, max], Dn2: [min, max]}}
    """

    ring_dct = {}
    for ring_atoms in rings_atoms:
        dct_label = '-'.join(str(atm+1) for atm in ring_atoms)
        ring_dct[dct_label] = ring_samp_ranges(zma, ring_atoms)

    return ring_dct


# Functions for a single ring
def ring_distances(zma, rng_atoms):
    """ Return the distances between each pair of ring atoms.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    dist_value_dct = {}
    for i, _ in enumerate(rng_atoms):
        dist_value_dct[i] = automol.zmat.distance(
            zma, rng_atoms[i-1], rng_atoms[i])

    return dist_value_dct


def ring_distances_reasonable(zma, rng_atoms, dist_value_dct):
    """ Are the distances between ring atoms reasonable?

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    condition = True
    for i, _ in enumerate(rng_atoms):
        chk_dist = (
            dist_value_dct[i] -
            automol.zmat.distance(zma, rng_atoms[i-1], rng_atoms[i])
        )
        if abs(chk_dist) > .3:
            condition = False

    return condition


def ring_dihedrals(zma, rng_atoms):
    """ Get ring dihedral names and their angle values

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    coos = automol.zmat.coordinates(zma)
    da_names = automol.zmat.dihedral_angle_names(zma)
    val_dct = automol.zmat.value_dictionary(zma)

    ring_value_dct = {}
    for da_name in da_names:
        da_idxs = list(coos[da_name])[0]
        if len(list(set(da_idxs) & set(rng_atoms))) == 4:
            ring_value_dct[da_name] = val_dct[da_name]

    return ring_value_dct


def ring_samp_ranges(zma, rng_atoms):
    """ Set sampling range for ring dihedrals.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    samp_range_dct = {}
    ring_value_dct = ring_dihedrals(zma, rng_atoms)
    for key, value in ring_value_dct.items():
        samp_range_dct[key] = [value - math.pi/4, value + math.pi/4]

    return samp_range_dct
