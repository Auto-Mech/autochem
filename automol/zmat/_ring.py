""" Level 4 Z-Matrix functions for generating ring information
"""

import numpy
import automol.graph


def ring_atoms(zma, zrxn=None):
    """ Get ring atoms.

        :param zma: Z-Matrix
        :type zma: automol.zmat object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
    """

    rngs_atoms = []
    if zrxn is None:
        rngs_atoms = automol.graph.rings_atom_keys(automol.zmat.graph(zma))
    else:
        rings = automol.reac.forming_rings_bond_keys(zrxn)
        for ring_bnds in rings:
            rng_atoms = []
            for ring_bnd in ring_bnds:
                atma, _ = ring_bnd
                if atma not in rng_atoms:
                    rng_atoms.append(atma)
            rngs_atoms.append(ring_atoms)

    return rngs_atoms


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
        samp_range_dct[key] = [value - numpy.pi/4, value + numpy.pi/4]

    return samp_range_dct


def distances_passes(samp_zma, rng_atoms, dist_value_dct):
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
            automol.zmat.distance(samp_zma, rng_atoms[i-1], rng_atoms[i])
        )
        if abs(chk_dist) > .3:
            condition = False

    return condition
