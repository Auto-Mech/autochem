""" functions for working with torsion axis (rotational bond)s and groups
"""
import automol.convert.zmat
from automol.graph._graph import bond_keys
from automol.zmat._zmat import key_matrix
from automol.zmat._zmat import name_matrix
from automol.zmat._zmat import string


def distance_coordinate_name(zma, key1, key2):
    """ get the name of a distance coordinate for a given bond

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key in the torsion axis (rotational bond)
        :type key1: int
        :param key2: the second key in the torsion axis (rotational bond)
        :type key2: int
        :rtype: str
    """

    key1, key2 = sorted([key1, key2])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert key_mat[key2][0] == key1, (
        "{:d}-{:d} is not a distance coordinate in this zmatrix:\n{}"
        .format(key1, key2, string(zma, one_indexed=False)))
    name = name_mat[key2][0]

    return name


def central_angle_coordinate_name(zma, key1, key2, key3):
    """ get the name of angle coordinate for a set of 3 atoms

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key
        :type key1: int
        :param key2: the second key (central atom)
        :type key2: int
        :param key3: the third key
        :type key3: int
        :rtype: str
    """

    key1, key3 = sorted([key1, key3])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert key_mat[key3][0] == key2 and key_mat[key3][1] == key1, (
        "{:d}-{:d}-{:d} is not a distance coordinate in this zmatrix:\n{}"
        .format(key1, key2, key3, string(zma, one_indexed=False)))
    name = name_mat[key3][1]

    return name


def dihedral_angle_coordinate_name(zma, key1, key2, key3, key4):
    """ get the name of dihedral coordinate for a set of 4 atoms

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key
        :type key1: int
        :param key2: the second key
        :type key2: int
        :param key3: the third key
        :type key3: int
        :param key4: the fourth key
        :type key4: int
        :rtype: str
    """

    if key1 > key4:
        key1, key2, key3, key4 = key4, key3, key2, key1

    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert (
        key_mat[key4][0] == key3 and key_mat[key4][1] == key2 and
        key_mat[key4][2] == key1
    ), (
        "{:d}-{:d}-{:d}-{:d} is not a dihedral coordinate in this zmatrix:\n{}"
        .format(key1, key2, key3, key4, string(zma, one_indexed=False)))

    name = name_mat[key4][2]

    return name


def torsion_coordinate_name(zma, key1, key2, zgra=None):
    """ Obtain the name for dihedral coordinate about a torsion axis
        (rotational bond).

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key in the torsion axis (rotational bond)
        :type key1: int
        :param key2: the second key in the torsion axis (rotational bond)
        :type key2: int
        :param gra: an automol graph data structure, aligned to the z-matrix;
            used to check connectivity when necessary
        :rtype: str
    """

    key = torsion_leading_atom(zma, key1, key2, zgra=zgra)
    name_mat = name_matrix(zma)
    name = name_mat[key][-1]

    return name


def torsion_leading_atom(zma, key1, key2, zgra=None):
    """ Obtain the leading atom for a torsion coordinate about a torsion axis.

        The leading atom is the atom whose dihedral defines the torsional
        coordinate, which must always be the first dihedral coordinate
        for this bond.

        A bond is properly decoupled if all other dihedrals along this
        bond depend on the leading atom.

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key in the torsion axis (rotational bond)
        :type key1: int
        :param key2: the second key in the torsion axis (rotational bond)
        :type key2: int
        :param gra: an automol graph data structure, aligned to the z-matrix;
            used to check connectivity when necessary
        :rtype: int
    """

    key_mat = key_matrix(zma)
    krs1 = [(key, row) for key, row in enumerate(key_mat)
            if row[:2] == (key1, key2)]
    krs2 = [(key, row) for key, row in enumerate(key_mat)
            if row[:2] == (key2, key1)]

    lead_key_candidates = []

    for krs in (krs1, krs2):
        if krs:
            keys, rows = zip(*krs)
            start_key = keys[0]
            assert all(row[-1] == start_key for row in rows[1:]), (
                "Torsion coordinate along bond {:d}-{:d} not decoupled:\n{}"
                .format(key1, key2, string(zma, one_indexed=False)))
            if rows[0][-1] is not None:
                lead_key_candidates.append(start_key)

    if not lead_key_candidates:
        lead_key = None
    elif len(lead_key_candidates) == 1:
        lead_key = lead_key_candidates[0]
    else:
        # If we get to this point, then the z-matrix includes dihedrals across
        # the key1-key2 bond in both directions and we have to choose which
        # dihedral to use. This mans there will be two lead_key_candidates.
        zgra = automol.convert.zmat.graph(zma) if zgra is None else zgra

        # Let key0 be the lead key and let (key1, key2, key3) be its key row in
        # the z-matrix. For the torsion coordinate, key0-key1-key2-key3 should
        # all be connected in a line. For subsidiary dihedral coordinates, key3
        # will be connected to key1 instead of key2.
        # A simple solution is therefore to choose the lead key based on
        # whether or not key2 and key3 are connected, which is what this code
        # does.
        bnd_keys = bond_keys(zgra)
        lead_key = next((k for k in lead_key_candidates if
                         frozenset(key_mat[k][-2:]) in bnd_keys), None)

        # If that fails, choose the key that appears earlier. It's possible
        # that it would be better to choose the later one, in which case we
        # would replace the min() here with a max().
        if lead_key is None:
            lead_key = min(lead_key_candidates)

    return lead_key
