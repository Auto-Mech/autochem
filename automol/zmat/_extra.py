""" functions for working with torsion axis (rotational bond)s and groups
"""
from automol.zmat._zmat import key_matrix
from automol.zmat._zmat import name_matrix
from automol.zmat._zmat import string


def distance_coordinate_name(zma, key1, key2):
    """ get the name of a distance coordinate for a given bond
    """
    key1, key2 = sorted([key1, key2])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert key_mat[key2][0] == key1, (
        "{:d}-{:d} is not a distance coordinate in this zmatrix:\n{}"
        .format(key1, key2, string(zma, one_indexed=False)))
    name = name_mat[key2][0]
    return name


def torsion_coordinate_name(zma, key1, key2):
    """ name for dihedral coordinate about a torsion axis (rotational bond)

    :param zma: the z-matrix
    :param key1: the first key in the torsion axis (rotational bond)
    :param key2: the second key in the torsion axis (rotational bond)
    """
    key = torsion_leading_atom(zma, key1, key2)
    name_mat = name_matrix(zma)
    name = name_mat[key][-1]
    return name


def torsion_leading_atom(zma, key1, key2):
    """ leading atom for a torsion coordinate about a torsion axis

    :param zma: the z-matrix
    :param key1: the first key in the torsion axis (rotational bond)
    :param key2: the second key in the torsion axis (rotational bond)

    The leading atom is the atom whose dihedral defines the torsional
    coordinate, which must always be the first dihedral coordinate for this
    bond

    A bond is properly decoupled if all other dihedrals along this bond depend
    on the leading atom
    """
    key_mat = key_matrix(zma)
    krs1 = [(key, row) for key, row in enumerate(key_mat)
            if row[:2] == (key1, key2)]
    krs2 = [(key, row) for key, row in enumerate(key_mat)
            if row[:2] == (key2, key1)]

    lead_key = None

    for krs in (krs1, krs2):
        if krs:
            keys, rows = zip(*krs)
            start_key = keys[0]
            assert all(row[-1] == start_key for row in rows[1:]), (
                "Torsion coordinate along bond {:d}-{:d} not decoupled:\n{}"
                .format(key1, key2, string(zma, one_indexed=False)))
            if rows[0][-1] is not None:
                assert lead_key is None
                lead_key = start_key

    return lead_key
