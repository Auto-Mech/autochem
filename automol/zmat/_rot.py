""" functions for working with rotational bonds and groups
"""
import automol.convert.zmat
import automol.graph
from automol.zmat._zmat import key_matrix
from automol.zmat._zmat import name_matrix
from automol.zmat._zmat import string


def torsion_leading_atom(zma, key1, key2):
    """ leading atom for a torsion coordinate about a rotational bond

    :param zma: the z-matrix
    :param key1: the first key in the rotational bond
    :param key2: the second key in the rotational bond

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


def torsion_coordinate_name(zma, key1, key2):
    """ name for a torsion coordinate about a rotational bond

    :param zma: the z-matrix
    :param key1: the first key in the rotational bond
    :param key2: the second key in the rotational bond
    """
    key = torsion_leading_atom(zma, key1, key2)
    name_mat = name_matrix(zma)
    name = name_mat[key][-1]
    return name


# def torsion_coordinates(zma):
#     """ torsion coordinates for this z-matrix
#     """
#     gra = automol.convert.zmat.graph(zma)
#     rot_bnd_keys = automol.graph.rotational_bond_keys(gra)
#
#     key_mat = key_matrix(zma)
#     for row in key_mat:
#         if frozenset
#     print(automol.graph.string(gra, one_indexed=False))


if __name__ == '__main__':
    import automol
    ICH = automol.smiles.inchi('C=C(CC)CC=CC(O)(N)CC1CCC1')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.convert.geom.zmatrix_x2z(GEO)
    GEO = automol.convert.zmat.geometry(ZMA)
    NAMES = automol.geom.zmatrix_torsion_coordinate_names(GEO)
    COO_DCT = automol.vmat.coordinates(ZMA)
    COOS = [x[0] for x in map(COO_DCT.__getitem__, NAMES)]
    print(automol.zmat.string(ZMA, one_indexed=False))
    # print(torsion_coordinates(ZMA))
    print(COOS)
    GRA = automol.geom.graph(GEO)
    ROT_BND_KEYS = sorted(map(sorted, automol.graph.rotational_bond_keys(GRA)))

    TORS_NAMES = tuple(torsion_coordinate_name(ZMA, *k) for k in ROT_BND_KEYS)
    print(TORS_NAMES)
    print(TORS_NAMES == NAMES)
