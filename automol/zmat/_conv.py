""" Level 4 Z-Matrix functions
"""
import itertools
import numpy
from automol import util
import automol.graph
import automol.geom
from automol.zmat.base import symbols
from automol.zmat.base import key_matrix
from automol.zmat.base import name_matrix
from automol.zmat.base import value_matrix
from automol.zmat.base import string
from automol.zmat.base import dummy_neighbor_keys


# # conversions
def graph(zma, stereo=True, dummy=False):
    """ Convert a Z-Matrix to a molecular graph.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param dummy: parameter to include dummy atoms
        :type dummy: bool
        :rtype: (automol molecular geometry data structure, dict[int, int])
    """

    geo = geometry_with_dummy_atoms(zma)
    if not dummy:
        geo = automol.geom.without_dummy_atoms(geo)
    gra = automol.geom.graph(geo, stereo=stereo)

    if dummy:
        bnd_keys = tuple(dummy_neighbor_keys(zma).items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = automol.graph.add_bonds(gra, bnd_keys, ord_dct=ord_dct)

    return gra


def connectivity_graph(zma, dummy=False,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ Convert a Z-Matrix to a molecular connectivitiy graph.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param rqq_bond_max: maximum distance between heavy atoms
        :type rqq_bond_max: float
        :param rqh_bond_max: maximum distance between heavy atoms and hydrogens
        :type rqh_bond_max: float
        :param rhh_bond_max: maximum distance between hydrogens
        :type rhh_bond_max: float
        :param dummy: parameter to include dummy atoms
        :type dummy: bool
        :rtype: automol molecular graph data structure
    """

    geo = geometry_with_dummy_atoms(zma)
    if not dummy:
        geo = automol.geom.without_dummy_atoms(geo)

    gra = automol.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)

    if dummy:
        bnd_keys = tuple(dummy_neighbor_keys(zma).items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = automol.graph.add_bonds(gra, bnd_keys, ord_dct=ord_dct)

    return gra


def geometry(zma, dummy=False):
    """ Convert a Z-Matrix to a molecular geometry, with or without dummy atoms

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param dummy: include dummy atoms in the geometry?
        :type dummy: bool
        :returns: automol molecular geometry data structure
    """
    if dummy:
        geo = geometry_with_dummy_atoms(zma)
    else:
        geo, _ = geometry_with_conversion_info(zma)

    return geo


def geometry_with_conversion_info(zma):
    """ Convert a Z-Matrix to a molecular geometry that precludes dummy atoms.
        A dictionary is generated showing how dummy atoms are added:

            dummy_key_dct = {atom_key: dummy_atom_key}

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :returns: automol molecular geometry data structure
    """
    geo = geometry_with_dummy_atoms(zma)
    gra = connectivity_graph(zma, dummy=True)
    gra, dummy_key_dct = automol.graph.shift_remove_dummy_atoms(gra)
    geo = automol.geom.without_dummy_atoms(geo)

    return geo, dummy_key_dct


def geometry_with_dummy_atoms(zma):
    """ Convert a Z-Matrix to a molecular geometry that includes dummy atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: automol molecular geometry data structure
    """

    syms = symbols(zma)

    natms = len(syms)
    key_mat = key_matrix(zma)
    val_mat = value_matrix(zma)

    xyzs = numpy.zeros((natms, 3))

    for key in range(1, natms):
        vals = val_mat[key][:min(key, 3)]
        keys = key_mat[key][:min(key, 3)]
        ref_xyzs = xyzs[list(keys)]
        xyz = util.vec.from_internals(*itertools.chain(*zip(vals, ref_xyzs)))
        xyzs[key] = xyz

    geo = automol.geom.from_data(syms, xyzs)

    return geo


# # derived properties
def distance(zma, key1, key2, angstrom=False):
    """ Measure the distance between two atoms defined in a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the pair to be measured
        :type key1: int
        :param key2: key of atom 2 in the pair to be measured
        :type key2: int
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
    """
    geo = geometry_with_dummy_atoms(zma)
    return automol.geom.distance(geo, key1, key2, angstrom=angstrom)


def central_angle(zma, key1, key2, key3, degree=False):
    """ Measure the angle inscribed by three atoms in a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the triplet to be measured
        :type key1: int
        :param key2: key of atom 2 in the triplet to be measured
        :type key2: int
        :param key3: key of atom 3 in the triplet to be measured
        :type key3: int
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
    """
    geo = geometry_with_dummy_atoms(zma)
    return automol.geom.central_angle(geo, key1, key2, key3, degree=degree)


def dihedral_angle(zma, key1, key2, key3, key4, degree=False):
    """ Measure the angle inscribed by four atoms in a molecular geometry.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the quartet to be measured
        :type key1: int
        :param key2: key of atom 2 in the quartet to be measured
        :type key2: int
        :param key3: key of atom 3 in the quartet to be measured
        :type key3: int
        :param key4: key of atom 4 in the quartet to be measured
        :type key4: int
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
    """
    geo = geometry_with_dummy_atoms(zma)
    return automol.geom.dihedral_angle(geo, key1, key2, key3, key4,
                                       degree=degree)


# # torsions
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
                "Torsion coordinate along bond "
                f"{key1:d}-{key2:d} not decoupled:\n"
                f"{string(zma, one_indexed=False)}")
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
        zgra = graph(zma) if zgra is None else zgra

        # Let key0 be the lead key and let (key1, key2, key3) be its key row in
        # the z-matrix. For the torsion coordinate, key0-key1-key2-key3 should
        # all be connected in a line. For subsidiary dihedral coordinates, key3
        # will be connected to key1 instead of key2.
        # A simple solution is therefore to choose the lead key based on
        # whether or not key2 and key3 are connected, which is what this code
        # does.
        bnd_keys = automol.graph.bond_keys(zgra)
        lead_key = next((k for k in lead_key_candidates if
                         frozenset(key_mat[k][-2:]) in bnd_keys), None)

        # If that fails, choose the key that appears earlier. It's possible
        # that it would be better to choose the later one, in which case we
        # would replace the min() here with a max().
        if lead_key is None:
            lead_key = min(lead_key_candidates)

    return lead_key
