""" z-matrix conversions
"""
import itertools
import numpy
from phydat import phycon
from automol import create
from automol import util
from automol.convert import _pyx2z
from automol.convert import _util
from automol.vmat import symbols
from automol.vmat import key_matrix
from automol.vmat import name_matrix
from automol.vmat import standard_name_matrix
from automol.vmat import standard_names
from automol.graph._graph import standard_keys_without_dummy_atoms
from automol.graph._graph import add_bonds
from automol.graph.geom import is_atom
import automol.graph.vmat
from automol.convert.geom import (
    without_dummy_atoms as geom_without_dummy_atoms)
from automol.convert.geom import graph as geom_graph
from automol.convert.geom import connectivity_graph as geom_connectivity_graph
from automol.convert.geom import from_subset
from automol.convert.geom import insert_dummies_on_linear_atoms
from automol.convert.geom import central_angle as _central_angle
from automol.convert.geom import dihedral_angle as _dihedral_angle
from automol.convert.geom import distance as _distance
from automol.graph.geom import symbols as geom_symbols


# z-matrix => geometry
def geometry(zma):
    """ Convert a Z-Matrix to a molecular geometry that precludes dummy atoms.
        A dictionary is generated showing how dummy atoms are added:

            dummy_key_dct = {atom_key: dummy_atom_key}

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: (automol molecular geometry data structure, dict[int, int])
    """
    geo = geometry_with_dummy_atoms(zma)
    gra = connectivity_graph(zma, dummy=True)
    gra, dummy_key_dct = standard_keys_without_dummy_atoms(gra)
    geo = geom_without_dummy_atoms(geo)

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

    geo = create.geom.from_data(syms, xyzs)

    return geo


# z-matrix => graph
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
        geo = geom_without_dummy_atoms(geo)
    gra = geom_graph(geo, stereo=stereo)

    if dummy:
        bnd_keys = tuple(dummy_neighbor_keys(zma).items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = add_bonds(gra, bnd_keys, ord_dct=ord_dct)

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
        geo = geom_without_dummy_atoms(geo)

    gra = geom_connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)

    if dummy:
        bnd_keys = tuple(dummy_neighbor_keys(zma).items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = add_bonds(gra, bnd_keys, ord_dct=ord_dct)

    return gra


# z-matrix => formula
def formula(zma):
    """ Generate a stoichiometric formula dictionary from a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :type: dict[str: int]
    """

    syms = symbols(zma)
    fml = _util.formula(syms)

    return fml
#
#
# if __name__ == '__main__':
#     ZMA = automol.zmat.from_string("""
#             O
#             C  0    R1
#             H  0    R2  1    A2
#             X  1    R3  0    A3  2    D3
#             H  1    R4  3    A4  0    D4
#             H  1    R5  3    A5  0    D5
#             H  1    R6  3    A6  0    D6
#             C  1    R7  3    A7  0    D7
#             C  7    R8  1    A8  3    D8
#             H  7    R9  1    A9  8    D9
#             H  7   R10  1   A10  8   D10
#             H  8   R11  7   A11  1   D11
#             H  8   R12  7   A12  11   D12
#             H  8   R13  7   A13  11   D13
#
#             R1   =   1.442918
#             R2   =   0.967297
#             R3   =   1.000000
#             R4   =   1.081914
#             R5   =   1.086005
#             R6   =   1.087359
#             R7   =   2.010115
#             R8   =   1.479275
#             R9   =   1.080107
#             R10  =   1.096435
#             R11  =   1.089456
#             R12  =   1.107331
#             R13  =   1.101396
#             A2   = 105.655170
#             A3   =  90.000000
#             A4   =  56.389133
#             A5   =  66.442680
#             A6   = 168.250153
#             A7   =  97.179291
#             A8   =  91.208948
#             A9   =  74.145648
#             A10  = 122.233661
#             A11  = 113.127728
#             A12  = 109.910642
#             A13  = 113.513335
#             D3   =   0.000000
#             D4   = 112.709994
#             D5   = 253.876418
#             D6   = 151.399711
#             D7   = 184.721216
#             D8   = 263.523904
#             D9   = 230.715641
#             D10  = 123.944057
#             D11  = 283.002678
#             D12  = 239.994535
#             D13  = 122.854079
#     """, one_indexed=False)
#
#     print(automol.zmat.string(ZMA, one_indexed=False))
#     print(automol.zmat.dummy_keys(ZMA))
#     print(automol.zmat.dummy_neighbor_keys(ZMA))
#     GRA1 = connectivity_graph(ZMA, dummy=False)
#     GRA2 = connectivity_graph(ZMA, dummy=True)
#     print(automol.graph.string(GRA1, one_indexed=False))
#     print(automol.graph.string(GRA2, one_indexed=False))
#     GEO, DUMMY_KEY_DCT = geometry(ZMA)
#     print(automol.geom.string(GEO))
#     print(DUMMY_KEY_DCT)


def value_matrix(zma, angstrom=False, degree=False):
    """ Obtain the value matrix of the Z-Matrix that contains the
        values of the coordinates by row and column.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: tuple(tuple(str))
    """

    if zma:
        val_mat = tuple(zip(*zma))[3]

        val_mat = [list(row) + [None]*(3-len(row)) for row in val_mat]
        val_mat = numpy.array(val_mat, dtype=numpy.object_)

        val_mat[1:, 0] *= phycon.BOHR2ANG if angstrom else 1
        val_mat[2:, 1] *= phycon.RAD2DEG if degree else 1
        val_mat[3:, 2] *= phycon.RAD2DEG if degree else 1
    else:
        val_mat = ()

    return tuple(map(tuple, val_mat))


def dummy_keys(zma):
    """ Obtain keys to dummy atoms in the Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: tuple[int]
    """
    keys = tuple(key for key, sym in enumerate(symbols(zma)) if sym == 'X')
    return keys


def dummy_neighbor_keys(zma, geom_indexing=False):
    """ Obtain keys to dummy atoms in the Z-Matrix, along with their
        neighboring atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :returns: a dictionary with the dummy atoms as keys and their neighbors
            as values
        :rtype: dict[int: int]
    """
    key_mat = key_matrix(zma)
    dum_keys = dummy_keys(zma)
    key_dct = {}
    for dum_key in dum_keys:
        ngb_key = key_mat[dum_key][0]
        if ngb_key is None:
            ngb_key = next(row for row, (k, _, _) in enumerate(key_mat)
                           if k == dum_key)

        key_dct[dum_key] = ngb_key
    return key_dct


# geometry => z-matrix
def zmatrix(geo, ts_bnds=()):
    """ Generate a corresponding Z-Matrix for a molecular geometry
        using internal autochem procedures.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param ts_bnds: keys for the breaking/forming bonds in a TS
        :type ts_bnds: tuple(frozenset(int))
    """

    if ts_bnds:
        raise NotImplementedError

    if automol.geom.is_atom(geo):
        symbs = symbols(geo)
        key_mat = [[None, None, None]]
        val_mat = [[None, None, None]]
        zma = create.zmat.from_data(symbs, key_mat, val_mat)
        zma_keys = [0]
        dummy_key_dct = {}
    else:
        geo, dummy_key_dct = insert_dummies_on_linear_atoms(geo)
        gra = geom_connectivity_graph(geo)
        bnd_keys = tuple(dummy_key_dct.items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = add_bonds(gra, bnd_keys, ord_dct=ord_dct)
        vma, zma_keys = automol.graph.vmat.vmatrix(gra)
        geo = from_subset(geo, zma_keys)
        zma = from_geometry(vma, geo)

    return zma, zma_keys, dummy_key_dct


def zmatrix_x2z(geo, ts_bnds=()):
    """ Generate a corresponding Z-Matrix for a molecular geometry
        using x2z interface.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param ts_bnds: keys for the breaking/forming bonds in a TS
        :type ts_bnds: tuple(frozenset(int))
    """

    if is_atom(geo):
        symbs = geom_symbols(geo)
        key_mat = [[None, None, None]]
        val_mat = [[None, None, None]]
        zma = create.zmat.from_data(symbs, key_mat, val_mat)
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        zma = _pyx2z.to_zmatrix(x2m)
    zma = standard_form(zma)

    return zma


def zmatrix_torsion_coordinate_names(geo, ts_bnds=()):
    """ Generate a list of torsional coordinates using x2z interface. These
        names corresond to the Z-Matrix generated using the same algorithm.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param ts_bnds: keys for the breaking/forming bonds in a TS
        :type ts_bnds: tuple(frozenset(int))
        :rtype: tuple(str)
    """

    symbs = geom_symbols(geo)
    if len(symbs) == 1:
        names = ()
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        names = _pyx2z.zmatrix_torsion_coordinate_names(x2m)

        zma = _pyx2z.to_zmatrix(x2m)
        name_dct = standard_names(zma)
        names = tuple(map(name_dct.__getitem__, names))

    return names


def zmatrix_atom_ordering(geo, ts_bnds=()):
    """ Generate a dictionary which maps the order of atoms from the input
        molecular geometry to the order of atoms of the resulting Z-Matrix
        that is generated by the x2z interface.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param ts_bnds: keys for the breaking/forming bonds in a TS
        :type ts_bnds: tuple(frozenset(int))
        :rtype: dict[int: int]
    """

    symbs = geom_symbols(geo)
    if len(symbs) == 1:
        idxs = {0: 0}
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        idxs = _pyx2z.zmatrix_atom_ordering(x2m)

    return idxs


def from_geometry(vma, geo):
    """  Build a Z-Matrix from a V-Matrix and a molecular geometry.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: automol Z-Matrix data structure
    """
    assert symbols(vma) == symbols(geo)

    syms = symbols(vma)
    key_mat = key_matrix(vma)
    name_mat = name_matrix(vma)
    val_mat = numpy.empty(numpy.shape(key_mat), dtype=numpy.object_)

    for row, key_row in enumerate(key_mat):
        if row > 0:
            val_mat[row, 0] = _distance(geo, row, *key_row[:1])
        if row > 1:
            val_mat[row, 1] = _central_angle(geo, row, *key_row[:2])
        if row > 2:
            val_mat[row, 2] = _dihedral_angle(geo, row, *key_row[:3])

    zma = create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


def standard_form(zma, shift=0):
    """ Build a Z-Matrix where all of the coordinate names of an input Z-Matrix
        have been put into standard form:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
    """
    name_mat = standard_name_matrix(zma, shift=shift)
    return set_name_matrix(zma, name_mat)


def set_name_matrix(zma, name_mat):
    """ Re-set the name matrix of a Z-Matrix using the input name matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param name_mat: matrix of Z-Matrix coordinate names
        :type name_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    syms = symbols(zma)
    val_mat = value_matrix(zma)
    key_mat = key_matrix(zma)
    zma = create.zmat.from_data(syms, key_mat, val_mat, name_mat)

    return zma
