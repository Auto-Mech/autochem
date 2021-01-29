""" z-matrix conversions
"""
import itertools
import numpy
from automol import create
from automol import util
from automol.convert import _util
import automol.zmat


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
    gra, dummy_key_dct = automol.graph.standard_keys_without_dummy_atoms(gra)
    geo = automol.geom.without_dummy_atoms(geo)

    return geo, dummy_key_dct


def geometry_with_dummy_atoms(zma):
    """ Convert a Z-Matrix to a molecular geometry that includes dummy atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: automol molecular geometry data structure
    """

    syms = automol.zmat.symbols(zma)

    natms = len(syms)
    key_mat = automol.zmat.key_matrix(zma)
    val_mat = automol.zmat.value_matrix(zma)

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
        geo = automol.geom.without_dummy_atoms(geo)

    gra = automol.convert.geom.graph(geo, stereo=stereo)

    if dummy:
        bnd_keys = tuple(automol.zmat.dummy_neighbor_keys(zma).items())
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

    gra = automol.convert.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)

    if dummy:
        bnd_keys = tuple(automol.zmat.dummy_neighbor_keys(zma).items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = automol.graph.add_bonds(gra, bnd_keys, ord_dct=ord_dct)

    return gra


# z-matrix => formula
def formula(zma):
    """ Generate a stoichiometric formula dictionary from a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :type: dict[str: int]
    """

    syms = automol.zmat.symbols(zma)
    fml = _util.formula(syms)

    return fml


if __name__ == '__main__':
    import automol
    ZMA = automol.zmat.from_string("""
            O
            C  0    R1
            H  0    R2  1    A2
            X  1    R3  0    A3  2    D3
            H  1    R4  3    A4  0    D4
            H  1    R5  3    A5  0    D5
            H  1    R6  3    A6  0    D6
            C  1    R7  3    A7  0    D7
            C  7    R8  1    A8  3    D8
            H  7    R9  1    A9  8    D9
            H  7   R10  1   A10  8   D10
            H  8   R11  7   A11  1   D11
            H  8   R12  7   A12  11   D12
            H  8   R13  7   A13  11   D13

            R1   =   1.442918
            R2   =   0.967297
            R3   =   1.000000
            R4   =   1.081914
            R5   =   1.086005
            R6   =   1.087359
            R7   =   2.010115
            R8   =   1.479275
            R9   =   1.080107
            R10  =   1.096435
            R11  =   1.089456
            R12  =   1.107331
            R13  =   1.101396
            A2   = 105.655170
            A3   =  90.000000
            A4   =  56.389133
            A5   =  66.442680
            A6   = 168.250153
            A7   =  97.179291
            A8   =  91.208948
            A9   =  74.145648
            A10  = 122.233661
            A11  = 113.127728
            A12  = 109.910642
            A13  = 113.513335
            D3   =   0.000000
            D4   = 112.709994
            D5   = 253.876418
            D6   = 151.399711
            D7   = 184.721216
            D8   = 263.523904
            D9   = 230.715641
            D10  = 123.944057
            D11  = 283.002678
            D12  = 239.994535
            D13  = 122.854079
    """, one_indexed=False)

    print(automol.zmat.string(ZMA, one_indexed=False))
    print(automol.zmat.dummy_keys(ZMA))
    print(automol.zmat.dummy_neighbor_keys(ZMA))
    GRA1 = connectivity_graph(ZMA, dummy=False)
    GRA2 = connectivity_graph(ZMA, dummy=True)
    print(automol.graph.string(GRA1, one_indexed=False))
    print(automol.graph.string(GRA2, one_indexed=False))
    GEO, DUMMY_KEY_DCT = geometry(ZMA)
    print(automol.geom.string(GEO))
    print(DUMMY_KEY_DCT)
