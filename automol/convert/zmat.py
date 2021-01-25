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
    gra = automol.convert.geom.connectivity_graph(geo, dummy=True)
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
def graph(zma, stereo=True):
    """ Convert a Z-Matrix to a molecular graph.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: (automol molecular geometry data structure, dict[int, int])
    """

    geo, _ = geometry(zma)
    gra = automol.convert.geom.graph(geo, stereo=stereo)

    return gra


def connectivity_graph(zma,
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
        :rtype: automol molecular graph data structure
    """

    geo = geometry(zma)
    gra = automol.convert.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)

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
