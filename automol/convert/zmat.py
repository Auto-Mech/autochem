""" z-matrix conversions
"""
import itertools
import numpy
from automol import create
from automol import util
from automol.convert import _util
import automol.zmat


# z-matrix => geometry
def geometry(zma, remove_dummy_atoms=False):
    """ z-matrix => geometry
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
    if remove_dummy_atoms:
        geo = automol.geom.without_dummy_atoms(geo)

    return geo


# z-matrix => graph
def graph(zma, remove_stereo=False):
    """ z-matrix => graph
    """
    geo = geometry(zma)
    gra = automol.convert.geom.graph(geo, remove_stereo=remove_stereo)
    return gra


def connectivity_graph(zma,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ z-matrix => connectivity graph
    """
    geo = geometry(zma)
    gra = automol.convert.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)
    return gra


# z-matrix => formula
def formula(zma):
    """ z-matrix => formula
    """
    syms = automol.zmat.symbols(zma)
    fml = _util.formula(syms)
    return fml
