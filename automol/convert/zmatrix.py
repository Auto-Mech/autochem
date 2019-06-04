""" z-matrix conversions
"""
import itertools
import numpy
from automol import create
from automol import cart
import automol.convert.geom
import automol.geom
import automol.zmatrix


# z-matrix => geometry
def geometry(zma, remove_ghost_atoms=None):
    """ z-matrix => geometry
    """
    syms = automol.zmatrix.symbols(zma)

    natms = len(syms)
    key_mat = automol.zmatrix.key_matrix(zma)
    val_mat = automol.zmatrix.value_matrix(zma)

    xyzs = numpy.zeros((natms, 3))

    for key in range(1, natms):
        vals = val_mat[key][:min(key, 3)]
        keys = key_mat[key][:min(key, 3)]
        ref_xyzs = xyzs[list(keys)]
        xyz = cart.vec.from_internals(*itertools.chain(*zip(vals, ref_xyzs)))
        xyzs[key] = xyz

    geo = create.geom.from_data(syms, xyzs)
    if remove_ghost_atoms:
        geo = automol.geom.without_ghost_atoms(geo)

    return geo
