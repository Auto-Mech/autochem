""" z-matrix to cartesian geometry conversion

borrows from https://github.com/tmpchem/computational_chemistry
"""
import numpy
from automol.constructors.geom import from_data as _geom_from_data
from automol.cart.vec import unit_direction as _unit_direction
from automol.cart.vec import unit_perpendicular as _unit_perpendicular
from automol.zmatrix._core import symbols as _symbols
from automol.zmatrix._core import key_matrix as _key_matrix
from automol.zmatrix._core import value_matrix as _value_matrix


def geometry(zma):
    """ geometry from a z-matrix
    """
    syms = _symbols(zma)

    natms = len(syms)
    key_mat = _key_matrix(zma)
    val_mat = _value_matrix(zma)

    xyzs = numpy.zeros((natms, 3))

    for key in range(1, natms):
        vals = val_mat[key][:min(key, 3)]
        keys = key_mat[key][:min(key, 3)]
        ref_xyzs = xyzs[list(keys)]
        ref_axes = _local_axes(*ref_xyzs)
        loc_xyz = _local_position(*vals)
        xyzs[key] = ref_xyzs[0] + numpy.dot(loc_xyz, ref_axes)

    geo = _geom_from_data(syms, xyzs)
    return geo


def _local_position(dist=0., ang=0., dih=0.):
    """ position by internal coordinates in the local axis frame
    """
    x_comp = dist * numpy.sin(ang) * numpy.sin(dih)
    y_comp = dist * numpy.sin(ang) * numpy.cos(dih)
    z_comp = dist * numpy.cos(ang)
    return (x_comp, y_comp, z_comp)


def _local_axes(xyz1=(0., 0., 0.), xyz2=(0., 0., 1.), xyz3=(0., 1., 0.)):
    """ local axes for defining bond, angle, dihedral from support atoms
    """
    uxyz12 = _unit_direction(xyz1, xyz2)
    uxyz23 = _unit_direction(xyz2, xyz3)
    uxyz123_perp = _unit_perpendicular(uxyz23, uxyz12)

    z_ax = tuple(uxyz12)
    y_ax = tuple(_unit_perpendicular(uxyz12, uxyz123_perp))
    x_ax = tuple(_unit_perpendicular(y_ax, z_ax))
    return (x_ax, y_ax, z_ax)
