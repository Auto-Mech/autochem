""" some math functions
"""
import numpy
from ._core import symbols as _symbols
from ._core import coordinate_matrix as _coordinate_matrix
from ._core import value_matrix as _value_matrix


def almost_equal(zma1, zma2):
    """ are these z-matrices almost equal?
    """
    ret = False
    if _symbols(zma1) == _symbols(zma2):
        if _coordinate_matrix(zma1) == _coordinate_matrix(zma2):
            for shift in [0., numpy.pi/3., numpy.pi/5., numpy.pi/7.]:
                val_mat1 = _array_mod_angles(_value_matrix(zma1), shift=shift)
                val_mat2 = _array_mod_angles(_value_matrix(zma2), shift=shift)
                if numpy.allclose(val_mat1, val_mat2, equal_nan=True):
                    ret = True
                    break
    return ret


def _array_mod_angles(val_mat, shift=0.):
    val_mat = numpy.array(val_mat, dtype=float)
    val_mat[:, 1:3] = numpy.mod(val_mat[:, 1:3]+shift, 2*numpy.pi)
    return val_mat
