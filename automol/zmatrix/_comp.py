""" some comparison functions
"""
import numpy
from ._core import symbols as _symbols
from ._core import key_matrix as _key_matrix
from ._core import name_matrix as _name_matrix
from ._core import values as _values
from ._core import variable_names as _variable_names
from ._core import distance_value_names as _distance_value_names
from ._core import angle_value_names as _angle_value_names


def almost_equal(zma1, zma2, rtol=2e-5):
    """ are these z-matrices numerically equal?
    """
    ret = False
    if (_symbols(zma1) == _symbols(zma2) and
            _key_matrix(zma1) == _key_matrix(zma2) and
            _name_matrix(zma1) == _name_matrix(zma2) and
            _variable_names(zma1) == _variable_names(zma2)):
        val_dct1 = _values(zma1)
        val_dct2 = _values(zma2)
        dist_names = _distance_value_names(zma1)
        dist_vals1 = tuple(map(val_dct1.__getitem__, dist_names))
        dist_vals2 = tuple(map(val_dct2.__getitem__, dist_names))
        if numpy.allclose(dist_vals1, dist_vals2, rtol=rtol):
            ang_names = _angle_value_names(zma1)
            ang_vals1 = tuple(map(val_dct1.__getitem__, ang_names))
            ang_vals2 = tuple(map(val_dct2.__getitem__, ang_names))
            for shift in (0., numpy.pi/10.):
                ang_vals1 = numpy.mod(numpy.add(ang_vals1, shift), 2*numpy.pi)
                ang_vals2 = numpy.mod(numpy.add(ang_vals2, shift), 2*numpy.pi)
                if numpy.allclose(ang_vals1, ang_vals2, rtol=rtol):
                    ret = True
                    break
    return ret
