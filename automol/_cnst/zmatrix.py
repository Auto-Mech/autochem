""" zmarix constructor
"""
import numpy
from .. import atom as _atom
from .. import _units


def from_data(symbols, distance_column, angle_column, torsion_column,
              angstroms=False, one_indexed=False):
    """ zmatrix data structure from the usual zmatrix columns
    """
    syms = list(map(_atom.standard_case, symbols))
    assert all(sym in _atom.SYMBOLS for sym in syms)
    natms = len(syms)

    dst_keys, dst_vals = _column_keys_and_values(
        distance_column, 1, natms, angstroms, one_indexed)
    ang_keys, ang_vals = _column_keys_and_values(
        angle_column, 2, natms, angstroms, one_indexed)
    tor_keys, tor_vals = _column_keys_and_values(
        torsion_column, 3, natms, angstroms, one_indexed)
    zma = tuple(zip(syms,
                    zip(dst_keys, ang_keys, tor_keys),
                    zip(dst_vals, ang_vals, tor_vals)))
    return zma


def _column_keys_and_values(column, num, natms, angstroms, one_indexed):
    if natms - num > 0:
        assert numpy.ndim(column) == 2
        assert numpy.shape(column) == (natms - num, 2)
        keys, vals = zip(*column)
        assert all(map(float.is_integer, map(float, keys)))
        idx = 1 if one_indexed else 0
        keys = tuple(int(key) - idx for key in keys)
        conv = (1. if not angstroms else _units.BOHR2ANG if num == 1 else
                _units.DEG2RAD)
        vals = tuple(float(val) * conv for val in vals)
        assert all(key < ref_key + num for ref_key, key in enumerate(keys))
        keys = (None,) * num + keys
        vals = (None,) * num + vals
    else:
        keys = vals = (None,) * min(natms, num)
    return keys, vals
