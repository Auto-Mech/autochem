""" multi-valued dictionary helpers

dictionary values must all be tuples of the same length
"""

from collections.abc import Mapping as _Mapping
from automol.util.dict_._dict_ import values_by_key as _values_by_key
from automol.util.dict_._dict_ import transform_values as _transform_values


def is_multidict(mdct):
    """ is this a multi-valued dictionary?
    """
    assert isinstance(mdct, _Mapping)
    vals = mdct.values()
    return (not vals or all(isinstance(val, tuple) for val in vals) and
            len(set(map(len, vals))) == 1)


def by_key_by_position(mdct, keys, pos):
    """ position values, as a dictionary with these keys
    """
    assert is_multidict(mdct)
    assert set(keys) <= set(mdct.keys())
    dct = {}
    if keys:
        keys = list(keys)
        vals = [row[pos] for row in _values_by_key(mdct, keys)]
        dct = dict(zip(keys, vals))
    return dct


def set_by_key_by_position(mdct, dct, pos):
    """ set values by position and key
    """
    assert is_multidict(mdct)
    assert set(dct.keys()) <= set(mdct.keys())
    if dct:
        mdct = _unfreeze(mdct)
        for key, val in dct.items():
            mdct[key][pos] = val
        mdct = _freeze(mdct)
    return mdct


def _unfreeze(mdct):
    return dict(zip(mdct.keys(), map(list, mdct.values())))


def _freeze(mdct):
    return _transform_values(mdct, func=tuple)
