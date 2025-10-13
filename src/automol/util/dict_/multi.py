"""multi-valued dictionary helpers.

dictionary values must all be tuples of the same length
"""


from collections.abc import Mapping as _Mapping

from ._dict_ import transform_values as _transform_values
from ._dict_ import values_by_key as _values_by_key

MultiDict = dict[object, tuple[object, ...]]


def is_multidict(mdct: MultiDict) -> bool:
    """Is this a multi-valued dictionary?.
    :param mdict: Dictionary
    :return: True if it is a multi-valued dictionary, False if not.
    """
    assert isinstance(mdct, _Mapping)
    vals = mdct.values()
    return (
        not vals
        or all(isinstance(val, tuple) for val in vals)
        and len(set(map(len, vals))) == 1
    )


def by_key_by_position(mdct: MultiDict, keys: object, pos: int) -> dict[object, object]:
    """Position values, as a dictionary with these keys.
    :param mdct: Multi-valued dictionary
    :param keys: Keys for sorting
    :param pos: Position to select column for dictionary.
    """
    assert is_multidict(mdct)
    assert set(keys) <= set(mdct.keys())
    dct = {}
    if keys:
        keys = list(keys)
        vals = [row[pos] for row in _values_by_key(mdct, keys)]
        dct = dict(zip(keys, vals, strict=True))
    return dct


def set_by_key_by_position(mdct: MultiDict, dct: dict, pos: int) -> MultiDict:
    """Set values by position and key.
    :param mdct: Multi-valued dictionary
    :param dct: Dictionary
    :param pos:Position to select column for dictionary.
    :return: New multi-valued dictionary set by key and position.
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
    return dict(zip(mdct.keys(), map(list, mdct.values()), strict=True))


def _freeze(mdct):
    return _transform_values(mdct, func=tuple)
