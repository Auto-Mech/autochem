""" Helper functions for working with Python dictionaries
"""

import itertools
from _collections_abc import Callable
from copy import deepcopy
from typing import Any

import numpy


# for dictionaries do dict[object, tuple[object]] [obj, tuple [obj, ...]]
def invert(dct:dict[object, tuple[object]][object, tuple[object]])-> dict[object, tuple[object]][object, tuple[object]]:
    """Transposes the keys and values in a dictionary

    :param dct: Dictionary to transpose
    :return: Transposed dictionary
    """
    return {val: key for key, val in dct.items()}


def empty_if_none(obj:object):
    """Returns an empty dictionary if input object is None,
    otherwise return the object.
    :param obj: Generic object
    :return: Object of none if none
    """
    return {} if obj is None else obj


def compose(dct1:dict[object, tuple[object]][object, tuple[object]], 
            dct2:dict[object, tuple[object]][object, tuple[object]])-> dict[object, tuple[object]]:
    """Get the composition of `dct1` with `dct2`

    That is, dct[k] = dct1[dct2[k]]

    :param dct1: The first dictionary
    :param dct2: The second dictionary
    :returns: The composed dictionary
    """
    return {k2: dct1[v2] for k2, v2 in dct2.items()}


def by_key(dct:dict[object, tuple[object]], keys, fill:bool=True, fill_val:bool=None):
    """dictionary on a set of keys, filling missing entries

    :param fill: Fill in missing values?
    :type fill: bool
    :param fill_val: If `fill` is `True`, fill missing entries with this value.
    """
    if fill:
        dct = dict[object, tuple[object]](zip(keys, values_by_key(dct, keys, fill_val=fill_val)))
    else:
        dct = {key: dct[key] for key in keys if key in dct}
    return dct


def by_value(dct:dict[object, tuple[object]], func: Callable=lambda x: x):
    """dictionary on a set of values, determined by a function.
    :param dct: Dictionary of values
    :return: Dictionary with keys.
    """
    keys = keys_by_value(dct, func)
    return by_key(dct, keys)


def values_by_key(dct:dict[object, tuple[object]], keys:object, fill_val:object=None):
    """Return dictionary values for specific keys, filling missing entries"""
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def value_by_unordered_key(dct:dict[object, tuple[object]], key:object, fill_val=None):
    """Return the first value matching a tuple key in any order"""

    key_perms = itertools.permutations(key, len(key))
    val = next((dct[k] for k in key_perms if k in dct), fill_val)
    return val


def value_in_floatkey_dct(dct:dict[object, tuple[object]], key:object, tol=1.0e-5):
    """Access value in a dictionary that may have floats with large numbers of
    decimal places or strings
    """

    if isinstance(key, float):
        val = None
        for dkey in dct:
            if isinstance(dkey, float):
                if numpy.isclose(key, dkey, atol=tol):
                    val = dct[dkey]
    else:
        val = dct.get(key, None)

    return val


def values_in_multilevel_dct(dct:dict[object, tuple[object]], key1:object,
                              key2:object, fill_val=None):
    """Obtain a dictionary value where
    dct[key1][key2]
    """

    dct2 = dct.get(key1, None)
    if dct2 is not None:
        val = dct2.get(key2, fill_val)
    else:
        val = fill_val

    return val


def keys_by_value(dct:dict[object, tuple[object]], func: Callable =lambda x: x):
    """return dictionary keys for specific values"""
    return frozenset(key for key, val in dct.items() if func(val))


def transform_keys(dct:dict[object, tuple[object]], func:Callable=lambda x: x):
    """apply a function to each key"""
    return dict[object, tuple[object]](zip(map(func, dct.keys()), dct.values()))


def transform_values(dct:dict[object, tuple[object]], func:Callable=lambda x: x):
    """apply a function to each value"""
    return dict[object, tuple[object]](zip(dct.keys(), map(func, dct.values())))


def transform_items_to_values(dct:dict[object, tuple[object]], func:Callable=lambda x: x):
    """apply a function to each value"""
    return dict[object, tuple[object]](zip(dct.keys(), itertools.starmap(func, dct.items())))


def keys_sorted_by_value(dct:dict[object, tuple[object]]):
    """dictionary keys sorted by their associated values"""
    return tuple(key for key, _ in sorted(dct.items(), key=lambda x: x[1]))


def values_sorted_by_key(dct:dict[object, tuple[object]]):
    """dictionary values sorted by their associated keys"""
    return tuple(val for _, val in sorted(dct.items()))


def filter_by_key(dct:dict[object, tuple[object]], func:Callable=lambda x: x):
    """filter dictionary entries by their values"""
    return {key: val for key, val in dct.items() if func(key)}


def filter_by_value(dct:dict[object, tuple[object]], func:Callable=lambda x: x):
    """filter dictionary entries by their values"""
    return {key: val for key, val in dct.items() if func(val)}


def filter_keys(dct_1:dict[object, tuple[object]], 
                dct_2:dict[object, tuple[object]])-> dict[object, tuple[object]]:
    """Given two dictionaries (dct1 bigger dct2), filter out
    from 1 all the entries present in 2.

    :param dct1:
    :param dct2:
    :return: filtered dct1
    """

    dct_ret = deepcopy(dct_1)

    keys_topop = list(dct_2.keys())

    for key in keys_topop:
        try:
            dct_ret.pop(key, None)
        except KeyError:
            continue

    return dct_ret


def merge_sequence(dcts:dict[object, tuple[object]]):
    """merge a sequence of dictionaries"""
    merged_dct = {}
    for dct in dcts:
        merged_dct.update(dct)
    return merged_dct


def sort_value_(dct:dict[object, tuple[object]], allow_missing: bool = True, missing_val: Any = None):
    """Generate a sort value function from a dictionary

    :param dct: A dictionary
    :param allow_missing: Allow missing values?, defaults to True
    :return:
    """

    def sort_value(key:object):
        assert allow_missing or key in dct, "No key {key} in dictionary:\n{dict[object, tuple[object]]}"
        return dct[key] if key in dct else missing_val

    return sort_value
