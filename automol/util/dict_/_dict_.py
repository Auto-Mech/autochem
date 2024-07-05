""" Helper functions for working with Python dictionaries
"""
import itertools
from copy import deepcopy
from typing import Any

import numpy


def invert(dct):
    """Transposes the keys and values in a dictionary

    :param dct: dictionary to transpose
    :tupe dct: dict
    :rtype: dict
    """
    return {val: key for key, val in dct.items()}


def empty_if_none(obj):
    """Returns an empty dictionary if input object is None,
    otherwise return the object.

    :param obj: generic object
    :type obj: any
    """
    return {} if obj is None else obj


def compose(dct1, dct2):
    """Get the composition of `dct1` with `dct2`

    That is, dct[k] = dct1[dct2[k]]

    :param dct1: The first dictionary
    :type dct1: dict
    :param dct2: The second dictionary
    :type dct2: dict
    :returns: The composed dictionary
    :rtype: dict
    """
    return {k2: dct1[v2] for k2, v2 in dct2.items()}


def right_update(dct1, dct2):
    """Updates the entries of `dct1` with those of `dct2`.

    :param dct1: dictionary1 that will be updated
    :type dct1: dict
    :param dct2: dictionary2 whose entries will override dct1
    :type dct2: dict
    :rtype: dict
    """

    dct = {}
    dct1 = empty_if_none(dct1)
    dct2 = empty_if_none(dct2)
    dct.update(dct1)
    dct.update(dct2)
    return dct


def by_key(dct, keys, fill=True, fill_val=None):
    """dictionary on a set of keys, filling missing entries

    :param fill: Fill in missing values?
    :type fill: bool
    :param fill_val: If `fill` is `True`, fill missing entries with this value.
    """
    if fill:
        dct = dict(zip(keys, values_by_key(dct, keys, fill_val=fill_val)))
    else:
        dct = {key: dct[key] for key in keys if key in dct}
    return dct


def by_value(dct, func=lambda x: x):
    """dictionary on a set of values, determined by a function"""
    keys = keys_by_value(dct, func)
    return by_key(dct, keys)


def values_by_key(dct, keys, fill_val=None):
    """return dictionary values for specific keys, filling missing entries"""
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def value_by_unordered_key(dct, key, fill_val=None):
    """return the first value matching a tuple key in any order"""

    key_perms = itertools.permutations(key, len(key))
    val = next((dct[k] for k in key_perms if k in dct), fill_val)
    return val


def value_in_floatkey_dct(dct, key, tol=1.0e-5):
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


def values_in_multilevel_dct(dct, key1, key2, fill_val=None):
    """Obtain a dictionary value where
    dct[key1][key2]
    """

    dct2 = dct.get(key1, None)
    if dct2 is not None:
        val = dct2.get(key2, fill_val)
    else:
        val = fill_val

    return val


def separate_subdct(dct, key="global"):
    """Pulls out a sub-dictionary indexed by the given key and returns it
    and the original dictioanry with the requested sub-dictionary removed
    """

    # Grab the sub-dictonary and
    sub_dct = dct.get(key, {})

    # Build a new copy of the input dictionary with the sub-dict removed
    dct2 = deepcopy(dct)
    dct2.pop(key, None)

    return dct2, sub_dct


def merge_subdct(dct, key="global", keep_subdct=False):
    """Obtain a sub-dictionary indexed by a given key and merge its contents
    with all of the other sub-dictionaries in the main dictionary.

    :param dct: dictionary containing several sub-dictionaries
    :type dct: dict[str: dict]
    :param key: key for the sub-dictionary to be merged
    :type key: str, int, float, tuple
    :param keep_subdct: keep/remove the sub-dictionary following the merge
    :type keep_subdct: bool
    """

    if list(dct.keys()) == [key]:
        new_dct = {key: dct[key]}
    else:
        new_dct, sub_dct = separate_subdct(dct, key=key)
        for new_key in new_dct:
            new_dct[new_key] = right_update(sub_dct, new_dct[new_key])
        if keep_subdct:
            new_dct[key] = sub_dct

    return new_dct


def keys_by_value(dct, func=lambda x: x):
    """return dictionary keys for specific values"""
    return frozenset(key for key, val in dct.items() if func(val))


def transform_keys(dct, func=lambda x: x):
    """apply a function to each key"""
    return dict(zip(map(func, dct.keys()), dct.values()))


def transform_values(dct, func=lambda x: x):
    """apply a function to each value"""
    return dict(zip(dct.keys(), map(func, dct.values())))


def transform_items_to_values(dct, func=lambda x: x):
    """apply a function to each value"""
    return dict(zip(dct.keys(), itertools.starmap(func, dct.items())))


def keys_sorted_by_value(dct):
    """dictionary keys sorted by their associated values"""
    return tuple(key for key, _ in sorted(dct.items(), key=lambda x: x[1]))


def values_sorted_by_key(dct):
    """dictionary values sorted by their associated keys"""
    return tuple(val for _, val in sorted(dct.items()))


def filter_by_key(dct, func=lambda x: x):
    """filter dictionary entries by their values"""
    return {key: val for key, val in dct.items() if func(key)}


def filter_by_value(dct, func=lambda x: x):
    """filter dictionary entries by their values"""
    return {key: val for key, val in dct.items() if func(val)}


def filter_keys(dct_1, dct_2):
    """Given two dictionaries (dct1 bigger dct2), filter out
    from 1 all the entries present in 2.

    :param dct1:
    :param dct2:
    :return: filtered dct1
    :rtype: dict[]
    """

    dct_ret = deepcopy(dct_1)

    keys_topop = list(dct_2.keys())

    for key in keys_topop:
        try:
            dct_ret.pop(key, None)
        except KeyError:
            continue

    return dct_ret


def merge_sequence(dcts):
    """merge a sequence of dictionaries"""
    merged_dct = {}
    for dct in dcts:
        merged_dct.update(dct)
    return merged_dct


def sort_value_(dct, allow_missing: bool = True, missing_val: Any=None):
    """Generate a sort value function from a dictionary

    :param dct: A dictionary
    :type dct: dict
    :param allow_missing: Allow missing values?, defaults to True
    :type allow_missing: bool, optional
    :param missing_val: Value to assign to missing values, defaults to None
    :type missing_val: Any, optional
    """
    def sort_value(key):
        assert allow_missing or key in dct, "No key {key} in dictionary:\n{dict}"
        return dct[key] if key in dct else missing_val
    return sort_value
