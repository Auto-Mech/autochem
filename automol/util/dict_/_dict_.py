""" Helper functions for working with Python dictionaries
"""

from itertools import permutations
from itertools import starmap as _starmap


def empty_if_none(obj):
    """ Returns an empty dictionary if input object is None,
        otherwise return the object.

        :param obj: generic object
        :type obj: any
    """
    return dict() if obj is None else obj


def right_update(dct1, dct2):
    """ Updates the entries of `dct1` with those of `dct2`.

        :param dct1: dictionary1 that will be updated
        :type dct1: dict
        :param dct2: dictionary2 whose entries will override dct1
        :type dct2: dict
        :rtype: dict
    """

    dct = {}
    dct.update(dct1)
    dct.update(dct2)
    return dct


def by_key(dct, keys, fill_val=None):
    """ dictionary on a set of keys, filling missing entries
    """
    return dict(zip(keys, values_by_key(dct, keys, fill_val=fill_val)))


def by_value(dct, func):
    """ dictionary on a set of values, determined by a function
    """
    keys = keys_by_value(dct, func)
    return by_key(dct, keys)


def values_by_key(dct, keys, fill_val=None):
    """ return dictionary values for specific keys, filling missing entries
    """
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def values_by_unordered_tuple(dct, key, fill_val=None):
    """ return dictionary values where keys are a tuple where either order
        of tuple will access element

        should really add a check if flipping key worder gives different vals
    """

    val = None
    nkeys = len(key)
    # vals = tuple(dct.get(key, None) for itertools.permutations(key, nkeys)))
    for _key in permutations(key, nkeys):
        val = dct.get(_key, None)
        if val is not None:
            break

    if val is None:
        val = fill_val

    return val


def values_in_multilevel_dct(dct, key1, key2, fill_val=None):
    """ Obtain a dictionary value where
        dct[key1][key2]
    """

    dct2 = dct.get(key1, None)
    if dct2 is not None:
        val = dct2.get(key2, fill_val)
    else:
        val = fill_val

    return val


def keys_by_value(dct, func):
    """ return dictionary keys for specific values
    """
    return frozenset(key for key, val in dct.items() if func(val))


def transform_keys(dct, func):
    """ apply a function to each key
    """
    return dict(zip(map(func, dct.keys()), dct.values()))


def transform_values(dct, func):
    """ apply a function to each value
    """
    return dict(zip(dct.keys(), map(func, dct.values())))


def transform_items_to_values(dct, func):
    """ apply a function to each value
    """
    return dict(zip(dct.keys(), _starmap(func, dct.items())))


def keys_sorted_by_value(dct):
    """ dictionary keys sorted by their associated values
    """
    return tuple(key for key, _ in sorted(dct.items(), key=lambda x: x[1]))


# def filter_by_key(dct, func):
#     """ filter dictionary entries by their keys
#     """
#     return {key: val for key, val in dct.items() if func(key)}


def filter_by_value(dct, func=lambda val: val):
    """ filter dictionary entries by their values
    """
    return {key: val for key, val in dct.items() if func(val)}


def merge_sequence(dcts):
    """ merge a sequence of dictionaries
    """
    merged_dct = {}
    for dct in dcts:
        merged_dct.update(dct)
    return merged_dct


def filter_keys(dct_1, dct_2):
    """ filters out from 1 all entries present in 2
    """

    keys_topop = list(dct_2.keys())

    for key in keys_topop:
        try:
            dct_1.pop(key, None)
        except KeyError:
            continue

    return dct_1
