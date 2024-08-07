"""Helper functions for working with Python dictionaries."""

import itertools
from collections.abc import Callable, Sequence
from copy import deepcopy
from typing import Any

import numpy


# for dictionaries do dict[object,object] [obj, tuple [obj, ...]]
def invert(dct: dict[object, object]) -> dict[object, object]:
    """Transposes the keys and values in a dictionary.

    :param dct: Dictionary to transpose
    :return: Transposed dictionary
    """
    return {val: key for key, val in dct.items()}


def empty_if_none(obj: object):
    """Returns an empty dictionary if input object is None,
    otherwise return the object.
    :param obj: Generic object
    :return: Object of none if none.
    """  # noqa: D401
    return {} if obj is None else obj


def compose(
    dct1: dict[object, object],
    dct2: dict[object, object],
) -> dict[object, object]:
    """Get the composition of `dct1` with `dct2`.

    That is, dct[k] = dct1[dct2[k]]

    :param dct1: The first dictionary
    :param dct2: The second dictionary
    :returns: The composed dictionary
    """
    return {k2: dct1[v2] for k2, v2 in dct2.items()}


def by_key(
    dct: dict[object, object],
    keys: object,
    fill: bool = True,
    fill_val: bool | None = None,
) -> dict[object, object]:
    """Dictionary on a set of keys, filling missing entries.

    :param fill: Fill in missing values
    :param keys: Keys for missing values
    :param fill_val: If `fill` is `True`, fill missing entries with this value.
    :return: 'fill; or dictionary with filled values
    """  # noqa: D401
    if fill:
        dct = dict[object, object](
            zip(keys, values_by_key(dct, keys, fill_val=fill_val), strict=False)
        )
    else:
        dct = {key: dct[key] for key in keys if key in dct}
    return dct


def by_value(dct: dict[object, object], func: Callable = lambda x: x):
    """Dictionary on a set of values, determined by a function.
    :param dct: Dictionary of values
    :param func: Function
    :return: Dictionary with keys.
    """  # noqa: D401
    keys: Sequence[object] = keys_by_value(dct, func)
    return by_key(dct, keys)


def values_by_key(dct: dict[object, object], keys: object, fill_val: object = None):
    """Return dictionary values for specific keys, filling missing entries.
    :param dct: Dictionary of values
    :param keys: Keys for entries
    :param fill_val: Fill where value is none
    :return: Dictionary with new keys.
    """
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def value_by_unordered_key(dct: dict[object, object], key: object, fill_val=None):
    """Return the first value matching a tuple key in any order.
    :param dct: Dictionary of values
    :param keys: Keys for entries
    :param fill_val: Values to fill the missing keys
    :return: Dictionary with new keys in order.
    """
    key_perms = itertools.permutations(key, len(key))
    val = next((dct[k] for k in key_perms if k in dct), fill_val)
    return val


def value_in_floatkey_dct(dct: dict[object, object], key: object, tol=1.0e-5):
    """Access value in a dictionary that may have floats with large numbers of
    decimal places or strings.
    :param dct: Dictionary of values
    :param keys: Keys for entries
    :param tol: Tolerance 1.0e-5
    :return: Dictionary with new values.
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


def values_in_multilevel_dct(
    dct: dict[object, object], key1: object, key2: object, fill_val=None
):
    """Obtain a dictionary value where
    dct[key1][key2].
    :param dict: Dictionary with two objects
    :param key1: First key
    :param key2: Second key
    :return: Value of the second entry.
    """
    dct2 = dct.get(key1, None)
    if dct2 is not None:
        val = dct2.get(key2, fill_val)
    else:
        val = fill_val

    return val


def keys_by_value(dct: dict[object, object], func: Callable = lambda x: x):
    """Return dictionary keys for specific values.
    :param dict: Dictionary
    :param func: Callable Function
    :return: Key that was called based on value.
    """
    return frozenset(key for key, val in dct.items() if func(val))


def transform_keys(dct: dict[object, object], func: Callable = lambda x: x):
    """Apply a function to each key.
    :param dct: Dictionary
    :param func: Callable function
    :return: Dictionary and new values.
    """
    return dict[object, object](zip(map(func, dct.keys()), dct.values(), strict=False))


def transform_values(dct: dict[object, object], func: Callable = lambda x: x):
    """Apply a function to each value.
    :param dct: Dictionary
    :param func: Callable function
    :return: New values in a dictionary.
    """
    return dict[object, object](zip(dct.keys(), map(func, dct.values()), strict=False))


def transform_items_to_values(dct: dict[object, object], func: Callable = lambda x: x):
    """Apply a function to each value.
    :param dct: Dictionary
    :param func: Callable function
    :return: Dictionary with new values.
    """
    return dict[object, object](
        zip(dct.keys(), itertools.starmap(func, dct.items()), strict=False)
    )


def keys_sorted_by_value(dct: dict[object, object]):
    """Dictionary keys sorted by their associated values.
    :param dct: Dictionary
    :param func: Callable function
    :return:Dictionary with sorted values.
    """  # noqa: D401
    return tuple(key for key, _ in sorted(dct.items(), key=lambda x: x[1]))


def values_sorted_by_key(dct: dict[object, object]):
    """Dictionary values sorted by their associated keys.
    :param dct: Dictionary
    :return: Dictionary with sorted values.
    """  # noqa: D401
    return tuple(val for _, val in sorted(dct.items()))


def filter_by_key(dct: dict[object, object], func: Callable = lambda x: x):
    """Filter dictionary entries by their values.
    :param dct: Dictionary
    :param func: Callable function
    :return:Dictionary with sorted values.
    """
    return {key: val for key, val in dct.items() if func(key)}


def filter_by_value(dct: dict[object, object], func: Callable = lambda x: x):
    """Filter dictionary entries by their values.
    :param dct: Dictionary
    :param func: Callable function
    :return:Dictionary with sorted values.
    """
    return {key: val for key, val in dct.items() if func(val)}


def filter_keys(
    dct_1: dict[object, object], dct_2: dict[object, object]
) -> dict[object, object]:
    """Given two dictionaries (dct1 bigger dct2), filter out
    from 1 all the entries present in 2.

    :param dct1: First dictionary
    :param dct2: Second dictionary
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


def merge_sequence(dcts: Sequence[dict[object, object]]):
    """Merge a sequence of dictionaries.
    :param dcts: Dictionary to merge
    :return: Merged dictionaries.
    """
    merged_dct = {}
    for dct in dcts:
        merged_dct.update(dct)
    return merged_dct


def sort_value_(
    dct: dict[object, object],
    allow_missing: bool = True,
    missing_val: Any = None,
):
    """Generate a sort value function from a dictionary.

    :param dct: A dictionary
    :param allow_missing: Allow missing values?, defaults to True
    :return:
    """

    def sort_value(key: object):
        assert (
            allow_missing or key in dct
        ), "No key {key} in dictionary:\n{dict[object,object]}"
        return dct[key] if key in dct else missing_val

    return sort_value
