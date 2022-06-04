""" Graph-specific utilities
"""
import numbers


def sort_keys(keys, func=lambda x: x):
    """ Sort a mixed list of atom and bond keys

        :param keys: a mixed list of atom and/or bond keys
        :param func: the sorting function to use
        :returns: the sorted list of keys
    """
    def _func(key):
        if isinstance(key, numbers.Number):
            ret = [func(key)]
        else:
            ret = sorted(map(func, key))
        return ret

    return tuple(sorted(keys, key=_func))


def sort_keys_representation(keys, func=lambda x: x):
    """ Generate the sortable representation used by the sort_keys function

        :param keys: a mixed list of atom and/or bond keys
        :param func: the transforming function to use
        :returns: the transformed list of keys
    """
    def _func(key):
        if isinstance(key, numbers.Number):
            ret = [func(key)]
        else:
            ret = sorted(map(func, key))
        return ret

    return tuple(map(_func, keys))
