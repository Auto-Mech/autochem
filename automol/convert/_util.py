""" miscellaneous conversion utilities
"""

from phydat import ptab


def formula(symbs):
    """ Build a molecular formula from a list of atomic symbols.

        (note: dummy atoms will be filtered out and cases will be standardized)

        :param symbs: atomic symbols
        :type symbs: tuple(str)
        :rtype: str
    """

    symbs = list(filter(ptab.to_number, map(ptab.to_symbol, symbs)))

    return _unique_item_counts(symbs)


def _unique_item_counts(iterable):
    """ Build a dictionary giving the count of each unique item in a sequence.

        :param iterable: sequence to obtain counts for
        :type iterable: iterable object
        :rtype: dict[obj: int]
    """

    items = tuple(iterable)

    return {item: items.count(item) for item in sorted(set(items))}
