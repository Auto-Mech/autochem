""" miscellaneous conversion utilities
"""

from phydat import ptab


def formula(symbs):
    """ molecular formula from a list of atomic symbols

    (note: dummy atoms will be filtered out and cases will be standardized)
    """
    symbs = list(filter(ptab.to_number, map(pt.to_symbol, symbs)))
    return _unique_item_counts(symbs)


def _unique_item_counts(iterable):
    """ a dictionary giving the count of each unique item in a sequence
    """
    items = tuple(iterable)
    return {item: items.count(item) for item in sorted(set(items))}
