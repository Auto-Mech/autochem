""" miscellaneous conversion utilities
"""
from qcelemental import periodictable as pt


def formula(syms):
    """ molecular formula from a list of atomic symbols

    (note: dummy atoms will be filtered out and cases will be standardized)
    """
    syms = list(filter(pt.to_Z, map(pt.to_E, syms)))
    return _unique_item_counts(syms)


def _unique_item_counts(iterable):
    """ a dictionary giving the count of each unique item in a sequence
    """
    items = tuple(iterable)
    return {item: items.count(item) for item in sorted(set(items))}
