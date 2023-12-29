""" miscellaneous utilities
"""
import itertools
from collections.abc import Iterable
from typing import List

from phydat import ptab


def flatten(lst):
    """Flatten an arbitrarily nested list of lists (iterator)

    Source: https://stackoverflow.com/a/2158532
    """
    for elem in lst:
        if isinstance(elem, Iterable) and not isinstance(elem, (str, bytes)):
            yield from flatten(elem)
        else:
            yield elem


def is_odd_permutation(seq1: List, seq2: List):
    """Determine whether a permutation of a sequence is odd.

    :param seq1: the first sequence
    :param seq2: the second sequence, which must be a permuation of the first
    :returns: True if the permutation is even, False if it is odd
    :rtype: bool
    """
    return not is_even_permutation(seq1, seq2)


def is_even_permutation(seq1: List, seq2: List):
    """Determine whether a permutation of a sequence is even or odd.

    :param seq1: the first sequence
    :param seq2: the second sequence, which must be a permuation of the first
    :returns: True if the permutation is even, False if it is odd
    :rtype: bool
    """
    size = len(seq1)
    assert set(seq1) == set(seq2) and len(set(seq1)) == size
    perm = [seq2.index(val) for val in seq1]

    sgn = 1
    for idx in range(size):
        if perm[idx] != idx:
            sgn *= -1
            swap_idx = perm.index(idx)
            perm[idx], perm[swap_idx] = perm[swap_idx], perm[idx]

    parity = sgn == 1

    return parity


def equivalence_partition(iterable, relation, perfect=False):
    """Partitions a set of objects into equivalence classes

    canned function taken from https://stackoverflow.com/a/38924631

    Args:
        iterable: collection of objects to be partitioned
        relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            if and only if o1 and o2 are equivalent
        perfect: is this a perfect equivalence relation, where a = c and b = c
            guarantees a = b? if not, an extra search is performed to make sure
            that a, b, and c still end up in the same class

    Returns: classes, partitions
        classes: A sequence of sets. Each one is an equivalence class
    """
    # 1. This part only works assuming it is a 'perfect' equivalence relation,
    # where a = c and b = c implies a = b
    classes = []
    for obj in iterable:  # for each object
        # find the class it is in
        found = False
        for cls in classes:
            # is it equivalent to this class?
            if relation(next(iter(cls)), obj):
                cls.add(obj)
                found = True
                break
        if not found:  # it is in a new class
            classes.append(set([obj]))

    # 2. Now, account for the possibility of 'imperfect' equivalence relations,
    # where the relation gives a = c and b = c, but not a = b, and yet we still
    # want a, b, and c to end up in the same class
    if not perfect:
        new_classes = []
        while True:
            new_classes = classes.copy()
            for cls1, cls2 in itertools.combinations(classes, r=2):
                if any(relation(o1, o2) for o1, o2 in itertools.product(cls1, cls2)):
                    if cls2 in new_classes:
                        new_classes.remove(cls2)
                        cls1 |= cls2

            if classes == new_classes:
                break

            classes = new_classes

    return classes


# Useful functions on Python objects
def move_item_to_front(lst, item):
    """Move an item to the front of a list.

    :param lst: the list
    :type lst: list or tuple
    :param item: the item, which must be in `lst`
    :returns: the list, with the item moved to front
    :rtype: tuple
    """
    lst = list(lst)
    lst.insert(0, lst.pop(lst.index(item)))
    return tuple(lst)


def move_item_to_end(lst, item):
    """Move an item to the end of a list.

    :param lst: the list
    :type lst: list or tuple
    :param item: the item, which must be in `lst`
    :returns: the list, with the item moved to end
    :rtype: tuple
    """
    lst = list(lst)
    lst.append(lst.pop(lst.index(item)))
    return tuple(lst)


def move_items_to_front(lst, items):
    """Move an item to the front of a list.

    :param lst: the list
    :type lst: list or tuple
    :param item: the item, which must be in `lst`
    :returns: the list, with the item moved to front
    :rtype: tuple
    """
    lst = list(lst)
    for item in reversed(items):
        lst.insert(0, lst.pop(lst.index(item)))
    return tuple(lst)


def breakby(lst, elem):
    """Break a list by element, dropping the element itself.

    Analogous to '<char>'.split('<string>') for strings.
    """
    lsts = tuple(
        tuple(g) for k, g in itertools.groupby(lst, lambda x: x == elem) if not k
    )
    return lsts


def separate_negatives(lst):
    """Seperate a list of numbers into negative and nonnegative (>= 0)"""

    neg_lst = tuple(val for val in lst if val < 0)
    pos_lst = tuple(val for val in lst if val >= 0)

    return neg_lst, pos_lst


def value_similar_to(val, lst, thresh):
    """Check if a value is close to some lst of values within some threshold"""
    return any(abs(val - vali) < thresh for vali in lst)


def scale_iterable(iterable, scale_factor):
    """Scale some type of iterable of floats by a scale factor"""

    if isinstance(iterable, list):
        scaled_iterable = list(val * scale_factor for val in iterable)
    elif isinstance(iterable, tuple):
        scaled_iterable = tuple(val * scale_factor for val in iterable)

    return scaled_iterable


def remove_duplicates_with_order(lst):
    """Remove all duplicates of a list while not reordering the list."""
    if isinstance(lst, list):
        lst = list(n for i, n in enumerate(lst) if n not in lst[:i])
    if isinstance(lst, tuple):
        lst = tuple(n for i, n in enumerate(lst) if n not in lst[:i])

    return lst


def sort_by_list(lst, ref_lst, include_missing=True):
    """Order the elements of the list by using the priorities given
    by some reference lst.

    if include_missing:
    a=[q, a, e, x, f, t], ref=[x, a, q, e] -> sort_a=[x, a, q, e, f, t]
    if not include_missing:
    a=[q, a, e, x, f], ref=[x, a, q, e] -> sort_a=[x, a, q, e]

    Note that any element in the original list not in original list is
    dropped if the user specifies not to include it.

    :param lst: list to sort
    :type lst: tuple
    :param ref_lst: list which sets the order of the previous list
    :type ref_lst: tuple
    :rtype: tuple
    """

    # Split input list by elements in and not in reference list
    x_in_ref = tuple(x for x in lst if x in ref_lst)
    x_missing = tuple(x for x in lst if x not in ref_lst)

    # Sorted list of elements in th reference
    sort_lst = tuple(sorted(list(x_in_ref), key=lambda x: ref_lst.index(x)))

    # If request append the missing elements
    if include_missing:
        sort_lst += x_missing

    return sort_lst


def formula_from_symbols(symbs):
    """Build a molecular formula from a list of atomic symbols.

    (note: dummy atoms will be filtered out and cases will be standardized)

    :param symbs: atomic symbols
    :type symbs: tuple(str)
    :rtype: str
    """

    symbs = list(filter(ptab.to_number, map(ptab.to_symbol, symbs)))

    return _unique_item_counts(symbs)


def _unique_item_counts(iterable):
    """Build a dictionary giving the count of each unique item in a sequence.

    :param iterable: sequence to obtain counts for
    :type iterable: iterable object
    :rtype: dict[obj: int]
    """

    items = tuple(iterable)

    return {item: items.count(item) for item in sorted(set(items))}
