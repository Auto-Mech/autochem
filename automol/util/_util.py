"""miscellaneous utilities."""

import itertools
from collections.abc import Callable, Collection, Iterable, Iterator, Sequence
from numbers import Number
from typing import Any

from phydat import ptab


def partner(pair: Collection, item: Any) -> Any:
    """Get the partner of an item in a pair.

    The two items must be distinct

    :param pair: An iterable of length 2
    :param item: One of two items
    :return: The other item
    """
    pair = set(pair)
    assert len(pair) == 2 and item in pair
    return next(iter(pair - {item}))


def flatten(lst: Collection) -> Iterator:
    """Flatten an arbitrarily nested list of lists (iterator).
    Source: https://stackoverflow.com/a/2158532 .
    :param lst: An arbitrarily nested list or tuple
    :return: An iterator over the flattened list of values.
    """
    for elem in lst:
        if isinstance(elem, Iterable) and not isinstance(elem, str | bytes):
            yield from flatten(elem)
        else:
            yield elem


def translate(
    seq: Collection, trans_dct: dict, drop: bool = False, item_typ: type = Number
) -> Collection:
    """Translate items in a nested sequence or collection with a dictionary.

    :param seq: An arbitrarily nested sequence or collection
    :param trans_dct: A translation dictionary
    :param drop: Drop values missing from translation dictionary?, defaults to False
    :return: Translated version of collection
    """

    def transform_(seq_in: Collection) -> Collection:
        """Recursively convert a nested list of z-matrix keys to geometry keys."""
        assert isinstance(seq_in, Collection), f"Cannot process non-sequence {seq_in}"
        type_ = type(seq_in)

        seq_out = []
        for item in seq_in:
            if isinstance(item, Collection) and not isinstance(item, item_typ):
                seq_out.append(transform_(item))
            elif not drop or item in trans_dct:
                seq_out.append(trans_dct.get(item))

        return type_(seq_out)

    return transform_(seq)


def is_odd_permutation(seq1: list, seq2: list) -> bool:
    """Determine whether a permutation of a sequence is odd.

    :param seq1: The first sequence
    :param seq2: The second sequence, which must be a permuation of the first
    :returns: True if the permutation is even, False if it is odd
    """
    return not is_even_permutation(seq1, seq2)


def is_even_permutation(seq1: Sequence, seq2: Sequence, check: bool = True) -> bool:
    """Determine whether a permutation of a sequence is even.

    :param seq1: The first sequence
    :param seq2: The second sequence, which must be a permuation of the first
    :returns: True if the permutation is even, False if it is odd
    """
    size = len(seq1)
    if check:
        assert (
            set(seq1) == set(seq2) and len(set(seq1)) == size
        ), f"No permutation between sequences:\n{seq1}\n{seq2}"

    perm = [seq2.index(val) for val in seq1]

    sgn = 1
    for idx in range(size):
        if perm[idx] != idx:
            sgn *= -1
            swap_idx = perm.index(idx)
            perm[idx], perm[swap_idx] = perm[swap_idx], perm[idx]

    parity = sgn == 1

    return parity


def equivalence_partition(
    iterable: Collection, relation: Callable[[Any, Any], bool], perfect: bool = False
) -> list:
    """Partitions a set of objects into equivalence classes.

    canned function taken from https://stackoverflow.com/a/38924631

    :param iterable: Collection of objects to be partitioned
    :param relation: Equivalence relation. I.e. relation(o1,o2) evaluates to True
        if and only if o1 and o2 are equivalent
    :param perfect: Is this a perfect equivalence relation, where a = c and b = c
            guarantees a = b? if not, an extra search is performed to make sure
            that a, b, and c still end up in the same class

    :returns:A sequence of sets. Each one is an equivalence class
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
            classes.append({obj})

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
def move_item_to_front(lst: Sequence, item) -> tuple:
    """Move an item to the front of a list.

    :param lst: The list
    :param item: The item, which must be in `lst`
    :returns: The list, with the item moved to front
    """
    lst = list(lst)
    lst.insert(0, lst.pop(lst.index(item)))
    return tuple(lst)


def move_item_to_end(lst: Sequence, item) -> tuple:
    """Move an item to the end of a list.

    :param lst: The list
    :param item: The item, which must be in `lst`
    :returns: The list, with the item moved to end
    """
    lst = list(lst)
    lst.append(lst.pop(lst.index(item)))
    return tuple(lst)


def move_items_to_front(lst: Sequence, items) -> tuple:
    """Move an item to the front of a list.

    :param lst: The list
    :param item: The item, which must be in `lst`
    :returns: The list, with the item moved to front
    """
    lst = list(lst)
    for item in reversed(items):
        lst.insert(0, lst.pop(lst.index(item)))
    return tuple(lst)


def breakby(lst: Sequence, elem) -> tuple[tuple, ...]:
    """Break a list by element, dropping the element itself.
    Analogous to '<char>'.split('<string>') for strings.
    :param lst: The list
    :param elem: The element to break the list by, gets deleted
    :return:The chunks between the break points of the input list.
    """
    lsts = tuple(
        tuple(g) for k, g in itertools.groupby(lst, lambda x: x == elem) if not k
    )
    return lsts


def separate_negatives(lst: Sequence) -> tuple[tuple, tuple]:
    """Seperate a list of numbers into negative and nonnegative (>= 0).
    :param lst: The list
    :return: Value for negatives and for non-negatives.
    """
    neg_lst = tuple(val for val in lst if val < 0)
    pos_lst = tuple(val for val in lst if val >= 0)

    return neg_lst, pos_lst


def value_similar_to(val: float, lst: Sequence[float], thresh: float) -> bool:
    """Check if a value is close to any one of a list of values.
    :param val: A number.
    :param lst: A collection of numbers to compare to.
    :param thresh: The comparison threshold.
    :return: 'True' if close, 'False' if not.
    """
    return any(abs(val - vali) < thresh for vali in lst)


def scale_iterable(
    iterable: Collection[float], scale_factor: float
) -> Collection[float]:
    """Scale some type of iterable of floats by a scale factor.
    :param iterable: A list of numbers.
    :param scale_factor: A factor to scale by.
    :return: The scaled list of numbers.
    """
    if isinstance(iterable, list):
        iterable = [val * scale_factor for val in iterable]
    elif isinstance(iterable, tuple):
        iterable = tuple(val * scale_factor for val in iterable)

    return iterable


def remove_duplicates_with_order(lst: Sequence) -> Sequence:
    """Remove all duplicates of a list while not reordering the list.
    :param lst: A list.
    :return: A list, without duplicates.
    Note: To be deprecated
    and replaced with calls to more_itertools.unique_justseen.
    """
    if isinstance(lst, list):
        lst = [n for i, n in enumerate(lst) if n not in lst[:i]]
    if isinstance(lst, tuple):
        lst = tuple(n for i, n in enumerate(lst) if n not in lst[:i])

    return lst


def sort_by_list(lst: tuple, ref_lst: tuple, include_missing: bool = True) -> tuple:
    """Order the elements of the list by using the priorities given
    by some reference lst.

    if include_missing:
    a=[q, a, e, x, f, t], ref=[x, a, q, e] -> sort_a=[x, a, q, e, f, t]
    if not include_missing:
    a=[q, a, e, x, f], ref=[x, a, q, e] -> sort_a=[x, a, q, e]

    Note that any element in the original list not in original list is
    dropped if the user specifies not to include it.

    :param lst: List to sort
    :param ref_lst: List which sets the order of the previous list
    """
    # Split input list by elements in and not in reference list
    x_in_ref = tuple(x for x in lst if x in ref_lst)
    x_missing = tuple(x for x in lst if x not in ref_lst)

    # Sorted list of elements in th reference
    sort_lst = tuple(sorted(x_in_ref, key=lambda x: ref_lst.index(x)))

    # If request append the missing elements
    if include_missing:
        sort_lst += x_missing

    return sort_lst


def formula_from_symbols(symbs: tuple[str]) -> str:
    """Build a molecular formula from a list of atomic symbols.

    (note: dummy atoms will be filtered out and cases will be standardized)

    :param symbs: Atomic symbols
    """
    symbs = list(filter(ptab.to_number, map(ptab.to_symbol, symbs)))

    return _unique_item_counts(symbs)


def _unique_item_counts(iterable: Iterable) -> dict[object:int]:
    """Build a dictionary giving the count of each unique item in a sequence.

    :param iterable: Sequence to obtain counts for
    :type iterable: Iterable object
    """
    items = tuple(iterable)

    return {item: items.count(item) for item in sorted(set(items))}
