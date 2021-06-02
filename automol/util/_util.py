""" miscellaneous utilities
"""
from phydat import ptab


def is_even_permutation(seq1, seq2):
    """ Determine whether a permutation of a sequence is even or odd.

    :param seq1: the first sequence
    :param seq2: the second sequence, which must be a permuation of the first
    :returns: True if the permutation is even, False if it is odd
    :rtype: bool
    """
    size = len(seq1)
    assert sorted(seq1) == sorted(seq2) and len(set(seq1)) == size
    perm = [seq2.index(val) for val in seq1]

    sgn = 1
    for idx in range(size):
        if perm[idx] != idx:
            sgn *= -1
            swap_idx = perm.index(idx)
            perm[idx], perm[swap_idx] = perm[swap_idx], perm[idx]

    parity = (sgn == 1)

    return parity


def equivalence_partition(iterable, relation):
    """Partitions a set of objects into equivalence classes

    canned function taken from https://stackoverflow.com/a/38924631

    Args:
        iterable: collection of objects to be partitioned
        relation: equivalence relation. I.e. relation(o1,o2) evaluates to True
            if and only if o1 and o2 are equivalent

    Returns: classes, partitions
        classes: A sequence of sets. Each one is an equivalence class
    """
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
    return classes


# Useful functions on Python objects
def separate_negatives(lst):
    """ Seperate a list of numbers into negative and nonnegative (>= 0)
    """

    neg_lst = tuple(val for val in lst if val < 0)
    pos_lst = tuple(val for val in lst if val >= 0)

    return neg_lst, pos_lst


def value_similar_to(val, lst, thresh):
    """ Check if a value is close to some lst of values within some threshold
    """
    return all(abs(val - vali) < thresh for vali in lst)


def scale_iterable(iterable, scale_factor):
    """ Scale some type of iterable of floats by a scale factor
    """

    if isinstance(iterable, list):
        scaled_iterable = list(val * scale_factor for val in iterable)
    elif isinstance(iterable, tuple):
        scaled_iterable = tuple(val * scale_factor for val in iterable)

    return scaled_iterable


def formula_from_symbols(symbs):
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
