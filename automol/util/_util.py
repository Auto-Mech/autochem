""" miscellaneous utilities
"""


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
