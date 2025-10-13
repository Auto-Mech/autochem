"""Functions for dealing with a list of items encoding a ring."""

import itertools
from collections.abc import Sequence

import more_itertools as mit


def normalize(items: Sequence[object]) -> tuple[object, ...]:
    """Normalize the ordering of items in a ring.

    :param items: The ring items
    :return: The ring items, in normalized order
    """
    nitems = len(items)
    start_item = min(items)
    next_item = items[(items.index(start_item) + 1) % nitems]
    prev_item = items[(items.index(start_item) - 1) % nitems]
    if next_item > prev_item:
        items = list(reversed(items))
    return cycle(items, count=items.index(start_item))


def cycle(items: Sequence[object], count: int = 1) -> tuple[object]:
    """Cycle ring items once.

    Example:    (1, 2, 3, 4) -> (2, 3, 4, 1)

    :param items: The ring items
    :return: The ring items, cycled once
    """
    nitems = len(items)
    cycler = itertools.cycle(items)
    mit.consume(cycler, n=count % nitems)
    return tuple(itertools.islice(cycler, nitems))


def edges(items: Sequence[object]) -> tuple[tuple[object, object], ...]:
    """Get the edge pairs of a ring.

    Example:    (1, 2, 3, 4) -> ((1, 2), (2, 3), (3, 4), (4, 1))

    :param items: The ring items
    :return: The ring edge pairs
    """
    return tuple(zip(items, cycle(items), strict=False))


def distance(
    items: Sequence[object], item1: object, item2: object, longest: bool = False
) -> int:
    """Find the distance between two items in a ring.

    By default, finds the shortest distance. Setting `longest=True` results in the
    longest distance.

    :param items: The ring items
    :param item1: The item to start measuring distance from
    :param item2: The item to measure distance to
    :param longest: Return the longest distance?, defaults to False
    :return: The number of ring items between these two
    """
    assert (
        item1 in items and item2 in items
    ), f"Items {item1} and {item2} must both be in the ring {items}"

    # If they are the same item and we are looking for the shortest distance, return 0
    if item1 == item2 and not longest:
        return 0

    # To find the short and long distances, we find the positions of the items around
    # the ring: The first position of item1, the next position item2 after that, and the
    # next position of item1 after that
    enum_cycler = enumerate(itertools.cycle(items))
    pos1 = next(ix for ix, it in enum_cycler if it == item1)
    pos2 = next(ix for ix, it in enum_cycler if it == item2)
    pos3 = next(ix for ix, it in enum_cycler if it == item1)

    # The short and long distances are the differences in position
    dists = [abs(pos2 - pos1), abs(pos3 - pos2)]

    return min(dists) if not longest else max(dists)


def cycle_item_to_front(
    items: Sequence[object], item: object, end_item: object = None
) -> tuple[object, ...]:
    """Cycle ring items until one is in front.

    Optionally, request one adjacent item to be at the end, reversing the ring order if
    necessary.

    :param items: The ring items
    :parm item: The item to cycle to the font
    :param end_item: optionally, ensure that this is the last item in the ring
    :returns: The ring items with `item` cycled to the front and `end_item` to the end
    """
    items = cycle_items_to_front(items, [item])

    if end_item is not None and items[-1] != end_item:
        items = cycle_items_to_front(list(reversed(items)), [item])

    assert (
        end_item is None or items[-1] == end_item
    ), f"Cannot have {item} at start and {end_item} at end of ring {items}"

    return items


def cycle_items_to_front(
    items: Sequence[object], front_items: Sequence[object]
) -> tuple[object, ...]:
    """Cycle ring items until a group of adjacent items is at the front of the list.

    :param items: The ring items
    :parm front_items: The item to cycle to the font
    :returns: The ring items with `item` cycled to the front and `end_item` to the end
    """
    nitems = len(items)

    cycler = itertools.cycle(items)

    # 1. Cycle through all of the items, to make sure we aren't splitting them up
    cycler = itertools.dropwhile(lambda x: x in front_items, cycler)

    # 2. Now, cycle until we hit any one of the items
    cycler = itertools.dropwhile(lambda x: x not in front_items, cycler)

    # 3. Now, return the current cycle
    return tuple(itertools.islice(cycler, nitems))


def cycle_items_to_back(
    items: Sequence[object], back_items: Sequence[object]
) -> tuple[object, ...]:
    """Cycle ring items until a group of adjacent items is at the end of the list.

    :param items: The ring items
    :parm back_items: The item to cycle to the end
    :returns: The ring items `back_items` cycled to the end
    """
    nitems = len(items)

    cycler = itertools.cycle(items)

    # 1. Cycle until we hit any one of the items
    cycler = itertools.dropwhile(lambda x: x not in back_items, cycler)

    # 2. Cycle through all of the items, to put them at the end
    cycler = itertools.dropwhile(lambda x: x in back_items, cycler)

    # 3. Now, return the current cycle
    return tuple(itertools.islice(cycler, nitems))


def cycle_to_split(
    items: Sequence[object], split_pair: tuple[object, object]
) -> tuple[object, ...]:
    """Cycle to split a pair of adjacent items, putting one on each end of the list.

    :param items: The ring items
    :type items: List[object]
    :param split_pair: The pair of items to split
    :type split_pair: Tuple[object, object]
    :return: The ring items with one item in the pair at the start and one at the end
    """
    nitems = len(items)
    split_pair = set(split_pair)

    # If it is already split at this point, return as-is
    if {items[0], items[-1]} == split_pair:
        return items

    # Otherwise, cycle until you hit either item in the pair
    cycler = itertools.cycle(items)
    cycler = itertools.dropwhile(lambda x: x not in split_pair, cycler)

    # Then go one further to split them
    next(cycler)

    items = tuple(itertools.islice(cycler, nitems))

    assert {
        items[0],
        items[-1],
    } == split_pair, f"Cannot split ring {items} on non-adjacent pair {split_pair}"

    # Now, return the current cycle
    return items


def cycle_to_optimal_split(
    items: Sequence[object],
    split_pairs: Sequence[tuple[object, object]],
    back_items: Sequence[object],
) -> tuple[object, ...]:
    """Cycle to find an "optimum split" that puts a subset of items as close as possible
    to the end of the list.

    :param items: The ring items
    :param split_pairs: Pairs where the cycle can be split
    :param back_items: Adjacent items that should be as close to the end as possible
    :return: The cycle with the optimal split
    """
    orig_items = items
    split_pairs = list(split_pairs)
    back_items = list(back_items)

    if not split_pairs and not back_items:
        return orig_items

    if not split_pairs:
        return cycle_items_to_back(orig_items, back_items)

    if not back_items:
        return cycle_to_split(orig_items, split_pairs[0])

    # Try each split in the forward and reverse direction and see which puts the
    # `back_items` closest to the end of the list
    opt_items = None
    max_dist = 0
    for split_pair in split_pairs:
        items = cycle_to_split(orig_items, split_pair)
        for items_ in [items, tuple(reversed(items))]:
            dist = next(ix for ix, it in enumerate(items_) if it in back_items)
            if dist > max_dist:
                opt_items = items_
                max_dist = dist

    return opt_items
