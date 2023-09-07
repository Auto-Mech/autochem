""" common utilities used by automol
"""

from automol.util import dict_, dummy_conv, heuristic, highd_mat, mat, ring, vec
from automol.util._util import (
    breakby,
    equivalence_partition,
    flatten,
    formula_from_symbols,
    is_even_permutation,
    is_odd_permutation,
    move_item_to_end,
    move_item_to_front,
    move_items_to_front,
    numpy_to_float,
    remove_duplicates_with_order,
    scale_iterable,
    separate_negatives,
    sort_by_list,
    value_similar_to,
)
from automol.util.dummy_conv import DummyConv

__all__ = [
    "dict_",
    "dummy_conv",
    "heuristic",
    "highd_mat",
    "mat",
    "ring",
    "vec",
    "breakby",
    "equivalence_partition",
    "flatten",
    "formula_from_symbols",
    "is_even_permutation",
    "is_odd_permutation",
    "move_item_to_end",
    "move_item_to_front",
    "move_items_to_front",
    "numpy_to_float",
    "remove_duplicates_with_order",
    "scale_iterable",
    "separate_negatives",
    "sort_by_list",
    "value_similar_to",
    "DummyConv",
]
