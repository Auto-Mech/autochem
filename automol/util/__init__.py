""" common utilities used by automol
"""

# functions
from automol.util._util import flatten
from automol.util._util import is_odd_permutation
from automol.util._util import is_even_permutation
from automol.util._util import equivalence_partition
from automol.util._util import move_item_to_front
from automol.util._util import move_item_to_end
from automol.util._util import move_items_to_front
from automol.util._util import breakby
from automol.util._util import separate_negatives
from automol.util._util import value_similar_to
from automol.util._util import scale_iterable
from automol.util._util import remove_duplicates_with_order
from automol.util._util import sort_by_list
from automol.util._util import formula_from_symbols
from automol.util._util import numpy_to_float
# submodules
from automol.util import vec
from automol.util import mat
from automol.util import highd_mat
from automol.util import dict_
from automol.util import heuristic
from automol.util import dummy_trans


__all__ = [
    # functions
    'flatten',
    'is_odd_permutation',
    'is_even_permutation',
    'equivalence_partition',
    'move_item_to_front',
    'move_item_to_end',
    'move_items_to_front',
    'breakby',
    'separate_negatives',
    'value_similar_to',
    'scale_iterable',
    'remove_duplicates_with_order',
    'sort_by_list',
    'formula_from_symbols',
    'numpy_to_float',
    # submodules
    'vec',
    'mat',
    'highd_mat',
    'dict_',
    'heuristic',
    'dummy_trans',
]
