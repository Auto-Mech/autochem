""" common utilities used by automol
"""

# functions
from automol.util._util import is_even_permutation
from automol.util._util import equivalence_partition
from automol.util._util import separate_negatives
from automol.util._util import value_similar_to
from automol.util._util import scale_iterable
from automol.util._util import remove_duplicates_with_order
from automol.util._util import formula_from_symbols
from automol.util._util import numpy_to_float
# submodules
from automol.util import vec
from automol.util import mat
from automol.util import highd_mat
from automol.util import dict_


__all__ = [
    # functions
    'is_even_permutation',
    'equivalence_partition',
    'separate_negatives',
    'value_similar_to',
    'scale_iterable',
    'remove_duplicates_with_order',
    'formula_from_symbols',
    'numpy_to_float',
    # submodules
    'vec',
    'mat',
    'highd_mat',
    'dict_',
]
