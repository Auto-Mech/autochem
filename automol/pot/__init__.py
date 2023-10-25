"""
  Libraries for constructing and manipulating potentials
"""

from automol.pot._pot import coords
from automol.pot._pot import scale
from automol.pot._pot import relax_scale
from automol.pot._pot import remove_empty_terms
from automol.pot._pot import by_index
from automol.pot._pot import is_nonempty
from automol.pot._read import find_max1d
from automol.pot._fit import setup_1d_potential


__all__ = [
    'coords',
    'scale',
    'relax_scale',
    'remove_empty_terms',
    'by_index',
    'is_nonempty',
    'find_max1d',
    'setup_1d_potential',
]
