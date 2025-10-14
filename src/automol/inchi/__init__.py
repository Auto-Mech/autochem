""" InChI strings
"""

# L3
# # "constructor"
from .base._core import from_data
# # recalculate/standardize
from .base._core import recalculate
from .base._core import standard_form
# # getters
from .base._core import version
from .base._core import formula_layer
from .base._core import main_layers
from .base._core import charge_layers
from .base._core import stereo_layers
from .base._core import isotope_layers
from .base._core import stereo_atoms
from .base._core import stereo_bonds
from .base._core import unassigned_stereo_bonds
from .base._core import is_enantiomer
from .base._core import are_enantiomers
from .base._core import are_diastereomers
from .base._core import reflect
# # conversions
from .base._core import inchi_key
from .base._core import smiles
from .base._core import formula
from .base._core import formula_string
from .base._core import connectivity
from .base._core import without_stereo
from .base._core import racemic
# # properties
from .base._core import is_standard_form
from .base._core import has_multiple_components
from .base._core import has_stereo
from .base._core import low_spin_multiplicity
# # comparisons
from .base._core import same_connectivity
from .base._core import equivalent
# # sort
from .base._core import sorted_
from .base._core import argsort
# # split/join
from .base._core import split
from .base._core import join
# L4
# # conversions
from ._conv import graph
from ._conv import geometry
from ._conv import zmatrix
from ._conv import amchi
from ._conv import rdkit_molecule
from ._conv import rdkit_reaction
from ._conv import display
from ._conv import display_reaction
# # derived properties
from ._conv import is_complete
from ._conv import is_bad
# # derived transformations
from ._conv import add_stereo
from ._conv import expand_stereo
# drawing tools
from ._draw import draw
from ._draw import draw_grid
# assessment tools
from ._assess import is_valid_inchi_multiplicity

# allow this as an alias
is_chiral = is_enantiomer


__all__ = [
    # L3
    # # "constructor"
    'from_data',
    # # recalculate/standardize
    'recalculate',
    'standard_form',
    # # getters
    'version',
    'formula_layer',
    'main_layers',
    'charge_layers',
    'stereo_layers',
    'isotope_layers',
    'stereo_atoms',
    'stereo_bonds',
    'unassigned_stereo_bonds',
    'is_enantiomer',
    'is_chiral',
    'are_enantiomers',
    'are_diastereomers',
    'reflect',
    # # conversions
    'inchi_key',
    'smiles',
    'formula',
    'formula_string',
    'connectivity',
    'without_stereo',
    'racemic',
    # # properties
    'is_standard_form',
    'has_multiple_components',
    'has_stereo',
    'low_spin_multiplicity',
    # # comparisons
    'same_connectivity',
    'equivalent',
    # # sort
    'sorted_',
    'argsort',
    # # split/join
    'split',
    'join',
    # L4
    # # conversions
    'graph',
    'geometry',
    'zmatrix',
    'amchi',
    'rdkit_molecule',
    'rdkit_reaction',
    'display',
    'display_reaction',
    # # derived properties
    'is_complete',
    'is_bad',
    # # derived transformations
    'add_stereo',
    'expand_stereo',
    # # ddrawing tools
    'draw',
    'draw_grid',
    # assessment tools
    'is_valid_inchi_multiplicity',
]
