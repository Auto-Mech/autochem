""" InChI strings
"""

# L3
# # "constructor"
from automol.inchi.base._core import from_data
# # recalculate/standardize
from automol.inchi.base._core import recalculate
from automol.inchi.base._core import standard_form
# # getters
from automol.inchi.base._core import version
from automol.inchi.base._core import formula_layer
from automol.inchi.base._core import main_layers
from automol.inchi.base._core import charge_layers
from automol.inchi.base._core import stereo_layers
from automol.inchi.base._core import isotope_layers
from automol.inchi.base._core import stereo_atoms
from automol.inchi.base._core import stereo_bonds
from automol.inchi.base._core import unassigned_stereo_bonds
from automol.inchi.base._core import is_enantiomer
from automol.inchi.base._core import are_enantiomers
from automol.inchi.base._core import are_diastereomers
from automol.inchi.base._core import reflect
# # conversions
from automol.inchi.base._core import inchi_key
from automol.inchi.base._core import smiles
from automol.inchi.base._core import formula
from automol.inchi.base._core import formula_string
from automol.inchi.base._core import connectivity
from automol.inchi.base._core import without_stereo
from automol.inchi.base._core import racemic
# # properties
from automol.inchi.base._core import is_standard_form
from automol.inchi.base._core import has_multiple_components
from automol.inchi.base._core import has_stereo
from automol.inchi.base._core import low_spin_multiplicity
# # comparisons
from automol.inchi.base._core import same_connectivity
from automol.inchi.base._core import equivalent
# # sort
from automol.inchi.base._core import sorted_
from automol.inchi.base._core import argsort
# # split/join
from automol.inchi.base._core import split
from automol.inchi.base._core import join
# L4
# # conversions
from automol.inchi._conv import graph
from automol.inchi._conv import geometry
from automol.inchi._conv import conformers
from automol.inchi._conv import zmatrix
from automol.inchi._conv import amchi
from automol.inchi._conv import rdkit_molecule
from automol.inchi._conv import rdkit_reaction
from automol.inchi._conv import display
from automol.inchi._conv import display_reaction
# # derived properties
from automol.inchi._conv import is_complete
from automol.inchi._conv import is_bad
# # derived transformations
from automol.inchi._conv import add_stereo
from automol.inchi._conv import expand_stereo
# drawing tools
from automol.inchi._draw import draw
from automol.inchi._draw import draw_grid
# assessment tools
from automol.inchi._assess import is_valid_inchi_multiplicity

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
    'conformers',
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
