""" Level 3 InChI functions (dependencies on extern, but not other types)
"""

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
from automol.inchi.base._core import without_stereo
from automol.inchi.base._core import racemic
from automol.inchi.base._core import connectivity
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


__all__ = [
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
    'are_enantiomers',
    'are_diastereomers',
    'reflect',
    # # conversions
    'inchi_key',
    'smiles',
    'formula',
    'formula_string',
    'without_stereo',
    'racemic',
    'connectivity',
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
]
