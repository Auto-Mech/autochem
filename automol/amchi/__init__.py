""" AMChI (AutoMech Chemical Identifier) strings

Closely follows the InChI format, although canonicalizations may differ.
Extends InChI to allow for resonance double-bond stereo.
"""

# L2
# # constructor
from automol.amchi.base._core import from_data
# # recalculate/standardize
from automol.amchi.base._core import standard_form
# # getters
from automol.amchi.base._core import prefix
from automol.amchi.base._core import version
from automol.amchi.base._core import formula_layer
from automol.amchi.base._core import main_layers
from automol.amchi.base._core import charge_layers
from automol.amchi.base._core import stereo_layers
from automol.amchi.base._core import isotope_layers
from automol.amchi.base._core import ts_layers
# # setters
from automol.amchi.base._core import with_inchi_prefix
from automol.amchi.base._core import reflect
from automol.amchi.base._core import canonical_enantiomer
from automol.amchi.base._core import reflect_reaction
from automol.amchi.base._core import canonical_enantiomer_reaction
# # conversions
from automol.amchi.base._core import formula
from automol.amchi.base._core import formula_string
from automol.amchi.base._core import connectivity
from automol.amchi.base._core import without_stereo
from automol.amchi.base._core import racemic
from automol.amchi.base._core import are_enantiomers
from automol.amchi.base._core import are_diastereomers
# # properties
# # # formula layer
from automol.amchi.base._core import symbols
from automol.amchi.base._core import canonical_indices
# # # main layers
from automol.amchi.base._core import bonds
from automol.amchi.base._core import hydrogen_valences
# # # charge layers
from automol.amchi.base._core import charge
# # # stereo layers
from automol.amchi.base._core import bond_stereo_parities
from automol.amchi.base._core import atom_stereo_parities
from automol.amchi.base._core import is_inverted_enantiomer
from automol.amchi.base._core import is_canonical_enantiomer
from automol.amchi.base._core import is_canonical_enantiomer_reaction
from automol.amchi.base._core import is_canonical_reaction_direction
from automol.amchi.base._core import is_enantiomer_list
from automol.amchi.base._core import is_enantiomer_reaction
# # # isotope layers
from automol.amchi.base._core import bond_isotope_stereo_parities
from automol.amchi.base._core import atom_isotope_stereo_parities
from automol.amchi.base._core import is_inverted_isotope_enantiomer
# # other properties
from automol.amchi.base._core import has_multiple_components
from automol.amchi.base._core import has_stereo
from automol.amchi.base._core import has_mobile_hydrogens
from automol.amchi.base._core import low_spin_multiplicity
from automol.amchi.base._core import is_enantiomer
# # comparisons
from automol.amchi.base._core import same_connectivity
from automol.amchi.base._core import equivalent
# # split/join
from automol.amchi.base._core import split
from automol.amchi.base._core import join
# # sort
from automol.amchi.base._core import sorted_
from automol.amchi.base._core import argsort
# # helpers
from automol.amchi.base._core import join_layers
from automol.amchi.base._core import split_layers
from automol.amchi.base._core import join_layer_strings
from automol.amchi.base._core import _split_layer_string
# L4
# # conversions
from automol.amchi._conv import amchi_key
from automol.amchi._conv import smiles
from automol.amchi._conv import graph
from automol.amchi._conv import geometry
from automol.amchi._conv import conformers
from automol.amchi._conv import zmatrix
from automol.amchi._conv import rdkit_molecule
from automol.amchi._conv import rdkit_reaction
from automol.amchi._conv import display
from automol.amchi._conv import display_reaction
# # derived properties
from automol.amchi._conv import is_complete
# # derived transformations
from automol.amchi._conv import add_stereo
from automol.amchi._conv import expand_stereo
# drawing tools
from automol.amchi._draw import draw
from automol.amchi._draw import draw_grid
# assessment tools
from automol.amchi._assess import is_valid_multiplicity

# allow this as an alias
is_chiral = is_enantiomer


__all__ = [
    # L2
    # # constructor
    'from_data',
    # # recalculate/standardize
    'standard_form',
    # # getters
    'prefix',
    'version',
    'formula_layer',
    'main_layers',
    'charge_layers',
    'stereo_layers',
    'isotope_layers',
    'ts_layers',
    # # setters
    'with_inchi_prefix',
    'reflect',
    'canonical_enantiomer',
    'reflect_reaction',
    'canonical_enantiomer_reaction',
    # # conversions
    'formula',
    'formula_string',
    'connectivity',
    'without_stereo',
    'racemic',
    'are_enantiomers',
    'are_diastereomers',
    # # properties
    # # # formula layer
    'symbols',
    'canonical_indices',
    # # # main layers
    'bonds',
    'hydrogen_valences',
    # # # charge layers
    'charge',
    # # # stereo layers
    'bond_stereo_parities',
    'atom_stereo_parities',
    'is_inverted_enantiomer',
    'is_canonical_enantiomer',
    'is_canonical_enantiomer_reaction',
    'is_canonical_reaction_direction',
    'is_enantiomer_list',
    'is_enantiomer_reaction',
    # # # isotope layers
    'bond_isotope_stereo_parities',
    'atom_isotope_stereo_parities',
    'is_inverted_isotope_enantiomer',
    # # other properties
    'has_multiple_components',
    'is_enantiomer',
    'is_chiral',
    'has_stereo',
    'has_mobile_hydrogens',
    'low_spin_multiplicity',
    # # comparisons
    'same_connectivity',
    'equivalent',
    # # split/join
    'split',
    'join',
    # # sort
    'sorted_',
    'argsort',
    'join_layers',
    'split_layers',
    'join_layer_strings',
    '_split_layer_string',
    # L4
    # # conversions
    'amchi_key',
    'smiles',
    'graph',
    'geometry',
    'conformers',
    'zmatrix',
    'rdkit_molecule',
    'rdkit_reaction',
    'display',
    'display_reaction',
    # # derived properties
    'is_complete',
    # # derived transformations
    'add_stereo',
    'expand_stereo',
    # drawing tools
    'draw',
    'draw_grid',
    # assessment tools
    'is_valid_multiplicity',
]
