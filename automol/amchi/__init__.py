""" AMChI (AutoMech Chemical Identifier) strings

Closely follows the InChI format, although canonicalizations may differ.
Extends InChI to allow for resonance double-bond stereo.
"""

# L2
# # constructor
from .base._core import from_data
# # recalculate/standardize
from .base._core import standard_form
# # getters
from .base._core import prefix
from .base._core import version
from .base._core import formula_layer
from .base._core import main_layers
from .base._core import charge_layers
from .base._core import stereo_layers
from .base._core import isotope_layers
from .base._core import ts_layers
# # setters
from .base._core import with_inchi_prefix
from .base._core import reflect
from .base._core import canonical_enantiomer
from .base._core import reflect_reaction
from .base._core import canonical_enantiomer_reaction
# # conversions
from .base._core import formula
from .base._core import formula_string
from .base._core import connectivity
from .base._core import without_stereo
from .base._core import racemic
from .base._core import are_enantiomers
from .base._core import are_diastereomers
# # properties
# # # prefix
from .base._core import is_amchi
from .base._core import is_inchi
# # # formula layer
from .base._core import symbols
from .base._core import canonical_indices
# # # main layers
from .base._core import bonds
from .base._core import hydrogen_valences
# # # charge layers
from .base._core import charge
# # # stereo layers
from .base._core import bond_stereo_parities
from .base._core import atom_stereo_parities
from .base._core import is_inverted_enantiomer
from .base._core import is_canonical_enantiomer
from .base._core import is_canonical_enantiomer_reaction
from .base._core import is_canonical_reaction_direction
from .base._core import is_enantiomer_list
from .base._core import is_enantiomer_reaction
# # # TS layers
from .base._core import breaking_bond_keys
from .base._core import forming_bond_keys
from .base._core import is_reversed_ts
# # # isotope layers
from .base._core import bond_isotope_stereo_parities
from .base._core import atom_isotope_stereo_parities
from .base._core import is_inverted_isotope_enantiomer
# # other properties
from .base._core import has_multiple_components
from .base._core import has_stereo
from .base._core import has_mobile_hydrogens
from .base._core import low_spin_multiplicity
from .base._core import is_enantiomer
from .base._core import is_racemic
# # comparisons
from .base._core import same_connectivity
from .base._core import equivalent
# # split/join
from .base._core import split
from .base._core import join
from .base._core import sorted_join
# # sort
from .base._core import sorted_
from .base._core import argsort
# # helpers
from .base._core import join_layers
from .base._core import split_layers
from .base._core import join_layer_strings
from .base._core import _split_layer_string
# L4
# # conversions
from ._conv import amchi_key
from ._conv import chemkin_name
from ._conv import connectivity_digest
from ._conv import stereo_digest
from ._conv import chi_
from ._conv import smiles
from ._conv import graph
from ._conv import geometry
from ._conv import zmatrix
from ._conv import rdkit_molecule
from ._conv import rdkit_reaction
from ._conv import display
from ._conv import display_reaction
# # derived properties
from ._conv import is_complete
from ._conv import is_valid_multiplicity
from ._conv import guess_spin
# # derived transformations
from ._conv import add_stereo
from ._conv import expand_stereo
# drawing tools
from ._draw import draw
from ._draw import draw_grid

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
    # # # prefix
    'is_amchi',
    'is_inchi',
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
    # # # TS layers
    'breaking_bond_keys',
    'forming_bond_keys',
    'is_reversed_ts',
    # # # isotope layers
    'bond_isotope_stereo_parities',
    'atom_isotope_stereo_parities',
    'is_inverted_isotope_enantiomer',
    # # other properties
    'has_multiple_components',
    'is_enantiomer',
    'is_racemic',
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
    'sorted_join',
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
    'chemkin_name',
    'connectivity_digest',
    'stereo_digest',
    'chi_',
    'smiles',
    'graph',
    'geometry',
    'zmatrix',
    'rdkit_molecule',
    'rdkit_reaction',
    'display',
    'display_reaction',
    # # derived properties
    'is_complete',
    'is_valid_multiplicity',
    'guess_spin',
    # # derived transformations
    'add_stereo',
    'expand_stereo',
    # drawing tools
    'draw',
    'draw_grid',
]
