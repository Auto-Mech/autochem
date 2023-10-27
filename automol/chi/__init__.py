""" ChI strings

A wrapper to redirect to either InChI or AMChI, as needed.

Currently only redirecting for resonance double bond stereo.
"""

# L3
# # "constructor"
from automol.chi.base._core import from_data
# # recalculate/standardize
from automol.chi.base._core import recalculate
from automol.chi.base._core import standard_form
# # getters
from automol.amchi.base import version
# from automol.chi.base._core import formula_sublayer  # defined below
from automol.amchi.base import formula_layer
from automol.amchi.base import main_layers as main_sublayers
from automol.amchi.base import charge_layers as charge_sublayers
from automol.amchi.base import stereo_layers as stereo_sublayers
from automol.amchi.base import isotope_layers as isotope_sublayers
# from automol.chi.base._core import stereo_atoms               # excluded
# from automol.chi.base._core import stereo_bonds               # excluded
# from automol.chi.base._core import unassigned_stereo_bonds    # excluded
from automol.chi.base._core import are_enantiomers
from automol.chi.base._core import are_diastereomers
from automol.chi.base._core import reflect
from automol.chi.base._core import racemic
# # setters
from automol.amchi.base._core import with_inchi_prefix
# # conversions
# from automol.chi.base._core import inchi_key  # below in L4
# from automol.chi.base._core import smiles     # below in L4
from automol.amchi.base import formula
from automol.amchi.base import formula_string
from automol.amchi.base import connectivity
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
# # # isotope layers
from automol.amchi.base._core import bond_isotope_stereo_parities
from automol.amchi.base._core import atom_isotope_stereo_parities
from automol.amchi.base._core import is_inverted_isotope_enantiomer
# # # other properties
from automol.chi.base._core import is_standard_form
from automol.amchi.base import has_multiple_components
from automol.amchi.base import is_enantiomer
from automol.amchi.base import is_enantiomer_list
from automol.amchi.base import is_enantiomer_reaction
from automol.amchi.base import has_stereo
from automol.amchi.base import has_mobile_hydrogens
from automol.amchi.base import low_spin_multiplicity
# # comparisons
from automol.amchi.base import same_connectivity
from automol.amchi.base import equivalent
# # split/join
from automol.chi.base._core import split
# from automol.chi.base._core import split    # below in L4
# # sort
from automol.chi.base._core import sorted_
from automol.chi.base._core import argsort
# L4
# # conversions
from automol.amchi._conv import amchi_key as inchi_key
from automol.amchi import smiles
from automol.amchi import graph
from automol.amchi import geometry
from automol.amchi import conformers
from automol.amchi import zmatrix
from automol.amchi import rdkit_molecule
from automol.amchi import rdkit_reaction
from automol.amchi import display
from automol.amchi import display_reaction
# # derived properties
from automol.chi._conv import is_complete
from automol.chi._conv import is_bad
from automol.chi._conv import chi_ as chi
from automol.chi._conv import without_stereo
# # derived transformations
from automol.chi._conv import join
from automol.chi._conv import add_stereo
from automol.chi._conv import expand_stereo
from automol.chi._conv import is_canonical_enantiomer
from automol.chi._conv import is_canonical_enantiomer_reaction
from automol.chi._conv import canonical_enantiomer
from automol.chi._conv import canonical_enantiomer_reaction
from automol.chi._conv import inchi_to_amchi
# drawing tools
from automol.amchi import draw
from automol.amchi import draw_grid
# assessment tools
from automol.amchi import is_valid_multiplicity as is_valid_inchi_multiplicity


def formula_sublayer(*args, **kwargs):
    """ Deprecated
    """
    return formula_layer(*args, **kwargs)


# keep this as an alias
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
    'formula_sublayer',
    'formula_layer',
    'main_sublayers',
    'charge_sublayers',
    'stereo_sublayers',
    'isotope_sublayers',
    # 'stereo_atoms',
    # 'stereo_bonds',
    # 'unassigned_stereo_bonds',
    'are_enantiomers',
    'are_diastereomers',
    'reflect',
    'racemic',
    # # setters
    'with_inchi_prefix',
    # # conversions
    # 'inchi_key',
    # 'smiles',
    'formula',
    'formula_string',
    'connectivity',
    'without_stereo',
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
    # # # isotope layers
    'bond_isotope_stereo_parities',
    'atom_isotope_stereo_parities',
    'is_inverted_isotope_enantiomer',
    # # # other properties
    'is_standard_form',
    'has_multiple_components',
    'is_enantiomer',
    'is_enantiomer_list',
    'is_enantiomer_reaction',
    'is_chiral',
    'has_stereo',
    'has_mobile_hydrogens',
    'low_spin_multiplicity',
    # # comparisons
    'same_connectivity',
    'equivalent',
    # # split/join
    'split',
    # 'join',
    # # sort
    'sorted_',
    'argsort',
    # L4
    # # conversions
    'inchi_key',
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
    'is_bad',
    'chi',
    'without_stereo',
    # # derived transformations
    'join',
    'add_stereo',
    'expand_stereo',
    'is_canonical_enantiomer',
    'is_canonical_enantiomer_reaction',
    'is_enantiomer_reaction',
    'canonical_enantiomer',
    'canonical_enantiomer_reaction',
    'inchi_to_amchi',
    # drawing tools
    'draw',
    'draw_grid',
    # assessment tools
    'is_valid_inchi_multiplicity',
]
