""" ChI strings

A wrapper to redirect to either InChI or AMChI, as needed.

Currently only redirecting for resonance double bond stereo.
"""

# L3
# # "constructor"
from .base._core import from_data
# # recalculate/standardize
from .base._core import recalculate
from .base._core import standard_form
# # getters
from ..amchi.base import version
# from .base._core import formula_sublayer  # defined below
from ..amchi.base import formula_layer
from ..amchi.base import main_layers as main_sublayers
from ..amchi.base import charge_layers as charge_sublayers
from ..amchi.base import stereo_layers as stereo_sublayers
from ..amchi.base import isotope_layers as isotope_sublayers
# from .base._core import stereo_atoms               # excluded
# from .base._core import stereo_bonds               # excluded
# from .base._core import unassigned_stereo_bonds    # excluded
from .base._core import are_enantiomers
from .base._core import are_diastereomers
from .base._core import reflect
from .base._core import racemic
# # setters
from ..amchi.base._core import with_inchi_prefix
# # conversions
# from .base._core import inchi_key  # below in L4
# from .base._core import smiles     # below in L4
from ..amchi.base import formula
from ..amchi.base import formula_string
from ..amchi.base import connectivity
# # properties
# # # formula layer
from ..amchi.base._core import symbols
from ..amchi.base._core import canonical_indices
# # # main layers
from ..amchi.base._core import bonds
from ..amchi.base._core import hydrogen_valences
# # # charge layers
from ..amchi.base._core import charge
# # # stereo layers
from ..amchi.base._core import bond_stereo_parities
from ..amchi.base._core import atom_stereo_parities
from ..amchi.base._core import is_inverted_enantiomer
# # # isotope layers
from ..amchi.base._core import bond_isotope_stereo_parities
from ..amchi.base._core import atom_isotope_stereo_parities
from ..amchi.base._core import is_inverted_isotope_enantiomer
# # # other properties
from .base._core import is_standard_form
from ..amchi.base import has_multiple_components
from ..amchi.base import is_enantiomer
from ..amchi.base import is_enantiomer_list
from ..amchi.base import is_enantiomer_reaction
from ..amchi.base import has_stereo
from ..amchi.base import has_mobile_hydrogens
from ..amchi.base import low_spin_multiplicity
# # comparisons
from ..amchi.base import same_connectivity
from ..amchi.base import equivalent
# # split/join
from .base._core import split
# from .base._core import split    # below in L4
# # sort
from .base._core import sorted_
from .base._core import argsort
# L4
# # conversions
from ..amchi._conv import amchi_key as inchi_key
from ..amchi import smiles
from ..amchi import graph
from ..amchi import geometry
from ..amchi import zmatrix
from ..amchi import rdkit_molecule
from ..amchi import rdkit_reaction
from ..amchi import display
from ..amchi import display_reaction
# # derived properties
from ._conv import is_complete
from ._conv import is_bad
from ._conv import chi_ as chi
from ._conv import without_stereo
# # derived transformations
from ._conv import join
from ._conv import add_stereo
from ._conv import expand_stereo
from ._conv import is_canonical_enantiomer
from ._conv import is_canonical_enantiomer_reaction
from ._conv import canonical_enantiomer
from ._conv import canonical_enantiomer_reaction
from ._conv import inchi_to_amchi
# drawing tools
from ..amchi import draw
from ..amchi import draw_grid
# assessment tools
from ..amchi import is_valid_multiplicity as is_valid_inchi_multiplicity


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
