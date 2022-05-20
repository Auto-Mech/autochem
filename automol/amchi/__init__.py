""" AMChI (AutoMech Chemical Identifier) strings
"""

# L2
# # constructor
from automol.amchi.base._core import from_data
# # getters
from automol.amchi.base._core import prefix
from automol.amchi.base._core import version
from automol.amchi.base._core import formula_string
from automol.amchi.base._core import main_layers
from automol.amchi.base._core import charge_layers
from automol.amchi.base._core import stereo_layers
from automol.amchi.base._core import isotope_layers
# # conversions
from automol.amchi.base._core import formula
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
# # other properties
from automol.amchi.base._core import has_stereo
from automol.amchi.base._core import has_multiple_components
from automol.amchi.base._core import low_spin_multiplicity
# # split/join
from automol.amchi.base._core import split
from automol.amchi.base._core import join
# L5
# # conversions
from automol.amchi._conv import connected_graph


__all__ = [
    # L2
    # # constructor
    'from_data',
    # # getters
    'prefix',
    'version',
    'formula_string',
    'main_layers',
    'charge_layers',
    'stereo_layers',
    'isotope_layers',
    # # conversions
    'formula',
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
    # # other properties
    'has_stereo',
    'has_multiple_components',
    'low_spin_multiplicity',
    # # split/join
    'split',
    'join',
    # L5
    # # conversions
    'connected_graph',
]
