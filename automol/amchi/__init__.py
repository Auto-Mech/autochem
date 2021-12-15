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
from automol.amchi.base._core import is_chiral
from automol.amchi.base._core import is_enantiomer
from automol.amchi.base._core import has_stereo
from automol.amchi.base._core import has_multiple_components
from automol.amchi.base._core import charge
from automol.amchi.base._core import low_spin_multiplicity
from automol.amchi.base._core import canonical_indices
from automol.amchi.base._core import symbols
from automol.amchi.base._core import bonds
from automol.amchi.base._core import hydrogen_valences
from automol.amchi.base._core import stereo_atoms
from automol.amchi.base._core import stereo_bonds
from automol.amchi.base._core import unassigned_stereo_bonds
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
    'is_chiral',
    'is_enantiomer',
    'has_stereo',
    'has_multiple_components',
    'charge',
    'low_spin_multiplicity',
    'canonical_indices',
    'symbols',
    'bonds',
    'hydrogen_valences',
    'stereo_atoms',
    'stereo_bonds',
    'unassigned_stereo_bonds',
    # # split/join
    'split',
    'join',
    # L5
    # # conversions
    'connected_graph',
]
