""" Level 3 ChI functions (depends on chi.base and extern)

A wrapper to redirect to either InChI or AMChI, as needed.

Excludes some functions that should no longer be used.

Currently only redirecting for resonance double bond stereo.
"""

# # "constructor"
from ._core import from_data
# # recalculate/standardize
from ._core import recalculate
from ._core import standard_form
# # getters
from ...amchi.base import version
# from ._core import formula_sublayer  # defined below
from ...amchi.base import formula_layer
from ...amchi.base import main_layers as main_sublayers
from ...amchi.base import charge_layers as charge_sublayers
from ...amchi.base import stereo_layers as stereo_sublayers
from ...amchi.base import isotope_layers as isotope_sublayers
# from ._core import stereo_atoms               # excluded
# from ._core import stereo_bonds               # excluded
# from ._core import unassigned_stereo_bonds    # excluded
from ._core import are_enantiomers
from ._core import are_diastereomers
from ._core import reflect
from ._core import racemic
# # setters
from ...amchi.base._core import with_inchi_prefix
# # conversions
# from ._core import inchi_key  # goes in L4
# from ._core import smiles     # goes in L4
from ...amchi.base import formula
from ...amchi.base import formula_string
# # # properties
from ._core import is_standard_form
from ...amchi.base import has_multiple_components
from ...amchi.base import is_enantiomer
from ...amchi.base import has_stereo
from ...amchi.base import has_mobile_hydrogens
# from ...amchi.base import low_spin_multiplicity          # excluded
# # comparisons
from ...amchi.base import same_connectivity
from ...amchi.base import equivalent
# # split/join
from ._core import split
# # sort
from ._core import sorted_
from ._core import argsort


def formula_sublayer(*args, **kwargs):
    """ Deprecated
    """
    return formula_layer(*args, **kwargs)


__all__ = [
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
    # # # properties
    'is_standard_form',
    'has_multiple_components',
    'is_enantiomer',
    'has_stereo',
    'has_mobile_hydrogens',
    # 'low_spin_multiplicity',
    # # comparisons
    'same_connectivity',
    'equivalent',
    # # split/join
    'split',
    # # sort
    'sorted_',
    'argsort',
]
