""" Level 2 RSMILES functions
"""

# # split/join
from automol.smiles.base._core import split
from automol.smiles.base._core import join
# # properties
from automol.smiles.base._core import parse_connected_molecule_properties

__all__ = [
    # # split/join
    'split',
    'join',
    # # properties
    'parse_connected_molecule_properties',
]
