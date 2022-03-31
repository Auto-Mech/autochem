""" RSMILES (Resonance Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""
# L2
# # split/join
from automol.rsmiles.base._core import split
from automol.rsmiles.base._core import join
# # properties
from automol.rsmiles.base._core import parse_connected_molecule_properties
# L4
# # conversions
from automol.rsmiles._conv import amchi
from automol.rsmiles._conv import inchi
from automol.rsmiles._conv import graph


__all__ = [
    # L2
    # # split/join
    'split',
    'join',
    # # properties
    'parse_connected_molecule_properties',
    # L4
    # # conversions
    'amchi',
    'inchi',
    'graph',
]
