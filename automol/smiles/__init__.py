""" RSMILES (Resonance Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""
# L2
# # split/join
from automol.smiles.base._core import split
from automol.smiles.base._core import join
# # properties
from automol.smiles.base._core import parse_connected_molecule_properties
# L4
# # conversions
from automol.smiles._conv import amchi
from automol.smiles._conv import inchi
from automol.smiles._conv import chi
from automol.smiles._conv import graph


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
    'chi',
    'graph',
]
