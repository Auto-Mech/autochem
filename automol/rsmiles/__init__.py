""" RSMILES (Resonance Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""
# L2
# # properties
from automol.rsmiles.base._core import parse_properties
# L4
# # conversions
from automol.rsmiles._conv import connected_graph


__all__ = [
    # L2
    # # properties
    'parse_properties',
    # L4
    # # conversions
    'connected_graph',
]
