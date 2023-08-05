""" RSMILES (Resonance Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""
# L2
# # conversions
from automol.smiles.base._core import without_resonance_stereo
from automol.smiles.base._core import without_stereo
from automol.smiles.base._core import reflect
# # split/join
from automol.smiles.base._core import split
from automol.smiles.base._core import join
from automol.smiles.base._core import reaction
# # properties
from automol.smiles.base._core import parse_connected_molecule_properties
# L4
# # conversions
from automol.smiles._conv import amchi
from automol.smiles._conv import inchi
from automol.smiles._conv import chi
from automol.smiles._conv import graph
from automol.smiles._conv import geometry
from automol.smiles._conv import formula_string
from automol.smiles._conv import recalculate_without_stereo
from automol.smiles._conv import rdkit_molecule
from automol.smiles._conv import svg_string
from automol.smiles._conv import rdkit_reaction
from automol.smiles._conv import display
from automol.smiles._conv import display_reaction


__all__ = [
    # L2
    # # conversions
    'without_resonance_stereo',
    'without_stereo',
    'reflect',
    # # split/join
    'split',
    'join',
    'reaction',
    # # properties
    'parse_connected_molecule_properties',
    # L4
    # # conversions
    'amchi',
    'inchi',
    'chi',
    'graph',
    'geometry',
    'formula_string',
    'recalculate_without_stereo',
    'rdkit_molecule',
    'svg_string',
    'rdkit_reaction',
    'display',
    'display_reaction',
]
