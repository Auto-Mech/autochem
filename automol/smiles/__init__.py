""" SMILES (Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""
# L2
# # conversions
from .base._core import without_resonance_stereo
from .base._core import without_stereo
from .base._core import reflect
# # split/join
from .base._core import split
from .base._core import join
from .base._core import reaction
from .base._core import is_reaction
from .base._core import reaction_reagents
from .base._core import reaction_reactant
from .base._core import reaction_product
from .base._core import reaction_reactants
from .base._core import reaction_products
from .base._core import reaction_reactants_and_products
# # properties
from .base._core import parse_connected_molecule_properties
# L4
# # conversions
from ._conv import amchi
from ._conv import inchi
from ._conv import chi
from ._conv import graph
from ._conv import geometry
from ._conv import formula_string
from ._conv import recalculate_without_stereo
from ._conv import rdkit_molecule
from ._conv import svg_string
from ._conv import reagent_svg_strings
from ._conv import reaction_reagent_svg_strings
from ._conv import rdkit_reaction
from ._conv import display
from ._conv import display_reaction


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
    'is_reaction',
    'reaction_reagents',
    'reaction_reactant',
    'reaction_product',
    'reaction_reactants',
    'reaction_products',
    'reaction_reactants_and_products',
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
    'reagent_svg_strings',
    'reaction_reagent_svg_strings',
    'rdkit_reaction',
    'display',
    'display_reaction',
]
