""" SMILES (Simplified Molecular Input Line Entry System) strings

SMILES, with an extension for resonance double-bond stereo.
"""

# L2
# # conversions

# L4
# # conversions
from ._conv import (
    amchi,
    chi,
    display,
    display_reaction,
    formula_string,
    geometry,
    graph,
    inchi,
    rdkit_molecule,
    rdkit_reaction,
    reaction_reagent_svg_strings,
    reagent_svg_strings,
    recalculate_without_stereo,
    svg_string,
)

# # split/join
# # properties
from .base._core import (
    is_reaction,
    join,
    parse_connected_molecule_properties,
    reaction,
    reaction_product,
    reaction_products,
    reaction_reactant,
    reaction_reactants,
    reaction_reactants_and_products,
    reaction_reagents,
    reflect,
    split,
    without_resonance_stereo,
    without_stereo,
)

__all__ = [
    # L2
    # # conversions
    "without_resonance_stereo",
    "without_stereo",
    "reflect",
    # # split/join
    "split",
    "join",
    "reaction",
    "is_reaction",
    "reaction_reagents",
    "reaction_reactant",
    "reaction_product",
    "reaction_reactants",
    "reaction_products",
    "reaction_reactants_and_products",
    # # properties
    "parse_connected_molecule_properties",
    # L4
    # # conversions
    "amchi",
    "inchi",
    "chi",
    "graph",
    "geometry",
    "formula_string",
    "recalculate_without_stereo",
    "rdkit_molecule",
    "svg_string",
    "reagent_svg_strings",
    "reaction_reagent_svg_strings",
    "rdkit_reaction",
    "display",
    "display_reaction",
]
