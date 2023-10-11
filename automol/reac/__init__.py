""" reaction-class specific functionality
"""
# base reaction class
from automol.reac._0core import Reaction
# # constructors
from automol.reac._0core import from_data
from automol.reac._0core import from_string
from automol.reac._0core import string
# # getters
from automol.reac._0core import ts_graph
from automol.reac._0core import reactants_keys
from automol.reac._0core import products_keys
from automol.reac._0core import class_
from automol.reac._0core import ts_structure
from automol.reac._0core import reactant_structures
from automol.reac._0core import product_structures
from automol.reac._0core import structure_type
# # setters
from automol.reac._0core import set_ts_graph
from automol.reac._0core import set_reactants_keys
from automol.reac._0core import set_products_keys
from automol.reac._0core import set_reaction_class
from automol.reac._0core import set_ts_structure
from automol.reac._0core import set_reactant_structures
from automol.reac._0core import set_product_structures
from automol.reac._0core import set_structure_type
from automol.reac._0core import set_structures
# # others
from automol.reac._0core import reverse
from automol.reac._0core import mapping
from automol.reac._0core import reaction_mapping
from automol.reac._0core import reagent_sort_order
from automol.reac._0core import reactant_graphs
from automol.reac._0core import product_graphs
from automol.reac._0core import reactants_graph
from automol.reac._0core import products_graph
from automol.reac._0core import standard_keys
from automol.reac._0core import standard_keys_with_sorted_geometries
from automol.reac._0core import relabel
from automol.reac._0core import without_stereo
from automol.reac._0core import apply_dummy_conversion
from automol.reac._0core import reverse_dummy_conversion
from automol.reac._0core import insert_dummy_atoms
from automol.reac._0core import without_dummy_atoms
from automol.reac._0core import relabel_for_geometry
from automol.reac._0core import is_radical_radical
from automol.reac._0core import unique
# stereo-specific reactions
from automol.reac._2stereo import expand_stereo
from automol.reac._2stereo import expand_stereo_for_reaction
from automol.reac._2stereo import from_old_string
from automol.reac._2stereo import reflect
# finders
from automol.reac._3find import trivial
from automol.reac._3find import hydrogen_migrations
from automol.reac._3find import beta_scissions
from automol.reac._3find import ring_forming_scissions
from automol.reac._3find import eliminations
from automol.reac._3find import hydrogen_abstractions
from automol.reac._3find import additions
from automol.reac._3find import insertions
from automol.reac._3find import substitutions
from automol.reac._3find import find
from automol.reac._3find import find_from_chi
from automol.reac._3find import intersystem_crossing
# TS geometries
from automol.reac._4geom import with_geom_structures
from automol.reac._4geom import ts_geometry_from_reactants
# TS zmatrices
from automol.reac._5zmat import ts_zmatrix
from automol.reac._5zmat import zmatrix_coordinate_names
# rotational bonds & torsions
from automol.reac._7rot import linear_atom_keys
from automol.reac._7rot import rotational_bond_keys
from automol.reac._7rot import rotational_groups
from automol.reac._7rot import rotational_symmetry_number
# conversions
from automol.reac._8conv import amchi
from automol.reac._8conv import ts_amchi
from automol.reac._8conv import inchi
from automol.reac._8conv import chi
from automol.reac._8conv import smiles
from automol.reac._8conv import reaction_smiles
from automol.reac._8conv import rdkit_reaction
from automol.reac._8conv import display
from automol.reac._8conv import is_canonical_enantiomer
from automol.reac._8conv import canonical_enantiomer
from automol.reac._8conv import with_structures_from_chi
from automol.reac._8conv import with_structures_from_smiles
from automol.reac._8conv import with_structures_from_geometry
# scan coordinates
from automol.reac._6scan import build_scan_info
from automol.reac._6scan import scan_coordinate
from automol.reac._6scan import constraint_coordinates
# reaction products
from automol.reac._enum import enumerate_reactions
# species instability transformations
from automol.reac._instab import instability_product_zmas
from automol.reac._instab import instability_product_inchis
from automol.reac._instab import instability_product_graphs
from automol.reac._instab import instability_transformation
# phase space theory
from automol.reac._pst import pst_kt
from automol.reac._pst import pst_cn
# tunneling treatments
from automol.reac import tunnel
# comp functions
from automol.reac._comp import similar_saddle_point_structure


__all__ = [
    # base reaction class
    'Reaction',
    # # constructors
    'from_data',
    'from_string',
    'string',
    # # getters
    'ts_graph',
    'reactants_keys',
    'products_keys',
    'class_',
    'ts_structure',
    'reactant_structures',
    'product_structures',
    'structure_type',
    # # setters
    'set_ts_graph',
    'set_reactants_keys',
    'set_products_keys',
    'set_reaction_class',
    'set_ts_structure',
    'set_reactant_structures',
    'set_product_structures',
    'set_structure_type',
    'set_structures',
    # # others
    'reverse',
    'mapping',
    'reaction_mapping',
    'reagent_sort_order',
    'reactant_graphs',
    'product_graphs',
    'reactants_graph',
    'products_graph',
    'standard_keys',
    'standard_keys_with_sorted_geometries',
    'relabel',
    'without_stereo',
    'apply_dummy_conversion',
    'reverse_dummy_conversion',
    'insert_dummy_atoms',
    'without_dummy_atoms',
    'relabel_for_geometry',
    'is_radical_radical',
    'unique',
    # stereo-specific reactions
    'expand_stereo',
    'expand_stereo_for_reaction',
    'from_old_string',
    'reflect',
    # finders
    'trivial',
    'hydrogen_migrations',
    'beta_scissions',
    'ring_forming_scissions',
    'eliminations',
    'hydrogen_abstractions',
    'additions',
    'insertions',
    'substitutions',
    'find',
    'find_from_chi',
    'intersystem_crossing',
    # TS geometries
    'with_geom_structures',
    'ts_geometry_from_reactants',
    # TS zmatrices
    'ts_zmatrix',
    'zmatrix_coordinate_names',
    # rotational bonds & torsions
    'linear_atom_keys',
    'rotational_bond_keys',
    'rotational_groups',
    'rotational_symmetry_number',
    # conversions
    'amchi',
    'ts_amchi',
    'inchi',
    'chi',
    'smiles',
    'reaction_smiles',
    'rdkit_reaction',
    'display',
    'is_canonical_enantiomer',
    'canonical_enantiomer',
    'with_structures_from_chi',
    'with_structures_from_smiles',
    'with_structures_from_geometry',
    # scan coordinates
    'build_scan_info',
    'scan_coordinate',
    'constraint_coordinates',
    # reaction products
    'enumerate_reactions',
    # species instability transformations
    'instability_product_zmas',
    'instability_product_inchis',
    'instability_product_graphs',
    'instability_transformation',
    # phase space theory
    'pst_kt',
    'pst_cn',
    # tunneling treatments
    'tunnel',
    # comp functions
    'similar_saddle_point_structure',
]
