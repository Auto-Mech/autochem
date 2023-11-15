""" reaction-class specific functionality
"""
# base reaction class
from automol.reac._0core import Reaction
# # constructors
from automol.reac._0core import from_data
from automol.reac._0core import from_forward_reverse
from automol.reac._0core import from_string
from automol.reac._2stereo import from_old_string
from automol.reac._2stereo import from_string_transitional
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
from automol.reac._0core import set_structures
from automol.reac._0core import update_structures
# # others
from automol.reac._0core import reverse_without_recalculating
from automol.reac._0core import mapping
from automol.reac._0core import reactant_mappings
from automol.reac._0core import product_mappings
from automol.reac._0core import reactant_graphs
from automol.reac._0core import product_graphs
from automol.reac._0core import reactants_graph
from automol.reac._0core import products_graph
from automol.reac._0core import relabel
from automol.reac._0core import without_stereo
from automol.reac._0core import without_structures
from automol.reac._0core import apply_zmatrix_conversion
from automol.reac._0core import undo_zmatrix_conversion
from automol.reac._0core import without_dummy_atoms
from automol.reac._0core import is_radical_radical
from automol.reac._0core import unique
# stereo-specific reactions
from automol.reac._2stereo import expand_stereo
from automol.reac._2stereo import expand_stereo_to_match_reagents
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
# TS geometries
from automol.reac._4struc import with_structures
from automol.reac._4struc import reverse
# TS zmatrices
from automol.reac._deprecated import zmatrix_coordinate_names
# conversions
# # constructors from data types
from automol.reac._5conv import from_graphs
from automol.reac._5conv import from_amchis
from automol.reac._5conv import from_inchis
from automol.reac._5conv import from_chis
from automol.reac._5conv import from_smiles
from automol.reac._5conv import from_geometries
from automol.reac._5conv import from_zmatrices
# # converters to various data types
from automol.reac._5conv import graphs
from automol.reac._5conv import amchis
from automol.reac._5conv import inchis
from automol.reac._5conv import chis
from automol.reac._5conv import smiles
from automol.reac._5conv import geometries
from automol.reac._5conv import zmatrices
# # additional data types
from automol.reac._5conv import ts_amchi
from automol.reac._5conv import reaction_smiles
from automol.reac._5conv import rdkit_reaction
from automol.reac._5conv import display
# # canonicity
from automol.reac._5conv import is_canonical_enantiomer
from automol.reac._5conv import canonical_enantiomer
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
    'from_forward_reverse',
    'from_string',
    'from_old_string',
    'from_string_transitional',
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
    'set_structures',
    'update_structures',
    # # others
    'reverse_without_recalculating',
    'mapping',
    'reactant_mappings',
    'product_mappings',
    'reactant_graphs',
    'product_graphs',
    'reactants_graph',
    'products_graph',
    'relabel',
    'without_stereo',
    'without_structures',
    'apply_zmatrix_conversion',
    'undo_zmatrix_conversion',
    'without_dummy_atoms',
    'is_radical_radical',
    'unique',
    # stereo-specific reactions
    'expand_stereo',
    'expand_stereo_to_match_reagents',
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
    # TS geometries
    'with_structures',
    'reverse',
    # TS zmatrices
    'zmatrix_coordinate_names',
    # conversions
    # # constructors from data types
    'from_graphs',
    'from_amchis',
    'from_inchis',
    'from_chis',
    'from_smiles',
    'from_geometries',
    'from_zmatrices',
    # # converters to various data types
    'graphs',
    'amchis',
    'inchis',
    'chis',
    'smiles',
    'geometries',
    'zmatrices',
    # # additional data types
    'ts_amchi',
    'reaction_smiles',
    'rdkit_reaction',
    'display',
    # # canonicity
    'is_canonical_enantiomer',
    'canonical_enantiomer',
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
