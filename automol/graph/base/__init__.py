""" Level 2 graph functions (no dependencies on extern or other types)

Import hierarchy:
    _core       no dependencies
    _networkx   dependencies: _core
    _algo       dependencies: _core, _networkx
    _resonance  dependencies: _core, _networkx, _algo
    _rot        dependencies: _core, _networkx, _algo, _resonance
    _stereo     dependencies: _core, _networkx, _algo, _resonance
    _canon      dependencies: _core, _networkx, _algo, _resonance
    _func_group dependencies: _core, _networkx, _algo, _resonance, _stereo
    _amchi      dependencies: _core, _networkx, _algo, _canon, _resonance
    _smiles     dependencies: _core, _networkx, _algo, _canon, _resonance
    ts          dependencies: _core, _networkx, _algo, _resonance, _stereo

Each next submodule in the hierarchy may depend on the ones before it, but
**not** the ones after. This avoids circular dependencies.

Please obey this hierarchy strictly and raise an issue to the group if
something needs to be updated or altered in it.
"""

# core functions:
# # constructors
from automol.graph.base._core import from_data
from automol.graph.base._core import atoms_from_data
from automol.graph.base._core import bonds_from_data
from automol.graph.base._core import from_atoms_and_bonds
# # getters
from automol.graph.base._core import atoms
from automol.graph.base._core import bonds
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
# # setters
from automol.graph.base._core import set_atom_symbols
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import set_atom_implicit_hydrogen_valences
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
# # I/O
from automol.graph.base._core import string
from automol.graph.base._core import yaml_dictionary
from automol.graph.base._core import from_string
from automol.graph.base._core import from_yaml_dictionary
# # conversions
from automol.graph.base._core import frozen
from automol.graph.base._core import formula
# # properties
from automol.graph.base._core import atom_count
from automol.graph.base._core import atom_count_by_type
from automol.graph.base._core import heavy_atom_count
from automol.graph.base._core import electron_count
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import has_stereo
from automol.graph.base._core import atomic_numbers
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import atom_element_valences
from automol.graph.base._core import atom_lone_pair_counts
from automol.graph.base._core import atom_van_der_waals_radius
from automol.graph.base._core import atom_bond_valences
from automol.graph.base._core import atom_unsaturated_valences
from automol.graph.base._core import atom_explicit_hydrogen_valences
from automol.graph.base._core import atom_hybridizations
from automol.graph.base._core import tetrahedral_atom_keys
from automol.graph.base._core import maximum_spin_multiplicity
from automol.graph.base._core import possible_spin_multiplicities
from automol.graph.base._core import atom_symbol_keys
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import atom_explicit_hydrogen_keys
from automol.graph.base._core import explicit_hydrogen_keys
from automol.graph.base._core import terminal_atom_keys
from automol.graph.base._core import terminal_heavy_atom_keys
from automol.graph.base._core import unsaturated_atom_keys
from automol.graph.base._core import lone_pair_atom_keys
from automol.graph.base._core import angle_keys
# # relabeling and changing keys
from automol.graph.base._core import relabel
from automol.graph.base._core import transform_keys
from automol.graph.base._core import standard_keys
from automol.graph.base._core import standard_keys_for_sequence
from automol.graph.base._core import relabel_for_zmatrix
from automol.graph.base._core import relabel_for_geometry
# # add/remove/insert/without
from automol.graph.base._core import add_atoms
from automol.graph.base._core import add_bonds
from automol.graph.base._core import remove_atoms
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import remove_atom_stereo_parities
from automol.graph.base._core import remove_bond_stereo_parities
from automol.graph.base._core import add_atom_implicit_hydrogen_valences
from automol.graph.base._core import add_atom_explicit_hydrogen_keys
from automol.graph.base._core import add_bonded_atom
from automol.graph.base._core import insert_bonded_atom
from automol.graph.base._core import add_dummy_atoms
from automol.graph.base._core import insert_dummy_atoms
from automol.graph.base._core import standard_keys_without_dummy_atoms
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import without_fractional_bonds
from automol.graph.base._core import without_dummy_bonds
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import explicit
from automol.graph.base._core import implicit
# # unions
from automol.graph.base._core import union
from automol.graph.base._core import union_from_sequence
# # subgraphs and neighborhoods
from automol.graph.base._core import subgraph
from automol.graph.base._core import bond_induced_subgraph
from automol.graph.base._core import atom_neighborhood
from automol.graph.base._core import atom_neighborhoods
from automol.graph.base._core import bond_neighborhood
from automol.graph.base._core import bond_neighborhoods
from automol.graph.base._core import atom_neighbor_atom_key
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atom_sorted_neighbor_atom_keys
from automol.graph.base._core import atoms_sorted_neighbor_atom_keys
from automol.graph.base._core import atoms_second_degree_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import dummy_atoms_neighbor_atom_key
from automol.graph.base._core import bonds_neighbor_atom_keys
from automol.graph.base._core import bonds_neighbor_bond_keys
# algorithm functions:
# # isomorphisms and equivalence
from automol.graph.base._algo import isomorphism
from automol.graph.base._algo import sequence_isomorphism
from automol.graph.base._algo import full_isomorphism
from automol.graph.base._algo import full_subgraph_isomorphism
from automol.graph.base._algo import backbone_isomorphism
from automol.graph.base._algo import backbone_isomorphic
from automol.graph.base._algo import backbone_unique
from automol.graph.base._algo import equivalent_atoms
from automol.graph.base._algo import equivalent_bonds
from automol.graph.base._algo import are_equivalent_atoms
from automol.graph.base._algo import are_equivalent_bonds
from automol.graph.base._algo import atom_equivalence_class_reps
from automol.graph.base._algo import bond_equivalence_class_reps
from automol.graph.base._algo import chem_unique_atoms_of_type
# # algorithms
from automol.graph.base._algo import connected_components
from automol.graph.base._algo import connected_components_atom_keys
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import atom_shortest_paths
from automol.graph.base._algo import shortest_path_between_atoms
from automol.graph.base._algo import shortest_path_between_groups
from automol.graph.base._algo import atom_longest_chains
from automol.graph.base._algo import atom_longest_chain
from automol.graph.base._algo import longest_chain
# # branches and groups
from automol.graph.base._algo import ring_atom_chirality
from automol.graph.base._algo import atom_groups
from automol.graph.base._algo import branch
from automol.graph.base._algo import branch_atom_keys
from automol.graph.base._algo import branch_bond_keys
from automol.graph.base._algo import is_branched
# # rings
from automol.graph.base._algo import rings
from automol.graph.base._algo import rings_atom_keys
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._algo import sorted_ring_atom_keys
from automol.graph.base._algo import sorted_ring_atom_keys_from_bond_keys
from automol.graph.base._algo import is_ring_key_sequence
from automol.graph.base._algo import cycle_ring_atom_key_to_front
from automol.graph.base._algo import ring_arc_complement_atom_keys
from automol.graph.base._algo import ring_systems
from automol.graph.base._algo import ring_systems_atom_keys
from automol.graph.base._algo import ring_systems_bond_keys
from automol.graph.base._algo import is_ring_system
from automol.graph.base._algo import ring_system_decomposed_atom_keys
from automol.graph.base._algo import ring_systems_decomposed_atom_keys
# resonance functions:
# # core functions
from automol.graph.base._resonance import dominant_resonance
from automol.graph.base._resonance import dominant_resonances
from automol.graph.base._resonance import resonances
from automol.graph.base._resonance import subresonances
from automol.graph.base._resonance import resonance_dominant_bond_orders
from automol.graph.base._resonance import one_resonance_dominant_bond_orders
from automol.graph.base._resonance import resonance_avg_bond_orders
# # derived properties
from automol.graph.base._resonance import linear_atom_keys
from automol.graph.base._resonance import linear_segments_atom_keys
from automol.graph.base._resonance import radical_atom_keys
from automol.graph.base._resonance import radical_atom_keys_from_resonance
from automol.graph.base._resonance import has_separated_radical_sites
from automol.graph.base._resonance import nonresonant_radical_atom_keys
from automol.graph.base._resonance import vinyl_radical_atom_keys
from automol.graph.base._resonance import sigma_radical_atom_keys
from automol.graph.base._resonance import resonance_dominant_radical_atom_keys
from automol.graph.base._resonance import sing_res_dom_radical_atom_keys
from automol.graph.base._resonance import radical_groups
from automol.graph.base._resonance import radical_group_dct
from automol.graph.base._resonance import sp2_bond_keys
from automol.graph.base._resonance import (
        resonance_dominant_atom_hybridizations)
from automol.graph.base._resonance import (
        resonance_dominant_atom_centered_cumulene_keys)
from automol.graph.base._resonance import (
        resonance_dominant_bond_centered_cumulene_keys)
# canonicalization functions:
# # canonical key functions
from automol.graph.base._canon import canonical
from automol.graph.base._canon import canonical_keys
# # canonical stereo functions
from automol.graph.base._canon import reflect
from automol.graph.base._canon import reflect_local_stereo
from automol.graph.base._canon import to_local_stereo
from automol.graph.base._canon import from_local_stereo
from automol.graph.base._canon import from_geometry
# # symmetry class functions
from automol.graph.base._canon import class_indices
from automol.graph.base._canon import class_indices_and_stereo_parities
# # parity evaluators
from automol.graph.base._canon import atom_parity_evaluator_from_geometry_
from automol.graph.base._canon import bond_parity_evaluator_from_geometry_
from automol.graph.base._canon import atom_parity_evaluator_to_local_stereo_
from automol.graph.base._canon import bond_parity_evaluator_to_local_stereo_
from automol.graph.base._canon import atom_parity_evaluator_from_local_stereo_
from automol.graph.base._canon import bond_parity_evaluator_from_local_stereo_
# functional groups code:
# # core functions
from automol.graph.base._func_group import FunctionalGroup
from automol.graph.base._func_group import functional_group_dct
from automol.graph.base._func_group import functional_group_count_dct
# # finders for overaching types
from automol.graph.base._func_group import hydrocarbon_species
from automol.graph.base._func_group import radical_species
# # finders for reactive sites and groups
from automol.graph.base._func_group import alkene_sites
from automol.graph.base._func_group import alkyne_sites
from automol.graph.base._func_group import alcohol_groups
from automol.graph.base._func_group import peroxy_groups
from automol.graph.base._func_group import hydroperoxy_groups
from automol.graph.base._func_group import ether_groups
from automol.graph.base._func_group import epoxy_groups
from automol.graph.base._func_group import aldehyde_groups
from automol.graph.base._func_group import ketone_groups
from automol.graph.base._func_group import ester_groups
from automol.graph.base._func_group import carboxylic_acid_groups
from automol.graph.base._func_group import amide_groups
from automol.graph.base._func_group import nitro_groups
from automol.graph.base._func_group import halide_groups
from automol.graph.base._func_group import thiol_groups
from automol.graph.base._func_group import radical_dissociation_products
# # helper functions
from automol.graph.base._func_group import bonds_of_type
from automol.graph.base._func_group import bonds_of_order
from automol.graph.base._func_group import two_bond_idxs
from automol.graph.base._func_group import neighbors_of_type
from automol.graph.base._func_group import radicals_of_type
# torsion/rotational bond functions:
from automol.graph.base._rot import rotational_bond_keys
from automol.graph.base._rot import rotational_groups
from automol.graph.base._rot import rotational_symmetry_number
from automol.graph.base._rot import bond_symmetry_numbers
# stereo functions:
# # core functions
from automol.graph.base._stereo import stereo_priority_vector
from automol.graph.base._stereo import stereogenic_atom_keys
from automol.graph.base._stereo import stereogenic_bond_keys
from automol.graph.base._stereo import stereomers
from automol.graph.base._stereo import substereomers
from automol.graph.base._stereo import to_index_based_stereo
from automol.graph.base._stereo import from_index_based_stereo
# # derived properties
from automol.graph.base._stereo import atom_stereo_sorted_neighbor_atom_keys
from automol.graph.base._stereo import bond_stereo_sorted_neighbor_atom_keys
from automol.graph.base._stereo import atoms_stereo_sorted_neighbor_atom_keys
# # stereo setting code
from automol.graph.base._stereo import set_stereo_from_geometry
# # stereo parity evaluation code
from automol.graph.base._stereo import atom_stereo_parity_from_geometry
from automol.graph.base._stereo import bond_stereo_parity_from_geometry
# AMChI functions:
from automol.graph.base._amchi import amchi
# SMILES functions:
from automol.graph.base._smiles import rsmiles
# TS graph submodule:
from automol.graph.base import ts


__all__ = [
    # core functions:
    # # constructors
    'from_data',
    'atoms_from_data',
    'bonds_from_data',
    'from_atoms_and_bonds',
    # # getters
    'atoms',
    'bonds',
    'atom_keys',
    'bond_keys',
    'atom_symbols',
    'bond_orders',
    'atom_implicit_hydrogen_valences',
    'atom_stereo_parities',
    'bond_stereo_parities',
    # # setters
    'set_atom_symbols',
    'set_bond_orders',
    'set_atom_implicit_hydrogen_valences',
    'set_atom_stereo_parities',
    'set_bond_stereo_parities',
    # # I/O
    'string',
    'yaml_dictionary',
    'from_string',
    'from_yaml_dictionary',
    # # conversions
    'frozen',
    'formula',
    # # properties
    'atom_count',
    'atom_count_by_type',
    'heavy_atom_count',
    'electron_count',
    'atom_stereo_keys',
    'bond_stereo_keys',
    'has_stereo',
    'atomic_numbers',
    'mass_numbers',
    'atom_element_valences',
    'atom_lone_pair_counts',
    'atom_van_der_waals_radius',
    'atom_bond_valences',
    'atom_unsaturated_valences',
    'atom_explicit_hydrogen_valences',
    'atom_hybridizations',
    'tetrahedral_atom_keys',
    'maximum_spin_multiplicity',
    'possible_spin_multiplicities',
    'atom_symbol_keys',
    'backbone_keys',
    'atom_explicit_hydrogen_keys',
    'explicit_hydrogen_keys',
    'terminal_atom_keys',
    'terminal_heavy_atom_keys',
    'unsaturated_atom_keys',
    'lone_pair_atom_keys',
    'angle_keys',
    # # relabeling and changing keys
    'relabel',
    'transform_keys',
    'standard_keys',
    'standard_keys_for_sequence',
    'relabel_for_zmatrix',
    'relabel_for_geometry',
    # # add/remove/insert/without
    'add_atoms',
    'add_bonds',
    'remove_atoms',
    'remove_bonds',
    'remove_atom_stereo_parities',
    'remove_bond_stereo_parities',
    'add_atom_implicit_hydrogen_valences',
    'add_atom_explicit_hydrogen_keys',
    'add_bonded_atom',
    'insert_bonded_atom',
    'add_dummy_atoms',
    'insert_dummy_atoms',
    'standard_keys_without_dummy_atoms',
    'without_bond_orders',
    'without_dummy_atoms',
    'without_fractional_bonds',
    'without_dummy_bonds',
    'without_stereo_parities',
    'explicit',
    'implicit',
    # # unions
    'union',
    'union_from_sequence',
    # # subgraphs and neighborhoods
    'subgraph',
    'bond_induced_subgraph',
    'atom_neighborhood',
    'atom_neighborhoods',
    'bond_neighborhood',
    'bond_neighborhoods',
    'atom_neighbor_atom_key',
    'atoms_neighbor_atom_keys',
    'atom_sorted_neighbor_atom_keys',
    'atoms_sorted_neighbor_atom_keys',
    'atoms_second_degree_neighbor_atom_keys',
    'atoms_bond_keys',
    'dummy_atoms_neighbor_atom_key',
    'bonds_neighbor_atom_keys',
    'bonds_neighbor_bond_keys',
    # algorithm functions:
    # # isomorphisms and equivalence
    'isomorphism',
    'sequence_isomorphism',
    'full_isomorphism',
    'full_subgraph_isomorphism',
    'backbone_isomorphism',
    'backbone_isomorphic',
    'backbone_unique',
    'equivalent_atoms',
    'equivalent_bonds',
    'are_equivalent_atoms',
    'are_equivalent_bonds',
    'atom_equivalence_class_reps',
    'bond_equivalence_class_reps',
    'chem_unique_atoms_of_type',
    # # algorithms
    'connected_components',
    'connected_components_atom_keys',
    'is_connected',
    'atom_shortest_paths',
    'shortest_path_between_atoms',
    'shortest_path_between_groups',
    'atom_longest_chains',
    'atom_longest_chain',
    'longest_chain',
    # # branches and groups
    'ring_atom_chirality',
    'atom_groups',
    'branch',
    'branch_atom_keys',
    'branch_bond_keys',
    'is_branched',
    # # rings
    'rings',
    'rings_atom_keys',
    'rings_bond_keys',
    'sorted_ring_atom_keys',
    'sorted_ring_atom_keys_from_bond_keys',
    'is_ring_key_sequence',
    'cycle_ring_atom_key_to_front',
    'ring_arc_complement_atom_keys',
    'ring_systems',
    'ring_systems_atom_keys',
    'ring_systems_bond_keys',
    'is_ring_system',
    'ring_system_decomposed_atom_keys',
    'ring_systems_decomposed_atom_keys',
    # resonance functions:
    # # core functions
    'dominant_resonance',
    'dominant_resonances',
    'resonances',
    'subresonances',
    'resonance_dominant_bond_orders',
    'one_resonance_dominant_bond_orders',
    'resonance_avg_bond_orders',
    # # derived properties
    'linear_atom_keys',
    'linear_segments_atom_keys',
    'radical_atom_keys',
    'radical_atom_keys_from_resonance',
    'has_separated_radical_sites',
    'nonresonant_radical_atom_keys',
    'vinyl_radical_atom_keys',
    'sigma_radical_atom_keys',
    'resonance_dominant_radical_atom_keys',
    'sing_res_dom_radical_atom_keys',
    'radical_groups',
    'radical_group_dct',
    'sp2_bond_keys',
    'resonance_dominant_atom_hybridizations',
    'resonance_dominant_atom_centered_cumulene_keys',
    'resonance_dominant_bond_centered_cumulene_keys',
    # canonicalization functions:
    # # canonical key functions
    'canonical',
    'canonical_keys',
    # # canonical stereo functions
    'reflect',
    'reflect_local_stereo',
    'to_local_stereo',
    'from_local_stereo',
    'from_geometry',
    # # symmetry class functions
    'class_indices',
    'class_indices_and_stereo_parities',
    # # parity evaluators
    'atom_parity_evaluator_from_geometry_',
    'bond_parity_evaluator_from_geometry_',
    'atom_parity_evaluator_to_local_stereo_',
    'bond_parity_evaluator_to_local_stereo_',
    'atom_parity_evaluator_from_local_stereo_',
    'bond_parity_evaluator_from_local_stereo_',
    # functional groups code:
    # # core functions
    'FunctionalGroup',
    'functional_group_dct',
    'functional_group_count_dct',
    # # finders for overaching types
    'hydrocarbon_species',
    'radical_species',
    # # finders for reactive sites and groups
    'alkene_sites',
    'alkyne_sites',
    'alcohol_groups',
    'peroxy_groups',
    'hydroperoxy_groups',
    'ether_groups',
    'epoxy_groups',
    'aldehyde_groups',
    'ketone_groups',
    'ester_groups',
    'carboxylic_acid_groups',
    'amide_groups',
    'nitro_groups',
    'halide_groups',
    'thiol_groups',
    'radical_dissociation_products',
    # # helper functions
    'bonds_of_type',
    'bonds_of_order',
    'two_bond_idxs',
    'neighbors_of_type',
    'radicals_of_type',
    # torsion/rotational bond functions:
    'rotational_bond_keys',
    'rotational_groups',
    'rotational_symmetry_number',
    'bond_symmetry_numbers',
    # stereo functions:
    # # core functions
    'stereo_priority_vector',
    'stereogenic_atom_keys',
    'stereogenic_bond_keys',
    'stereomers',
    'substereomers',
    'to_index_based_stereo',
    'from_index_based_stereo',
    # # derived properties
    'atom_stereo_sorted_neighbor_atom_keys',
    'bond_stereo_sorted_neighbor_atom_keys',
    'atoms_stereo_sorted_neighbor_atom_keys',
    # # stereo setting code
    'set_stereo_from_geometry',
    # # stereo parity evaluation code
    'atom_stereo_parity_from_geometry',
    'bond_stereo_parity_from_geometry',
    # AMChI functions:
    'amchi',
    # SMILES functions:
    'rsmiles',
    # TS graph submodule:
    'ts',
]
