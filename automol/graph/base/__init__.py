""" Level 2 graph functions (no dependencies on extern or other types)

Import hierarchy:
    _core       no dependencies
    _networkx   dependencies: _core
    _algo       dependencies: _core, _networkx
    _kekule     dependencies: _core, _networkx, _algo
    _rot        dependencies: _core, _networkx, _algo, _kekule
    _geom       dependencies: _core, _networkx, _algo, _kekule
    _canon      dependencies: _core, _networkx, _algo, _kekule
    _amchi      dependencies: _core, _networkx, _algo, _canon, _kekule
    _smiles     dependencies: _core, _networkx, _algo, _canon, _kekule
    _stereo     dependencies: _core, _networkx, _algo, _canon, _kekule, _amchi
    _func_group dependencies: _core, _networkx, _algo, _kekule, _stereo
    ts          dependencies: _core, _networkx, _algo, _kekule, _stereo

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
from automol.graph.base._core import atom_implicit_hydrogens
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import ts_graph
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import stereo_parities
# # setters
from automol.graph.base._core import set_atom_symbols
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import set_atom_implicit_hydrogens
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import set_stereo_parities
# # I/O
from automol.graph.base._core import string
from automol.graph.base._core import yaml_dictionary
from automol.graph.base._core import from_string
from automol.graph.base._core import from_yaml_dictionary
from automol.graph.base._core import from_old_yaml_dictionary
# # conversions
from automol.graph.base._core import frozen
from automol.graph.base._core import formula
# # properties
from automol.graph.base._core import atom_count
from automol.graph.base._core import electron_count
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import stereo_keys
from automol.graph.base._core import has_stereo
from automol.graph.base._core import has_atom_stereo
from automol.graph.base._core import has_bond_stereo
from automol.graph.base._core import has_pi_bonds
from automol.graph.base._core import is_ts_graph
from automol.graph.base._core import atomic_numbers
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import atomic_valences
from automol.graph.base._core import atom_lone_pairs
from automol.graph.base._core import atom_van_der_waals_radius
from automol.graph.base._core import atom_bond_counts
from automol.graph.base._core import atom_unpaired_electrons
from automol.graph.base._core import bond_unpaired_electrons
from automol.graph.base._core import tetrahedral_atom_keys
from automol.graph.base._core import maximum_spin_multiplicity
from automol.graph.base._core import possible_spin_multiplicities
from automol.graph.base._core import atom_symbol_keys
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import backbone_bond_keys
from automol.graph.base._core import backbone_hydrogen_keys
from automol.graph.base._core import nonbackbone_hydrogen_keys
from automol.graph.base._core import atom_backbone_hydrogen_keys
from automol.graph.base._core import atom_nonbackbone_hydrogen_keys
from automol.graph.base._core import terminal_atom_keys
from automol.graph.base._core import unsaturated_atom_keys
from automol.graph.base._core import unsaturated_bond_keys
from automol.graph.base._core import lone_pair_atom_keys
from automol.graph.base._core import angle_keys
# # relabeling and changing keys
from automol.graph.base._core import relabel
from automol.graph.base._core import standard_keys
from automol.graph.base._core import standard_keys_for_sequence
from automol.graph.base._core import relabel_for_zmatrix
from automol.graph.base._core import relabel_for_geometry
# # add/remove/insert/without
from automol.graph.base._core import add_atoms
from automol.graph.base._core import add_bonds
from automol.graph.base._core import remove_atoms
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import change_implicit_hydrogens
from automol.graph.base._core import add_atom_explicit_hydrogens
from automol.graph.base._core import add_bonded_atom
from automol.graph.base._core import shift_insert_bonded_atom
from automol.graph.base._core import add_dummy_atoms
from automol.graph.base._core import shift_insert_dummy_atoms
from automol.graph.base._core import shift_remove_dummy_atoms
from automol.graph.base._core import without_pi_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import ts_reverse
from automol.graph.base._core import ts_reagents_graph_without_stereo
from automol.graph.base._core import without_bonds_by_orders
from automol.graph.base._core import without_stereo
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
from automol.graph.base._core import atom_neighbor_atom_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atom_sorted_neighbor_atom_keys
from automol.graph.base._core import atoms_sorted_neighbor_atom_keys
from automol.graph.base._core import atom_bond_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import dummy_atoms_neighbor_atom_key
from automol.graph.base._core import bonds_neighbor_atom_keys
from automol.graph.base._core import bonds_neighbor_bond_keys
# algorithm functions:
# # isomorphisms and equivalence
from automol.graph.base._algo import isomorphism
from automol.graph.base._algo import isomorphic
from automol.graph.base._algo import unique
from automol.graph.base._algo import sequence_isomorphism
from automol.graph.base._algo import subgraph_isomorphism
from automol.graph.base._algo import equivalent_atoms
from automol.graph.base._algo import equivalent_bonds
from automol.graph.base._algo import are_equivalent_atoms
from automol.graph.base._algo import are_equivalent_bonds
from automol.graph.base._algo import atom_equivalence_class_reps
from automol.graph.base._algo import bond_equivalence_class_reps
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
from automol.graph.base._algo import branches
from automol.graph.base._algo import branch_dict
from automol.graph.base._algo import branch_atom_keys
from automol.graph.base._algo import branch
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
# kekule functions:
# # core functions
from automol.graph.base._kekule import kekule
from automol.graph.base._kekule import kekules
from automol.graph.base._kekule import kekule_bond_orders
from automol.graph.base._kekule import kekules_bond_orders
from automol.graph.base._kekule import kekules_bond_orders_collated
from automol.graph.base._kekule import kekules_bond_orders_averaged
# # derived properties
from automol.graph.base._kekule import linear_atom_keys
from automol.graph.base._kekule import linear_segments_atom_keys
from automol.graph.base._kekule import atom_hybridizations
from automol.graph.base._kekule import atom_hybridizations_from_kekule
from automol.graph.base._kekule import radical_atom_keys
from automol.graph.base._kekule import radical_atom_keys_from_kekule
from automol.graph.base._kekule import nonresonant_radical_atom_keys
from automol.graph.base._kekule import vinyl_radical_atom_keys
from automol.graph.base._kekule import sigma_radical_atom_keys
from automol.graph.base._kekule import has_separated_radical_sites
from automol.graph.base._kekule import resonance_bond_stereo_keys
from automol.graph.base._kekule import vinyl_bond_stereo_keys
from automol.graph.base._kekule import has_resonance_bond_stereo
from automol.graph.base._kekule import has_vinyl_bond_stereo
from automol.graph.base._kekule import has_nonkekule_bond_stereo
from automol.graph.base._kekule import has_noninchi_stereo
from automol.graph.base._kekule import radical_groups
from automol.graph.base._kekule import radical_group_dct
from automol.graph.base._kekule import rigid_planar_bond_keys
from automol.graph.base._kekule import atom_centered_cumulene_keys
from automol.graph.base._kekule import bond_centered_cumulene_keys
# torsion/rotational bond functions:
from automol.graph.base._rot import rotational_bond_keys
from automol.graph.base._rot import rotational_groups
from automol.graph.base._rot import rotational_symmetry_number
# geometry functions:
# # stereo parity evaluations
from automol.graph.base._geom import geometry_atom_parity
from automol.graph.base._geom import geometry_bond_parity
from automol.graph.base._geom import geometry_local_parity
from automol.graph.base._geom import geometries_have_matching_parities
from automol.graph.base._geom import geometries_parity_mismatches
# # corrections
from automol.graph.base._geom import linear_vinyl_corrected_geometry
from automol.graph.base._geom import geometry_rotate_bond
from automol.graph.base._geom import geometry_dihedrals_near_value
from automol.graph.base._geom import perturb_geometry_planar_dihedrals
# canonicalization functions:
# # canonical key functions
from automol.graph.base._canon import canonical_enantiomer
from automol.graph.base._canon import canonical_enantiomer_with_keys
from automol.graph.base._canon import canonical
from automol.graph.base._canon import canonical_keys
# # canonical stereo functions
from automol.graph.base._canon import stereogenic_atom_keys
from automol.graph.base._canon import stereogenic_bond_keys
from automol.graph.base._canon import stereogenic_keys
from automol.graph.base._canon import stereogenic_atom_keys_from_priorities
from automol.graph.base._canon import stereogenic_bond_keys_from_priorities
from automol.graph.base._canon import reflect
from automol.graph.base._canon import reflect_local_stereo
from automol.graph.base._canon import to_local_stereo
from automol.graph.base._canon import from_local_stereo
from automol.graph.base._canon import set_stereo_from_geometry
# # symmetry class functions
from automol.graph.base._canon import canonical_priorities
from automol.graph.base._canon import calculate_priorities_and_assign_stereo
# # parity evaluators
from automol.graph.base._canon import parity_evaluator_from_geometry_
from automol.graph.base._canon import parity_evaluator_read_canonical_
from automol.graph.base._canon import parity_evaluator_flip_local_
# AMChI functions:
from automol.graph.base._amchi import amchi
from automol.graph.base._amchi import amchi_with_indices
from automol.graph.base._amchi import inchi_is_bad
# SMILES functions:
from automol.graph.base._smiles import smiles
# stereo functions:
# # core functions
from automol.graph.base._stereo import expand_stereo
from automol.graph.base._stereo import expand_stereo_with_priorities_and_amchis
# # stereo correction
from automol.graph.base._stereo import stereo_corrected_geometry
# functional groups code:
# # core functions
from automol.graph.base._func_group import FunctionalGroup
from automol.graph.base._func_group import functional_group_dct
from automol.graph.base._func_group import functional_group_count_dct
from automol.graph.base._func_group import ring_substituents
# # finders for overaching types
from automol.graph.base._func_group import is_hydrocarbon_species
from automol.graph.base._func_group import is_radical_species
# # finders for reactive sites and groups
from automol.graph.base._func_group import alkene_sites
from automol.graph.base._func_group import alkyne_sites
from automol.graph.base._func_group import alcohol_groups
from automol.graph.base._func_group import peroxy_groups
from automol.graph.base._func_group import hydroperoxy_groups
from automol.graph.base._func_group import ether_groups
from automol.graph.base._func_group import cyclic_ether_groups
from automol.graph.base._func_group import aldehyde_groups
from automol.graph.base._func_group import ketone_groups
from automol.graph.base._func_group import ester_groups
from automol.graph.base._func_group import carboxylic_acid_groups
from automol.graph.base._func_group import amide_groups
from automol.graph.base._func_group import nitro_groups
from automol.graph.base._func_group import halide_groups
from automol.graph.base._func_group import thiol_groups
from automol.graph.base._func_group import methyl_groups
from automol.graph.base._func_group import radical_dissociation_products
# # helper functions
from automol.graph.base._func_group import bonds_of_type
from automol.graph.base._func_group import bonds_of_order
from automol.graph.base._func_group import two_bond_idxs
from automol.graph.base._func_group import neighbors_of_type
from automol.graph.base._func_group import radicals_of_type
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
    'atom_implicit_hydrogens',
    'atom_stereo_parities',
    'bond_stereo_parities',
    'stereo_parities',
    # # setters
    'set_atom_symbols',
    'set_bond_orders',
    'set_atom_implicit_hydrogens',
    'set_atom_stereo_parities',
    'set_bond_stereo_parities',
    'set_stereo_parities',
    # # I/O
    'string',
    'yaml_dictionary',
    'from_string',
    'from_yaml_dictionary',
    'from_old_yaml_dictionary',
    # # conversions
    'frozen',
    'formula',
    # # properties
    'atom_count',
    'electron_count',
    'atom_stereo_keys',
    'bond_stereo_keys',
    'stereo_keys',
    'has_stereo',
    'has_atom_stereo',
    'has_bond_stereo',
    'has_pi_bonds',
    'is_ts_graph',
    'atomic_numbers',
    'mass_numbers',
    'atomic_valences',
    'atom_lone_pairs',
    'atom_van_der_waals_radius',
    'atom_bond_counts',
    'atom_unpaired_electrons',
    'bond_unpaired_electrons',
    'tetrahedral_atom_keys',
    'maximum_spin_multiplicity',
    'possible_spin_multiplicities',
    'atom_symbol_keys',
    'backbone_keys',
    'backbone_bond_keys',
    'backbone_hydrogen_keys',
    'nonbackbone_hydrogen_keys',
    'atom_backbone_hydrogen_keys',
    'atom_nonbackbone_hydrogen_keys',
    'terminal_atom_keys',
    'unsaturated_atom_keys',
    'unsaturated_bond_keys',
    'lone_pair_atom_keys',
    'angle_keys',
    # # relabeling and changing keys
    'relabel',
    'standard_keys',
    'standard_keys_for_sequence',
    'relabel_for_zmatrix',
    'relabel_for_geometry',
    # # add/remove/insert/without
    'add_atoms',
    'add_bonds',
    'remove_atoms',
    'remove_bonds',
    'change_implicit_hydrogens',
    'add_atom_explicit_hydrogens',
    'add_bonded_atom',
    'shift_insert_bonded_atom',
    'add_dummy_atoms',
    'shift_insert_dummy_atoms',
    'shift_remove_dummy_atoms',
    'without_pi_bonds',
    'without_dummy_atoms',
    'ts_reverse',
    'ts_reagents_graph_without_stereo',
    'without_bonds_by_orders',
    'without_stereo',
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
    'atom_neighbor_atom_keys',
    'atoms_neighbor_atom_keys',
    'atom_sorted_neighbor_atom_keys',
    'atoms_sorted_neighbor_atom_keys',
    'atom_bond_keys',
    'atoms_bond_keys',
    'dummy_atoms_neighbor_atom_key',
    'bonds_neighbor_atom_keys',
    'bonds_neighbor_bond_keys',
    # algorithm functions:
    # # isomorphisms and equivalence
    'isomorphism',
    'isomorphic',
    'unique',
    'sequence_isomorphism',
    'subgraph_isomorphism',
    'equivalent_atoms',
    'equivalent_bonds',
    'are_equivalent_atoms',
    'are_equivalent_bonds',
    'atom_equivalence_class_reps',
    'bond_equivalence_class_reps',
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
    'branches',
    'branch_dict',
    'branch_atom_keys',
    'branch',
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
    # kekule functions:
    # # core functions
    'kekule',
    'kekules',
    'kekule_bond_orders',
    'kekules_bond_orders',
    'kekules_bond_orders_collated',
    'kekules_bond_orders_averaged',
    # # derived properties
    'linear_atom_keys',
    'linear_segments_atom_keys',
    'atom_hybridizations',
    'atom_hybridizations_from_kekule',
    'radical_atom_keys',
    'radical_atom_keys_from_kekule',
    'nonresonant_radical_atom_keys',
    'vinyl_radical_atom_keys',
    'sigma_radical_atom_keys',
    'has_separated_radical_sites',
    'resonance_bond_stereo_keys',
    'vinyl_bond_stereo_keys',
    'has_resonance_bond_stereo',
    'has_vinyl_bond_stereo',
    'has_nonkekule_bond_stereo',
    'has_noninchi_stereo',
    'radical_groups',
    'radical_group_dct',
    'rigid_planar_bond_keys',
    'atom_centered_cumulene_keys',
    'bond_centered_cumulene_keys',
    # torsion/rotational bond functions:
    'rotational_bond_keys',
    'rotational_groups',
    'rotational_symmetry_number',
    # geometry functions:
    # # stereo parity evaluations
    'geometry_atom_parity',
    'geometry_bond_parity',
    'geometry_local_parity',
    'geometries_have_matching_parities',
    'geometries_parity_mismatches',
    # # corrections
    'linear_vinyl_corrected_geometry',
    'geometry_rotate_bond',
    'geometry_dihedrals_near_value',
    'perturb_geometry_planar_dihedrals',
    # canonicalization functions:
    # # canonical key functions
    'canonical_enantiomer',
    'canonical_enantiomer_with_keys',
    'canonical',
    'canonical_keys',
    # # canonical stereo functions
    'stereogenic_atom_keys',
    'stereogenic_bond_keys',
    'stereogenic_keys',
    'stereogenic_atom_keys_from_priorities',
    'stereogenic_bond_keys_from_priorities',
    'reflect',
    'reflect_local_stereo',
    'to_local_stereo',
    'from_local_stereo',
    'set_stereo_from_geometry',
    # # symmetry class functions
    'canonical_priorities',
    'calculate_priorities_and_assign_stereo',
    # # parity evaluators
    'parity_evaluator_from_geometry_',
    'parity_evaluator_read_canonical_',
    'parity_evaluator_flip_local_',
    # AMChI functions:
    'amchi',
    'amchi_with_indices',
    'inchi_is_bad',
    # SMILES functions:
    'smiles',
    # stereo functions:
    # # core functions
    'expand_stereo',
    'expand_stereo_with_priorities_and_amchis',
    # # stereo correction
    'stereo_corrected_geometry',
    # functional groups code:
    # # core functions
    'FunctionalGroup',
    'functional_group_dct',
    'functional_group_count_dct',
    'ring_substituents',
    # # finders for overaching types
    'is_hydrocarbon_species',
    'is_radical_species',
    # # finders for reactive sites and groups
    'alkene_sites',
    'alkyne_sites',
    'alcohol_groups',
    'peroxy_groups',
    'hydroperoxy_groups',
    'ether_groups',
    'cyclic_ether_groups',
    'aldehyde_groups',
    'ketone_groups',
    'ester_groups',
    'carboxylic_acid_groups',
    'amide_groups',
    'nitro_groups',
    'halide_groups',
    'thiol_groups',
    'methyl_groups',
    'radical_dissociation_products',
    # # helper functions
    'bonds_of_type',
    'bonds_of_order',
    'two_bond_idxs',
    'neighbors_of_type',
    'radicals_of_type',
    # TS graph submodule:
    'ts_graph',
    'ts',
]
