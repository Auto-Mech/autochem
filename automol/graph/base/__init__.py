""" Level 2 graph functions (no dependencies on extern or other types)

Import hierarchy:
    _core       no dependencies
    _networkx   dependencies: _core
    _algo       dependencies: _core, _networkx
    _kekule     dependencies: _core, _networkx, _algo
    _heur       dependencies: _core, _networkx, _algo, _kekule
    _geom       dependencies: _core, _networkx, _algo, _kekule
    _canon      dependencies: _core, _networkx, _algo, _kekule
    _amchi      dependencies: _core, _networkx, _algo, _canon, _kekule
    _smiles     dependencies: _core, _networkx, _algo, _canon, _kekule
    _stereo     dependencies: _core, _networkx, _algo, _canon, _kekule, _amchi
    _func_group dependencies: _core, _networkx, _algo, _kekule, _stereo
    ts          dependencies: _core, _networkx, _algo, _kekule, _stereo
    vmat        dependencies: ts

Each next submodule in the hierarchy may depend on the ones before it, but
**not** the ones after. This avoids circular dependencies.

Please obey this hierarchy strictly and raise an issue to the group if
something needs to be updated or altered in it.
"""

# core functions:
# # constructors
from ._00core import from_data
from ._00core import atoms_from_data
from ._00core import bonds_from_data
from ._00core import from_atoms_and_bonds

# # getters
from ._00core import atoms
from ._00core import bonds
from ._00core import atom_keys
from ._00core import bond_keys
from ._00core import atom_symbols
from ._00core import bond_orders
from ._00core import atom_implicit_hydrogens
from ._00core import atom_stereo_parities
from ._00core import bond_stereo_parities
from ._00core import stereo_parities
from ._00core import ts_graph

# # setters
from ._00core import set_atom_symbols
from ._00core import set_bond_orders
from ._00core import set_atom_implicit_hydrogens
from ._00core import set_atom_stereo_parities
from ._00core import set_bond_stereo_parities
from ._00core import set_stereo_parities

# # I/O
from ._00core import string
from ._00core import yaml_data
from ._00core import from_string
from ._00core import from_yaml_data
from ._00core import from_old_yaml_data

# # conversions
from ._00core import frozen
from ._00core import formula
from ._00core import symbols

# # sorting
from ._00core import argsort_by_size
from ._00core import sort_by_size

# # properties
from ._00core import atom_count
from ._00core import count
from ._00core import electron_count
from ._00core import atom_stereo_keys
from ._00core import bond_stereo_keys
from ._00core import stereo_keys
from ._00core import has_stereo
from ._00core import has_atom_stereo
from ._00core import has_bond_stereo
from ._00core import has_dummy_atoms
from ._00core import has_pi_bonds
from ._00core import is_ts_graph
from ._00core import atomic_numbers
from ._00core import mass_numbers
from ._00core import van_der_waals_radii
from ._00core import covalent_radii
from ._00core import atomic_valences
from ._00core import atom_lone_pairs
from ._00core import atom_electron_pairs
from ._00core import atom_van_der_waals_radius
from ._00core import atom_bond_counts
from ._00core import atom_unpaired_electrons
from ._00core import bond_unpaired_electrons
from ._00core import atom_hypervalencies
from ._00core import tetrahedral_atoms
from ._00core import tetrahedral_atom_keys
from ._00core import vinyl_radical_bond_candidates
from ._00core import maximum_spin_multiplicity
from ._00core import possible_spin_multiplicities
from ._00core import atom_symbol_keys
from ._00core import backbone_keys
from ._00core import backbone_bond_keys
from ._00core import backbone_hydrogen_keys
from ._00core import nonbackbone_hydrogen_keys
from ._00core import atom_backbone_hydrogen_keys
from ._00core import atom_nonbackbone_hydrogen_keys
from ._00core import terminal_atom_keys
from ._00core import terminal_atom_neighbors
from ._00core import unsaturated_atom_keys
from ._00core import unsaturated_bond_keys
from ._00core import lone_pair_atom_keys
from ._00core import distance_keys
from ._00core import central_angle_keys
from ._00core import dihedral_angle_keys

# # relabeling and changing keys
from ._00core import relabel
from ._00core import standard_keys
from ._00core import standard_keys_for_sequence
from ._00core import zmatrix_conversion_info
from ._00core import apply_zmatrix_conversion
from ._00core import undo_zmatrix_conversion
from ._00core import align_with_geometry

# # add/remove/insert/without
from ._00core import add_atoms
from ._00core import add_bonds
from ._00core import remove_atoms
from ._00core import remove_bonds
from ._00core import change_implicit_hydrogens
from ._00core import add_atom_explicit_hydrogens
from ._00core import add_bonded_atom
from ._00core import without_pi_bonds
from ._00core import without_reacting_bonds
from ._00core import without_dummy_atoms
from ._00core import ts_reverse
from ._00core import ts_reactants_graph_without_stereo
from ._00core import ts_products_graph_without_stereo
from ._00core import ts_reagents_graphs_without_stereo
from ._00core import ts_reagents_graph_without_stereo
from ._00core import without_bonds_by_orders
from ._00core import without_stereo
from ._00core import with_explicit_stereo_hydrogens
from ._00core import explicit
from ._00core import implicit

# # unions
from ._00core import union
from ._00core import union_from_sequence

# # subgraphs and neighborhoods
from ._00core import subgraph
from ._00core import bond_induced_subgraph
from ._00core import atom_neighborhood
from ._00core import atom_neighborhoods
from ._00core import bond_neighborhood
from ._00core import bond_neighborhoods
from ._00core import atom_neighbor_atom_key
from ._00core import atom_neighbor_atom_keys
from ._00core import atoms_neighbor_atom_keys
from ._00core import atom_sorted_neighbor_atom_keys
from ._00core import local_stereo_priorities
from ._00core import atoms_sorted_neighbor_atom_keys
from ._00core import atom_bond_keys
from ._00core import atoms_bond_keys
from ._00core import dummy_source_dict
from ._00core import bond_neighbor_atom_keys
from ._00core import bond_neighbor_bond_keys
from ._00core import bonds_neighbor_atom_keys
from ._00core import bonds_neighbor_bond_keys

# algorithm functions:
# # isomorphisms and equivalence
from ._02algo import isomorphism
from ._02algo import isomorphic
from ._02algo import unique
from ._02algo import sequence_isomorphism
from ._02algo import subgraph_isomorphism
from ._02algo import equivalent_atoms
from ._02algo import equivalent_bonds
from ._02algo import are_equivalent_atoms
from ._02algo import are_equivalent_bonds
from ._02algo import atom_equivalence_class_reps
from ._02algo import bond_equivalence_class_reps

# # algorithms
from ._02algo import connected_components
from ._02algo import connected_components_atom_keys
from ._02algo import is_connected
from ._02algo import dfs_
from ._02algo import dfs_atom_keys
from ._02algo import dfs_bond_keys
from ._02algo import dfs_children
from ._02algo import dfs_parents
from ._02algo import dfs_missing_bond_keys
from ._02algo import dfs_missing_children
from ._02algo import atom_shortest_paths
from ._02algo import shortest_path_between_atoms
from ._02algo import shortest_path_between_groups
from ._02algo import atom_longest_chains
from ._02algo import atom_longest_chain
from ._02algo import longest_chain

# # branches and groups
from ._02algo import ring_atom_chirality
from ._02algo import branches
from ._02algo import branch_dict
from ._02algo import branch_atom_keys
from ._02algo import branch
from ._02algo import is_branched

# # rings
from ._02algo import rings
from ._02algo import rings_atom_keys
from ._02algo import rings_bond_keys
from ._02algo import sorted_ring_atom_keys
from ._02algo import sorted_ring_atom_keys_from_bond_keys
from ._02algo import is_ring_key_sequence
from ._02algo import ring_arc_complement_atom_keys
from ._02algo import ring_systems
from ._02algo import ring_systems_atom_keys
from ._02algo import ring_systems_bond_keys
from ._02algo import spiros
from ._02algo import spiro_atom_keys
from ._02algo import is_ring_system
from ._02algo import ring_system_decomposed_atom_keys
from ._02algo import ring_systems_decomposed_atom_keys

# kekule functions:
# # core functions
from ._03kekule import kekule
from ._03kekule import kekules
from ._03kekule import kekule_bond_orders
from ._03kekule import kekules_bond_orders
from ._03kekule import kekules_bond_orders_collated
from ._03kekule import kekules_bond_orders_averaged

# # derived properties
from ._03kekule import linear_atom_keys
from ._03kekule import linear_segment_cap_keys
from ._03kekule import linear_segments_atom_keys
from ._03kekule import unneeded_dummy_atom_keys
from ._03kekule import atom_hybridizations
from ._03kekule import atom_hybridizations_from_kekule
from ._03kekule import bad_stereo_bond_keys_from_kekule
from ._03kekule import good_stereo_bond_keys_from_kekule
from ._03kekule import radical_atom_keys
from ._03kekule import radical_atom_keys_from_kekule
from ._03kekule import nonresonant_radical_atom_keys
from ._03kekule import vinyl_radical_atom_bond_keys
from ._03kekule import sigma_radical_atom_bond_keys
from ._03kekule import vinyl_radical_atom_keys
from ._03kekule import sigma_radical_atom_keys
from ._03kekule import has_separated_radical_sites
from ._03kekule import addition_atom_keys
from ._03kekule import beta_scission_bond_keys
from ._03kekule import beta_scission_bond_keys_from_kekule
from ._03kekule import resonance_bond_stereo_keys
from ._03kekule import vinyl_bond_stereo_keys
from ._03kekule import has_resonance_bond_stereo
from ._03kekule import has_vinyl_bond_stereo
from ._03kekule import has_nonkekule_bond_stereo
from ._03kekule import has_noninchi_stereo
from ._03kekule import radical_groups
from ._03kekule import radical_group_dct
from ._03kekule import rigid_planar_bonds
from ._03kekule import rigid_planar_bond_keys
from ._03kekule import strict_rigid_planar_bond_keys
from ._03kekule import possible_rigid_planar_bond_keys
from ._05stereo import stereocenter_candidates
from ._03kekule import atom_centered_cumulene_keys
from ._03kekule import bond_centered_cumulene_keys

# structural heuristics:
from ._06heur import heuristic_bond_distance
from ._06heur import heuristic_bond_distance_limit
from ._06heur import heuristic_bond_angle
from ._06heur import rotational_bond_keys
from ._06heur import rotational_segment_keys
from ._06heur import rotational_coordinates
from ._06heur import rotational_groups
from ._06heur import rotational_symmetry_number

# geometry functions:
# # stereo parity evaluations
from ._05stereo import geometry_atom_parity
from ._05stereo import geometry_bond_parity
from ._07geom import geometry_local_parity
from ._07geom import geometries_parity_mismatches

# # corrections
from ._07geom import geometry_correct_linear_vinyls
from ._11stereo import geometry_pseudorotate_atom
from ._07geom import geometry_rotate_bond
from ._07geom import geometry_dihedrals_near_value

# canonicalization functions:
# # canonical key functions
from ._08canon import canonical
from ._08canon import canonical_keys
from ._08canon import smiles_graph

# # canonical stereo functions
from ._11stereo import unassigned_stereocenter_keys
from ._05stereo import unassigned_stereocenter_keys_from_candidates
from ._11stereo import reflect
from ._00core import invert_atom_stereo_parities
from ._08canon import to_local_stereo
from ._08canon import from_local_stereo
from ._11stereo import set_stereo_from_geometry

# # symmetry class functions
from ._08canon import canonical_priorities
from ._08canon import calculate_stereo

# # parity evaluators
from ._05stereo import parity_evaluator_measure_from_geometry_
from ._05stereo import parity_evaluator_read_from_graph
from ._05stereo import parity_evaluator_flip_from_graph

# AMChI functions:
from ._09amchi import amchi
from ._09amchi import amchi_with_numbers
from ._09amchi import inchi_is_bad

# SMILES functions:
from ._10smiles import smiles

# stereo functions:
# # core functions
from ._11stereo import expand_stereo
from ._11stereo import expand_reaction_stereo

# # stereo correction
from ._11stereo import has_fleeting_atom_or_bond_stereo
from ._11stereo import stereo_corrected_geometry

# functional groups code:
# # core functions
from ._func_group import FunctionalGroup
from ._func_group import functional_group_dct
from ._func_group import functional_group_count_dct
from ._func_group import ring_substituents

# # finders for overaching types
from ._func_group import is_hydrocarbon_species
from ._func_group import is_radical_species

# # finders for reactive sites and groups
from ._func_group import alkene_sites
from ._func_group import alkyne_sites
from ._func_group import alcohol_groups
from ._func_group import peroxy_groups
from ._func_group import hydroperoxy_groups
from ._func_group import ether_groups
from ._func_group import cyclic_ether_groups
from ._func_group import aldehyde_groups
from ._func_group import ketone_groups
from ._func_group import ester_groups
from ._func_group import carboxylic_acid_groups
from ._func_group import amide_groups
from ._func_group import nitro_groups
from ._func_group import halide_groups
from ._func_group import thiol_groups
from ._func_group import methyl_groups
from ._func_group import aromatic_groups
from ._func_group import phenyl_groups
from ._func_group import benzyl_groups
from ._func_group import cyclopentadienyl_groups
from ._func_group import radical_dissociation_products

# # helper functions
from ._func_group import bonds_of_type
from ._func_group import bonds_of_order
from ._func_group import two_bond_idxs
from ._func_group import neighbors_of_type
from ._func_group import radicals_of_type

# submodules:
from . import ts
from . import vmat


__all__ = [
    # core functions:
    # # constructors
    "from_data",
    "atoms_from_data",
    "bonds_from_data",
    "from_atoms_and_bonds",
    # # getters
    "atoms",
    "bonds",
    "atom_keys",
    "bond_keys",
    "atom_symbols",
    "bond_orders",
    "atom_implicit_hydrogens",
    "atom_stereo_parities",
    "bond_stereo_parities",
    "stereo_parities",
    "ts_graph",
    # # setters
    "set_atom_symbols",
    "set_bond_orders",
    "set_atom_implicit_hydrogens",
    "set_atom_stereo_parities",
    "set_bond_stereo_parities",
    "set_stereo_parities",
    # # I/O
    "string",
    "yaml_data",
    "from_string",
    "from_yaml_data",
    "from_old_yaml_data",
    # # conversions
    "frozen",
    "formula",
    "symbols",
    # # sorting
    "argsort_by_size",
    "sort_by_size",
    # # properties
    "atom_count",
    "count",
    "electron_count",
    "atom_stereo_keys",
    "bond_stereo_keys",
    "stereo_keys",
    "has_stereo",
    "has_atom_stereo",
    "has_bond_stereo",
    "has_dummy_atoms",
    "has_pi_bonds",
    "is_ts_graph",
    "atomic_numbers",
    "mass_numbers",
    "van_der_waals_radii",
    "covalent_radii",
    "atomic_valences",
    "atom_lone_pairs",
    "atom_electron_pairs",
    "atom_van_der_waals_radius",
    "atom_bond_counts",
    "atom_unpaired_electrons",
    "bond_unpaired_electrons",
    "atom_hypervalencies",
    "tetrahedral_atoms",
    "tetrahedral_atom_keys",
    "vinyl_radical_bond_candidates",
    "maximum_spin_multiplicity",
    "possible_spin_multiplicities",
    "atom_symbol_keys",
    "backbone_keys",
    "backbone_bond_keys",
    "backbone_hydrogen_keys",
    "nonbackbone_hydrogen_keys",
    "atom_backbone_hydrogen_keys",
    "atom_nonbackbone_hydrogen_keys",
    "terminal_atom_keys",
    "terminal_atom_neighbors",
    "unsaturated_atom_keys",
    "unsaturated_bond_keys",
    "lone_pair_atom_keys",
    "distance_keys",
    "central_angle_keys",
    "dihedral_angle_keys",
    # # relabeling and changing keys
    "relabel",
    "standard_keys",
    "standard_keys_for_sequence",
    "zmatrix_conversion_info",
    "apply_zmatrix_conversion",
    "undo_zmatrix_conversion",
    "align_with_geometry",
    # # add/remove/insert/without
    "add_atoms",
    "add_bonds",
    "remove_atoms",
    "remove_bonds",
    "change_implicit_hydrogens",
    "add_atom_explicit_hydrogens",
    "add_bonded_atom",
    "without_pi_bonds",
    "without_reacting_bonds",
    "without_dummy_atoms",
    "ts_reverse",
    "ts_reactants_graph_without_stereo",
    "ts_products_graph_without_stereo",
    "ts_reagents_graphs_without_stereo",
    "ts_reagents_graph_without_stereo",
    "without_bonds_by_orders",
    "without_stereo",
    "with_explicit_stereo_hydrogens",
    "explicit",
    "implicit",
    # # unions
    "union",
    "union_from_sequence",
    # # subgraphs and neighborhoods
    "subgraph",
    "bond_induced_subgraph",
    "atom_neighborhood",
    "atom_neighborhoods",
    "bond_neighborhood",
    "bond_neighborhoods",
    "atom_neighbor_atom_key",
    "atom_neighbor_atom_keys",
    "atoms_neighbor_atom_keys",
    "atom_sorted_neighbor_atom_keys",
    "local_stereo_priorities",
    "atoms_sorted_neighbor_atom_keys",
    "atom_bond_keys",
    "atoms_bond_keys",
    "dummy_source_dict",
    "bond_neighbor_atom_keys",
    "bond_neighbor_bond_keys",
    "bonds_neighbor_atom_keys",
    "bonds_neighbor_bond_keys",
    # algorithm functions:
    # # isomorphisms and equivalence
    "isomorphism",
    "isomorphic",
    "unique",
    "sequence_isomorphism",
    "subgraph_isomorphism",
    "equivalent_atoms",
    "equivalent_bonds",
    "are_equivalent_atoms",
    "are_equivalent_bonds",
    "atom_equivalence_class_reps",
    "bond_equivalence_class_reps",
    # # algorithms
    "connected_components",
    "connected_components_atom_keys",
    "is_connected",
    "dfs_",
    "dfs_atom_keys",
    "dfs_bond_keys",
    "dfs_children",
    "dfs_parents",
    "dfs_missing_bond_keys",
    "dfs_missing_children",
    "atom_shortest_paths",
    "shortest_path_between_atoms",
    "shortest_path_between_groups",
    "atom_longest_chains",
    "atom_longest_chain",
    "longest_chain",
    # # branches and groups
    "ring_atom_chirality",
    "branches",
    "branch_dict",
    "branch_atom_keys",
    "branch",
    "is_branched",
    # # rings
    "rings",
    "rings_atom_keys",
    "rings_bond_keys",
    "sorted_ring_atom_keys",
    "sorted_ring_atom_keys_from_bond_keys",
    "is_ring_key_sequence",
    "ring_arc_complement_atom_keys",
    "ring_systems",
    "ring_systems_atom_keys",
    "ring_systems_bond_keys",
    "spiros",
    "spiro_atom_keys",
    "is_ring_system",
    "ring_system_decomposed_atom_keys",
    "ring_systems_decomposed_atom_keys",
    # kekule functions:
    # # core functions
    "kekule",
    "kekules",
    "kekule_bond_orders",
    "kekules_bond_orders",
    "kekules_bond_orders_collated",
    "kekules_bond_orders_averaged",
    # # derived properties
    "linear_atom_keys",
    "linear_segment_cap_keys",
    "linear_segments_atom_keys",
    "unneeded_dummy_atom_keys",
    "atom_hybridizations",
    "atom_hybridizations_from_kekule",
    "bad_stereo_bond_keys_from_kekule",
    "good_stereo_bond_keys_from_kekule",
    "radical_atom_keys",
    "radical_atom_keys_from_kekule",
    "nonresonant_radical_atom_keys",
    "vinyl_radical_atom_bond_keys",
    "sigma_radical_atom_bond_keys",
    "vinyl_radical_atom_keys",
    "sigma_radical_atom_keys",
    "has_separated_radical_sites",
    "addition_atom_keys",
    "beta_scission_bond_keys",
    "beta_scission_bond_keys_from_kekule",
    "resonance_bond_stereo_keys",
    "vinyl_bond_stereo_keys",
    "has_resonance_bond_stereo",
    "has_vinyl_bond_stereo",
    "has_nonkekule_bond_stereo",
    "has_noninchi_stereo",
    "radical_groups",
    "radical_group_dct",
    "rigid_planar_bonds",
    "rigid_planar_bond_keys",
    "strict_rigid_planar_bond_keys",
    "possible_rigid_planar_bond_keys",
    "stereocenter_candidates",
    "atom_centered_cumulene_keys",
    "bond_centered_cumulene_keys",
    # structural heuristics:
    "heuristic_bond_distance",
    "heuristic_bond_distance_limit",
    "heuristic_bond_angle",
    "rotational_bond_keys",
    "rotational_segment_keys",
    "rotational_coordinates",
    "rotational_groups",
    "rotational_symmetry_number",
    # geometry functions:
    # # stereo parity evaluations
    "geometry_atom_parity",
    "geometry_bond_parity",
    "geometry_local_parity",
    "geometries_parity_mismatches",
    # # corrections
    "geometry_correct_linear_vinyls",
    "geometry_pseudorotate_atom",
    "geometry_rotate_bond",
    "geometry_dihedrals_near_value",
    # canonicalization functions:
    # # canonical key functions
    "canonical",
    "canonical_keys",
    "smiles_graph",
    # # canonical stereo functions
    "unassigned_stereocenter_keys",
    "unassigned_stereocenter_keys_from_candidates",
    "reflect",
    "invert_atom_stereo_parities",
    "to_local_stereo",
    "from_local_stereo",
    "set_stereo_from_geometry",
    # # symmetry class functions
    "canonical_priorities",
    "calculate_stereo",
    # # parity evaluators
    "parity_evaluator_measure_from_geometry_",
    "parity_evaluator_read_from_graph",
    "parity_evaluator_flip_from_graph",
    # AMChI functions:
    "amchi",
    "amchi_with_numbers",
    "inchi_is_bad",
    # SMILES functions:
    "smiles",
    # stereo functions:
    # # core functions
    "expand_stereo",
    "expand_reaction_stereo",
    # # stereo correction
    "has_fleeting_atom_or_bond_stereo",
    "stereo_corrected_geometry",
    # functional groups code:
    # # core functions
    "FunctionalGroup",
    "functional_group_dct",
    "functional_group_count_dct",
    "ring_substituents",
    # # finders for overaching types
    "is_hydrocarbon_species",
    "is_radical_species",
    # # finders for reactive sites and groups
    "alkene_sites",
    "alkyne_sites",
    "alcohol_groups",
    "peroxy_groups",
    "hydroperoxy_groups",
    "ether_groups",
    "cyclic_ether_groups",
    "aldehyde_groups",
    "ketone_groups",
    "ester_groups",
    "carboxylic_acid_groups",
    "amide_groups",
    "nitro_groups",
    "halide_groups",
    "thiol_groups",
    "methyl_groups",
    "aromatic_groups",
    "phenyl_groups",
    "benzyl_groups",
    "cyclopentadienyl_groups",
    "radical_dissociation_products",
    # # helper functions
    "bonds_of_type",
    "bonds_of_order",
    "two_bond_idxs",
    "neighbors_of_type",
    "radicals_of_type",
    # submodules:
    "ts",
    "vmat",
]
