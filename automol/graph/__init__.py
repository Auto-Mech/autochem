""" molecular graphs

Import hierarchy:
    vmat        no dependencies
    embed       no dependencies
    _conv       dependencies: vmat

Level 2 graph functions belong in graph/base/*.py. These have no dependencies
on automol.extern or on other automol types (geom, zmat, etc.)

Level 4 graph functions, which do depend on these things, belong in graph/*.py

"""

# L2
# core functions:
# # constructors
# L4
# embedding functions:
from automol.graph._0embed import (
    clean_geometry,
    embed_geometry,
    geometry_matches,
    zmatrix_matches,
)

# conversion functions:
# # RMG conversions
from automol.graph._1rmg import (
    RMG_ADJACENCY_LIST,
    from_parsed_rmg_adjacency_list,
    from_rmg_adjacency_list,
)

# # conversions
from automol.graph._2conv import (
    chi,
    display,
    display_reaction,
    geometry,
    inchi,
    ipywidget,
    perturb_geometry_planar_dihedrals,
    rdkit_molecule,
    rdkit_reaction,
    svg_string,
    ts_geometry_from_reactants,
)

# submodules:
from automol.graph.base import ts, vmat

# # getters
# # setters
# # I/O
# # conversions
# # sorting
# # properties
# # relabeling and changing keys
# # add/remove/insert/without
# # unions
# # subgraphs and neighborhoods
from automol.graph.base._00core import (
    add_atom_explicit_hydrogens,
    add_atoms,
    add_bonded_atom,
    add_bonds,
    angle_keys,
    apply_zmatrix_conversion,
    argsort_by_size,
    atom_backbone_hydrogen_keys,
    atom_bond_counts,
    atom_bond_keys,
    atom_count,
    atom_electron_pairs,
    atom_implicit_hydrogens,
    atom_keys,
    atom_lone_pairs,
    atom_neighbor_atom_key,
    atom_neighbor_atom_keys,
    atom_neighborhood,
    atom_neighborhoods,
    atom_nonbackbone_hydrogen_keys,
    atom_sorted_neighbor_atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    atom_symbol_keys,
    atom_symbols,
    atom_unpaired_electrons,
    atom_van_der_waals_radius,
    atomic_numbers,
    atomic_valences,
    atoms,
    atoms_bond_keys,
    atoms_from_data,
    atoms_neighbor_atom_keys,
    atoms_sorted_neighbor_atom_keys,
    backbone_bond_keys,
    backbone_hydrogen_keys,
    backbone_keys,
    bond_induced_subgraph,
    bond_keys,
    bond_neighbor_atom_keys,
    bond_neighbor_bond_keys,
    bond_neighborhood,
    bond_neighborhoods,
    bond_orders,
    bond_stereo_keys,
    bond_stereo_parities,
    bond_unpaired_electrons,
    bonds,
    bonds_from_data,
    bonds_neighbor_atom_keys,
    bonds_neighbor_bond_keys,
    change_implicit_hydrogens,
    count,
    covalent_radii,
    dummy_source_dict,
    electron_count,
    explicit,
    formula,
    from_atoms_and_bonds,
    from_data,
    from_old_yaml_data,
    from_string,
    from_yaml_data,
    frozen,
    has_atom_stereo,
    has_bond_stereo,
    has_dummy_atoms,
    has_pi_bonds,
    has_stereo,
    implicit,
    invert_atom_stereo_parities,
    is_ts_graph,
    local_stereo_priorities,
    lone_pair_atom_keys,
    mass_numbers,
    maximum_spin_multiplicity,
    nonbackbone_hydrogen_keys,
    possible_spin_multiplicities,
    relabel,
    remove_atoms,
    remove_bonds,
    set_atom_implicit_hydrogens,
    set_atom_stereo_parities,
    set_atom_symbols,
    set_bond_orders,
    set_bond_stereo_parities,
    set_stereo_parities,
    sort_by_size,
    standard_keys,
    standard_keys_for_sequence,
    stereo_keys,
    stereo_parities,
    string,
    subgraph,
    symbols,
    terminal_atom_keys,
    terminal_atom_neighbors,
    tetrahedral_atom_keys,
    tetrahedral_atoms,
    ts_graph,
    ts_products_graph_without_stereo,
    ts_reactants_graph_without_stereo,
    ts_reagents_graph_without_stereo,
    ts_reagents_graphs_without_stereo,
    ts_reverse,
    undo_zmatrix_conversion,
    union,
    union_from_sequence,
    unsaturated_atom_keys,
    unsaturated_bond_keys,
    van_der_waals_radii,
    vinyl_radical_candidates,
    with_explicit_stereo_hydrogens,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_pi_bonds,
    without_reacting_bonds,
    without_stereo,
    yaml_data,
    zmatrix_conversion_info,
)

# algorithm functions:
# # isomorphisms and equivalence
# # algorithms
# # branches and groups
# # rings
from automol.graph.base._02algo import (
    are_equivalent_atoms,
    are_equivalent_bonds,
    atom_equivalence_class_reps,
    atom_longest_chain,
    atom_longest_chains,
    atom_shortest_paths,
    bond_equivalence_class_reps,
    branch,
    branch_atom_keys,
    branch_dict,
    branches,
    connected_components,
    connected_components_atom_keys,
    equivalent_atoms,
    equivalent_bonds,
    is_branched,
    is_connected,
    is_ring_key_sequence,
    is_ring_system,
    isomorphic,
    isomorphism,
    longest_chain,
    ring_arc_complement_atom_keys,
    ring_atom_chirality,
    ring_system_decomposed_atom_keys,
    ring_systems,
    ring_systems_atom_keys,
    ring_systems_bond_keys,
    ring_systems_decomposed_atom_keys,
    rings,
    rings_atom_keys,
    rings_bond_keys,
    sequence_isomorphism,
    shortest_path_between_atoms,
    shortest_path_between_groups,
    sorted_ring_atom_keys,
    sorted_ring_atom_keys_from_bond_keys,
    spiro_atom_keys,
    spiro_atoms_grouped_neighbor_keys,
    subgraph_isomorphism,
    unique,
)

# kekule functions:
# # core functions
# # derived properties
from automol.graph.base._03kekule import (
    atom_centered_cumulene_keys,
    atom_hybridizations,
    atom_hybridizations_from_kekule,
    bad_stereo_bond_keys_from_kekule,
    bond_centered_cumulene_keys,
    good_stereo_bond_keys_from_kekule,
    has_noninchi_stereo,
    has_nonkekule_bond_stereo,
    has_resonance_bond_stereo,
    has_separated_radical_sites,
    has_vinyl_bond_stereo,
    kekule,
    kekule_bond_orders,
    kekules,
    kekules_bond_orders,
    kekules_bond_orders_averaged,
    kekules_bond_orders_collated,
    linear_atom_keys,
    linear_segment_cap_keys,
    linear_segments_atom_keys,
    nonresonant_radical_atom_keys,
    possible_rigid_planar_bond_keys,
    radical_atom_keys,
    radical_atom_keys_from_kekule,
    radical_group_dct,
    radical_groups,
    resonance_bond_stereo_keys,
    rigid_planar_bond_keys,
    rigid_planar_bonds,
    sigma_radical_atom_bond_keys,
    sigma_radical_atom_keys,
    strict_rigid_planar_bond_keys,
    unneeded_dummy_atom_keys,
    vinyl_bond_stereo_keys,
    vinyl_radical_atom_bond_keys,
    vinyl_radical_atom_keys,
)

# geometry functions:
# # stereo parity evaluations
# # parity evaluators
from automol.graph.base._05stereo import (
    geometry_atom_parity,
    geometry_bond_parity,
    parity_evaluator_flip_from_graph,
    parity_evaluator_measure_from_geometry_,
    parity_evaluator_read_from_graph,
    stereocenter_candidates,
    unassigned_stereocenter_keys_from_candidates,
)

# structural heuristics:
from automol.graph.base._06heur import (
    heuristic_bond_angle,
    heuristic_bond_distance,
    heuristic_bond_distance_limit,
    rotational_bond_keys,
    rotational_coordinates,
    rotational_groups,
    rotational_segment_keys,
    rotational_symmetry_number,
)

# # corrections
from automol.graph.base._07geom import (
    geometries_parity_mismatches,
    geometry_correct_linear_vinyls,
    geometry_dihedrals_near_value,
    geometry_local_parity,
    geometry_pseudorotate_atom,
    geometry_rotate_bond,
)

# canonicalization functions:
# # canonical key functions
# # symmetry class functions
from automol.graph.base._08canon import (
    calculate_stereo,
    canonical,
    canonical_keys,
    canonical_priorities,
)

# AMChI functions:
from automol.graph.base._09amchi import amchi, amchi_with_numbers, inchi_is_bad

# SMILES functions:
from automol.graph.base._10smiles import smiles

# # canonical stereo functions
# stereo functions:
# # core functions
# # stereo correction
from automol.graph.base._11stereo import (
    expand_reaction_stereo,
    expand_stereo,
    from_local_stereo,
    has_fleeting_atom_or_bond_stereo,
    reflect,
    set_stereo_from_geometry,
    stereo_corrected_geometry,
    to_local_stereo,
    unassigned_stereocenter_keys,
)

# functional groups code:
# # core functions
# # finders for overaching types
# # finders for reactive sites and groups
# # helper functions
from automol.graph.base._func_group import (
    FunctionalGroup,
    alcohol_groups,
    aldehyde_groups,
    alkene_sites,
    alkyne_sites,
    amide_groups,
    bonds_of_order,
    bonds_of_type,
    carboxylic_acid_groups,
    cyclic_ether_groups,
    ester_groups,
    ether_groups,
    functional_group_count_dct,
    functional_group_dct,
    halide_groups,
    hydroperoxy_groups,
    is_hydrocarbon_species,
    is_radical_species,
    ketone_groups,
    methyl_groups,
    neighbors_of_type,
    nitro_groups,
    peroxy_groups,
    radical_dissociation_products,
    radicals_of_type,
    ring_substituents,
    thiol_groups,
    two_bond_idxs,
)

__all__ = [
    # L2
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
    "tetrahedral_atoms",
    "tetrahedral_atom_keys",
    "vinyl_radical_candidates",
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
    "angle_keys",
    # # relabeling and changing keys
    "relabel",
    "standard_keys",
    "standard_keys_for_sequence",
    "zmatrix_conversion_info",
    "apply_zmatrix_conversion",
    "undo_zmatrix_conversion",
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
    "spiro_atoms_grouped_neighbor_keys",
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
    "radical_dissociation_products",
    # # helper functions
    "bonds_of_type",
    "bonds_of_order",
    "two_bond_idxs",
    "neighbors_of_type",
    "radicals_of_type",
    # TS graph submodule:
    "ts_graph",
    "ts",
    "vmat",
    # L4
    # embedding functions:
    "embed_geometry",
    "clean_geometry",
    "geometry_matches",
    "zmatrix_matches",
    # conversion functions:
    # # RMG conversions
    "RMG_ADJACENCY_LIST",
    "from_rmg_adjacency_list",
    "from_parsed_rmg_adjacency_list",
    # # conversions
    "geometry",
    "ts_geometry_from_reactants",
    "inchi",
    "chi",
    "rdkit_molecule",
    "rdkit_reaction",
    "display",
    "display_reaction",
    "svg_string",
    "ipywidget",
    "perturb_geometry_planar_dihedrals",
]
