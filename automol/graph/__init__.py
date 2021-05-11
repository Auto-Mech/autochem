""" molecular graph
"""
# graphbase
from automol.graph._graph_dep import atoms
from automol.graph._graph_dep import bonds
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import bond_keys
from automol.graph._graph_dep import atom_symbols
from automol.graph._graph_dep import atom_implicit_hydrogen_valences
from automol.graph._graph_dep import atom_stereo_parities
from automol.graph._graph_dep import bond_orders
from automol.graph._graph_dep import bond_stereo_parities
from automol.graph._graph_dep import set_atom_implicit_hydrogen_valences
from automol.graph._graph_dep import set_atom_stereo_parities
from automol.graph._graph_dep import set_bond_orders
from automol.graph._graph_dep import set_bond_stereo_parities
from automol.graph._graph_dep import string
from automol.graph._graph_dep import yaml_dictionary
from automol.graph._graph_dep import relabel
# getters
from automol.graph._graph_dep import without_bond_orders
from automol.graph._graph_dep import without_stereo_parities
from automol.graph._graph_dep import atom_explicit_hydrogen_keys
from automol.graph._graph_dep import add_atom_explicit_hydrogen_keys
from automol.graph._graph_dep import without_dummy_atoms
from automol.graph._graph_dep import without_fractional_bonds
from automol.graph._graph_dep import without_dummy_bonds
from automol.graph._graph_dep import add_bonds
from automol.graph._graph_dep import remove_atoms
from automol.graph._graph_dep import remove_bonds
from automol.graph._graph_dep import atom_bond_valences
from automol.graph._graph_dep import implicit
from automol.graph._graph_dep import explicit
from automol.graph._graph_dep import add_atoms
from automol.graph._graph_dep import backbone_keys
from automol.graph._graph_dep import atom_unsaturated_valences
from automol.graph._graph_dep import maximum_spin_multiplicity
from automol.graph._graph_dep import explicit_hydrogen_keys
from automol.graph._graph_dep import atom_neighborhood
from automol.graph._graph_dep import atom_neighborhoods
from automol.graph._graph_dep import bond_neighborhoods
from automol.graph._graph_dep import atom_sorted_neighbor_atom_keys
from automol.graph._graph_dep import atom_element_valences
from automol.graph._graph_dep import atom_neighbor_atom_key
from automol.graph._graph_dep import atoms_neighbor_atom_keys
from automol.graph._graph_dep import subgraph
from automol.graph._graph_dep import bond_induced_subgraph
from automol.graph._graph_dep import add_atom_implicit_hydrogen_valences
from automol.graph._graph_dep import atom_stereo_sorted_neighbor_atom_keys
from automol.graph._graph_dep import atom_lone_pair_counts
# stereo
from automol.graph._graph_dep import has_stereo
from automol.graph._graph_dep import atom_stereo_keys
from automol.graph._graph_dep import bond_stereo_keys
from automol.graph._graph_dep import stereo_priority_vector
from automol.graph._graph_dep import sp2_bond_keys
from automol.graph._graph_dep import dominant_resonance
from automol.graph._graph_dep import dominant_resonances
from automol.graph._graph_dep import resonances
from automol.graph._graph_dep import subresonances
from automol.graph._graph_dep import resonance_dominant_bond_orders
from automol.graph._graph_dep import atom_hybridizations
from automol.graph._graph_dep import resonance_dominant_atom_hybridizations
from automol.graph._graph_dep import atoms_sorted_neighbor_atom_keys
from automol.graph._graph_dep import dummy_atoms_neighbor_atom_key
# dep
from automol.graph._embed_dep import fake_stereo_geometry
from automol.graph._embed_dep import distance_bounds_matrices
from automol.graph._embed_dep import chirality_constraint_bounds
from automol.graph._embed_dep import planarity_constraint_bounds
from automol.graph._embed_dep import qualitative_convergence_checker_
from automol.graph._embed_dep import van_der_waals_radius
from automol.graph._embed_dep import path_distance_bounds_
from automol.graph._embed_dep import heuristic_bond_distance
from automol.graph._embed_dep import heuristic_bond_angle
from automol.graph._embed_dep import heuristic_bond_angle_distance
from automol.graph._embed_dep import heuristic_torsion_angle_distance
from automol.graph._embed_dep import shared_ring_size
from automol.graph._embed_dep import atom_shortest_paths
from automol.graph._embed_dep import closest_approach
from automol.graph._embed_dep import backbone_isomorphic
from automol.graph._embed_dep import backbone_isomorphism
from automol.graph._embed_dep import transform_keys
from automol.graph._embed_dep import union
# rings
from automol.graph._embed_dep import rings_atom_keys
from automol.graph._embed_dep import sorted_ring_atom_keys
from automol.graph._embed_dep import sorted_ring_atom_keys_from_bond_keys
from automol.graph._embed_dep import rings_bond_keys
# setters
from automol.graph._graph_base import remove_atom_stereo_parities
from automol.graph._graph_base import remove_bond_stereo_parities
# I/O
from automol.graph._graph_base import from_string
from automol.graph._graph_base import from_yaml_dictionary
# graph theory library
# setters
from automol.graph._graph import standard_keys
from automol.graph._graph import standard_keys_for_sequence
from automol.graph._graph import relabel_for_zmatrix
from automol.graph._graph import relabel_for_geometry
from automol.graph._graph import add_bonded_atom
from automol.graph._graph import add_dummy_atoms
from automol.graph._graph import insert_bonded_atom
from automol.graph._graph import insert_dummy_atoms
from automol.graph._graph import standard_keys_without_dummy_atoms
from automol.graph._graph import move_idx_to_top
# # atom properties
from automol.graph._graph import electron_count
from automol.graph._graph import atom_count
from automol.graph._graph import atom_count_by_type
from automol.graph._graph import heavy_atom_count
from automol.graph._graph import atoms_second_degree_neighbor_atom_keys
from automol.graph._graph import atoms_bond_keys
from automol.graph._graph import angle_keys
from automol.graph._graph import atom_groups
from automol.graph._res import radical_groups
from automol.graph._res import radical_group_dct
from automol.graph._res import radical_dissociation_prods
# # bond properties
from automol.graph._graph import bonds_neighbor_atom_keys
from automol.graph._graph import bonds_neighbor_bond_keys
# # other properties
from automol.graph._graph import terminal_heavy_atom_keys
from automol.graph._graph import branch
from automol.graph._graph import branch_atom_keys
from automol.graph._graph import branch_bond_keys
from automol.graph._graph import is_connected
from automol.graph._graph import connected_components
from automol.graph._graph import connected_components_atom_keys
from automol.graph._graph import shortest_path_between_atoms
from automol.graph._graph import shortest_path_between_groups
from automol.graph._graph import longest_chain
from automol.graph._graph import is_branched
from automol.graph._graph import atom_longest_chain
from automol.graph._graph import atom_longest_chains
from automol.graph._graph import union_from_sequence
# implicit/explicit hydrogen functions
# # atom properties
from automol.graph._graph import atom_explicit_hydrogen_valences
# # comparisons
from automol.graph._graph import isomorphism
from automol.graph._graph import full_isomorphism
from automol.graph._graph import full_subgraph_isomorphism
from automol.graph._graph import backbone_unique
from automol.graph._graph import unsaturated_atom_keys

# chemistry library
# # other properties
from automol.graph._graph import possible_spin_multiplicities
# # radical properies
from automol.graph._rad import isomorphic_radical_graphs
from automol.graph._rad import nonisomorphic_radical_graphs

# miscellaneous
# # bond properties
from automol.graph._graph import bond_symmetry_numbers

# resonance graph library
# # atom properties
from automol.graph._res import linear_atom_keys
from automol.graph._res import resonance_dominant_atom_centered_cumulene_keys
from automol.graph._res import resonance_dominant_bond_centered_cumulene_keys
from automol.graph._res import nonresonant_radical_atom_keys
from automol.graph._res import sigma_radical_atom_keys
from automol.graph._res import resonance_dominant_radical_atom_keys
from automol.graph._res import sing_res_dom_radical_atom_keys
# # bond properties
from automol.graph._res import one_resonance_dominant_bond_orders
from automol.graph._res import resonance_avg_bond_orders
# stereo graph library
from automol.graph._stereo import stereogenic_atom_keys
from automol.graph._stereo import stereogenic_bond_keys
from automol.graph._stereo import stereomers
from automol.graph._stereo import substereomers
from automol.graph._stereo import atoms_stereo_sorted_neighbor_atom_keys
from automol.graph._stereo import to_index_based_stereo
from automol.graph._stereo import from_index_based_stereo

# ring graph library
from automol.graph._ring import rings
from automol.graph._ring import is_ring_key_sequence
from automol.graph._ring import cycle_ring_atom_key_to_front
from automol.graph._ring import ring_arc_complement_atom_keys
from automol.graph._ring import ring_systems
from automol.graph._ring import ring_systems_atom_keys
from automol.graph._ring import ring_systems_bond_keys
from automol.graph._ring import is_ring_system
from automol.graph._ring import ring_system_decomposed_atom_keys
from automol.graph._ring import ring_systems_decomposed_atom_keys

# rotational group library
from automol.graph._rot import rotational_bond_keys
from automol.graph._rot import rotational_groups
from automol.graph._rot import rotational_symmetry_number
from automol.graph._rot import linear_segments_atom_keys

# functional group library
from automol.graph._func_group import FunctionalGroup
from automol.graph._func_group import functional_group_dct
from automol.graph._func_group import hydrocarbon_species
from automol.graph._func_group import radical_species
from automol.graph._func_group import chem_unique_atoms_of_type

# stereo geometry library
from automol.graph._stereo_geom import set_stereo_from_geometry

# embedding library (should replace graph <=> geometry library)
from automol.graph import embed

# submodules
from automol.graph import ts
from automol.graph import vmat

# constructors
import automol.create as _create
# conversions
from automol.convert.graph import inchi as _inchi
from automol.convert.graph import geometry as _geometry
from automol.convert.graph import formula as _formula

# util fxns
from automol.graph._util import ring_idxs


def from_data(atm_sym_dct, bnd_keys, atm_imp_hyd_vlc_dct=None,
              atm_ste_par_dct=None, bnd_ord_dct=None, bnd_ste_par_dct=None):
    """ construct a molecular graph from data

    format:
        gra = (atm_dct, bnd_dct)
        atm_dct := {atm_key: (atm_sym, atm_imp_hyd_vlc, atm_ste_par), ...}
        bnd_dct := {bnd_key: (bnd_ord, bnd_ste_par), ...}
        [where bnd_key := frozenset({atm1_key, atm2_key})]

    :param atm_sym_dct: atomic symbols, by atom key
    :type atm_sym_dct: dict
    :param bnd_keys: bond keys
    :type bnd_keys: set
    :param atm_imp_hyd_vlc_dct: the number of implicit hydrogens associated
        with each atom, by atom key
    :type atm_imp_hyd_vlc_dct: dict
    :param atm_ste_par_dct: atom stereo parities, by atom key
    :type atm_ste_par_dct: dict
    :param bnd_ord_dct: bond orders, by bond key
    :type bnd_ord_dct: dict
    :param bnd_ste_par_dct: bond stereo parities, by bond key
    :type bnd_ste_par_dct: dict
    """
    return _create.graph.from_data(
        atom_symbols=atm_sym_dct, bond_keys=bnd_keys,
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct
    )


def inchi(gra, stereo=False):
    """ graph => inchi
    """
    return _inchi(gra, stereo=stereo)


def stereo_inchi(gra):
    """ graph => inchi
    """
    return _inchi(gra, stereo=True)


def geometry(gra):
    """ graph => geometry
    """
    return _geometry(gra)


def formula(gra):
    """ graph => formula
    """
    return _formula(gra)


__all__ = [
    'without_bond_orders',
    'without_stereo_parities',
    'atom_explicit_hydrogen_keys',
    'add_atom_explicit_hydrogen_keys',
    'without_dummy_atoms',
    'without_fractional_bonds',
    'without_dummy_bonds',
    'add_bonds',
    'remove_atoms',
    'remove_bonds',
    'atom_bond_valences',
    'implicit',
    'explicit',
    'add_atoms',
    'backbone_keys',
    'atom_unsaturated_valences',
    'maximum_spin_multiplicity',
    'explicit_hydrogen_keys',
    'atom_neighborhood',
    'atom_neighborhoods',
    'bond_neighborhoods',
    'atom_sorted_neighbor_atom_keys',
    'atom_element_valences',
    'atom_neighbor_atom_key',
    'atoms_neighbor_atom_keys',
    'subgraph',
    'bond_induced_subgraph',
    'add_atom_implicit_hydrogen_valences',
    'atom_stereo_sorted_neighbor_atom_keys',
    'atom_lone_pair_counts',

    'atoms',
    'bonds',
    'atom_keys',
    'bond_keys',
    'atom_symbols',
    'atom_implicit_hydrogen_valences',
    'atom_stereo_parities',
    'bond_orders',
    'bond_stereo_parities',
    'set_atom_implicit_hydrogen_valences',
    'set_atom_stereo_parities',
    'set_bond_orders',
    'set_bond_stereo_parities',
    'string',
    'yaml_dictionary',
    'relabel',

    'has_stereo',
    'atom_stereo_keys',
    'bond_stereo_keys',
    'stereo_priority_vector',
    'sp2_bond_keys',
    'dominant_resonance',
    'dominant_resonances',
    'resonances',
    'subresonances',
    'resonance_dominant_bond_orders',
    'atom_hybridizations',
    'resonance_dominant_atom_hybridizations',
    'atoms_sorted_neighbor_atom_keys',
    'dummy_atoms_neighbor_atom_key',

    'fake_stereo_geometry',
    'distance_bounds_matrices',
    'chirality_constraint_bounds',
    'planarity_constraint_bounds',
    'qualitative_convergence_checker_',
    'van_der_waals_radius',
    'path_distance_bounds_',
    'heuristic_bond_distance',
    'heuristic_bond_angle',
    'heuristic_bond_angle_distance',
    'heuristic_torsion_angle_distance',
    'shared_ring_size',
    'atom_shortest_paths',
    'closest_approach',
    'backbone_isomorphic',
    'backbone_isomorphism',
    'transform_keys',
    'union',

    'rings_atom_keys',
    'sorted_ring_atom_keys',
    'sorted_ring_atom_keys_from_bond_keys',
    # constructors
    'from_data',
    # getters
    'atoms',
    'bonds',
    'atom_keys',
    'bond_keys',
    'atom_symbols',
    'atom_implicit_hydrogen_valences',
    'atom_stereo_parities',
    'bond_orders',
    'bond_stereo_parities',
    # setters
    'remove_atom_stereo_parities',
    'remove_bond_stereo_parities',
    # I/O
    'from_string',
    'from_yaml_dictionary',

    # setters
    'relabel',
    'standard_keys',
    'standard_keys_for_sequence',
    'relabel_for_zmatrix',
    'relabel_for_geometry',
    'add_bonded_atom',
    'add_dummy_atoms',
    'insert_bonded_atom',
    'insert_dummy_atoms',
    'standard_keys_without_dummy_atoms',
    'move_idx_to_top',

    # graph theory library
    # # atom properties
    'electron_count',
    'atom_count',
    'atom_count_by_type',
    'heavy_atom_count',
    'atoms_second_degree_neighbor_atom_keys',
    'atoms_bond_keys',
    'angle_keys',
    'atom_groups',
    'radical_groups',
    'radical_group_dct',
    'radical_dissociation_prods',
    # # bond properties
    'bonds_neighbor_atom_keys',
    'bonds_neighbor_bond_keys',
    # # other properties
    'terminal_heavy_atom_keys',
    'branch',
    'branch_atom_keys',
    'branch_bond_keys',
    'rings',
    'rings_bond_keys',
    'is_connected',
    'connected_components',
    'connected_components_atom_keys',
    'shortest_path_between_atoms',
    'shortest_path_between_groups',
    'longest_chain',
    'is_branched',
    'atom_longest_chain',
    'atom_longest_chains',
    'union_from_sequence',
    # implicit/explicit hydrogen functions
    # # atom properties
    'atom_explicit_hydrogen_valences',
    'atom_explicit_hydrogen_keys',
    # # comparisons
    'full_subgraph_isomorphism',
    'isomorphism',
    'full_isomorphism',
    'full_subgraph_isomorphism',
    'backbone_unique',

    # chemistry library
    # # atom properties
    'unsaturated_atom_keys',
    # # other properties
    'possible_spin_multiplicities',
    'bond_symmetry_numbers',
    # # radical properies
    'isomorphic_radical_graphs',
    'nonisomorphic_radical_graphs',

    # resonance library
    # # atom properties
    'linear_atom_keys',
    'resonance_dominant_atom_centered_cumulene_keys',
    'resonance_dominant_bond_centered_cumulene_keys',
    'nonresonant_radical_atom_keys',
    'sigma_radical_atom_keys',
    'resonance_dominant_radical_atom_keys',
    'sing_res_dom_radical_atom_keys',
    # # bond properties
    'one_resonance_dominant_bond_orders',
    'resonance_avg_bond_orders',

    # stereo graph library
    'stereogenic_atom_keys',
    'stereogenic_bond_keys',
    'stereomers',
    'substereomers',
    'atoms_stereo_sorted_neighbor_atom_keys',
    'to_index_based_stereo',
    'from_index_based_stereo',

    # ring graph library
    'rings',
    'is_ring_key_sequence',
    'cycle_ring_atom_key_to_front',
    'ring_arc_complement_atom_keys',
    'ring_systems',
    'ring_systems_atom_keys',
    'ring_systems_bond_keys',
    'is_ring_system',
    'ring_system_decomposed_atom_keys',
    'ring_systems_decomposed_atom_keys',

    # rotational group library
    'rotational_bond_keys',
    'rotational_groups',
    'rotational_symmetry_number',
    'linear_segments_atom_keys',

    # functional group library
    'FunctionalGroup',
    'functional_group_dct',
    'hydrocarbon_species',
    'radical_species',
    'chem_unique_atoms_of_type',

    # stereo geometry library
    'set_stereo_from_geometry',

    # conversions,
    'inchi',
    'formula',

    # embedding library
    'embed',

    # submodules
    'ts',
    'vmat',

    # util fxns
    'ring_idxs'
]
