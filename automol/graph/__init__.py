""" molecular graph
"""
# getters
from automol.graph._graph_base import atoms
from automol.graph._graph_base import bonds
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_implicit_hydrogen_valences
from automol.graph._graph_base import atom_stereo_parities
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import bond_stereo_parities
# setters
from automol.graph._graph_base import set_atom_implicit_hydrogen_valences
from automol.graph._graph_base import set_atom_stereo_parities
from automol.graph._graph_base import set_bond_orders
from automol.graph._graph_base import set_bond_stereo_parities
# I/O
from automol.graph._graph_base import string
from automol.graph._graph_base import from_string

# graph theory library
# setters
from automol.graph._graph import relabel
from automol.graph._graph import standard_keys
from automol.graph._graph import standard_keys_for_sequence
from automol.graph._graph import transform_keys
from automol.graph._graph import relabel_for_zmatrix
from automol.graph._graph import relabel_for_geometry
from automol.graph._graph import add_atom_implicit_hydrogen_valences
from automol.graph._graph import without_bond_orders
from automol.graph._graph import without_stereo_parities
from automol.graph._graph import add_atoms
from automol.graph._graph import add_bonds
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
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import atom_neighbor_key
from automol.graph._graph import sorted_atom_neighbor_keys
from automol.graph._graph import atom_second_degree_neighbor_keys
from automol.graph._graph import dummy_atom_neighbor_keys
from automol.graph._graph import atom_bond_keys
from automol.graph._graph import angle_keys
from automol.graph._graph import atom_neighborhoods
from automol.graph._graph import atom_neighborhood
from automol.graph._graph import atom_groups
from automol.graph._res import radical_groups
from automol.graph._res import radical_group_dct
from automol.graph._res import radical_dissociation_prods
# # bond properties
from automol.graph._graph import bond_neighbor_keys
from automol.graph._graph import bond_neighbor_bonds
from automol.graph._graph import bond_neighborhoods
# # other properties
from automol.graph._graph import terminal_heavy_atom_keys
from automol.graph._graph import branch
from automol.graph._graph import branch_atom_keys
from automol.graph._graph import branch_bond_keys
from automol.graph._graph import is_connected
from automol.graph._graph import connected_components
from automol.graph._graph import connected_components_atom_keys
from automol.graph._graph import atom_shortest_paths
from automol.graph._graph import shortest_path_between_groups
from automol.graph._graph import longest_chain
from automol.graph._graph import atom_longest_chain
from automol.graph._graph import atom_longest_chains
from automol.graph._graph import union
from automol.graph._graph import union_from_sequence
from automol.graph._graph import subgraph
from automol.graph._graph import bond_induced_subgraph
# # transformations
from automol.graph._graph import remove_atoms
from automol.graph._graph import remove_bonds
from automol.graph._graph import without_dummy_atoms
from automol.graph._graph import add_ts_bonds
from automol.graph._graph import ts_bond_keys
from automol.graph._graph import without_ts_bonds
# implicit/explicit hydrogen functions
# # atom properties
from automol.graph._graph import atom_explicit_hydrogen_valences
from automol.graph._graph import atom_explicit_hydrogen_keys
# # other properties
from automol.graph._graph import backbone_keys
from automol.graph._graph import explicit_hydrogen_keys
# # transformations
from automol.graph._graph import add_atom_explicit_hydrogen_keys
from automol.graph._graph import implicit
from automol.graph._graph import explicit
# # comparisons
from automol.graph._graph import full_isomorphism
from automol.graph._graph import full_subgraph_isomorphism
from automol.graph._graph import backbone_isomorphic
from automol.graph._graph import backbone_isomorphism
from automol.graph._graph import backbone_unique

# chemistry library
# # atom properties
from automol.graph._graph import atom_element_valences
from automol.graph._graph import atom_lone_pair_counts
from automol.graph._graph import atom_bond_valences
from automol.graph._graph import atom_unsaturated_valences
from automol.graph._graph import unsaturated_atom_keys
# # other properties
from automol.graph._graph import maximum_spin_multiplicity
from automol.graph._graph import possible_spin_multiplicities
# # radical properies
from automol.graph._rad import isomorphic_radical_graphs

# miscellaneous
# # bond properties
from automol.graph._graph import bond_symmetry_numbers

# resonance graph library
# # atom properties
from automol.graph._res import atom_hybridizations
from automol.graph._res import resonance_dominant_atom_hybridizations
from automol.graph._res import linear_atom_keys
from automol.graph._res import resonance_dominant_atom_centered_cumulene_keys
from automol.graph._res import resonance_dominant_bond_centered_cumulene_keys
from automol.graph._res import nonresonant_radical_atom_keys
from automol.graph._res import resonance_dominant_radical_atom_keys
from automol.graph._res import sing_res_dom_radical_atom_keys
# # bond properties
from automol.graph._res import resonance_dominant_bond_orders
from automol.graph._res import one_resonance_dominant_bond_orders
from automol.graph._res import resonance_avg_bond_orders
# # transformations
from automol.graph._res import resonances
from automol.graph._res import subresonances
from automol.graph._res import dominant_resonances
from automol.graph._res import dominant_resonance

# stereo graph library
from automol.graph._stereo import has_stereo
from automol.graph._stereo import atom_stereo_keys
from automol.graph._stereo import bond_stereo_keys
from automol.graph._stereo import stereo_priority_vector
from automol.graph._stereo import stereogenic_atom_keys
from automol.graph._stereo import stereogenic_bond_keys
from automol.graph._stereo import sp2_bond_keys
from automol.graph._stereo import stereomers
from automol.graph._stereo import substereomers
from automol.graph._stereo import stereo_sorted_atom_neighbor_keys

# ring graph library
from automol.graph._ring import rings
from automol.graph._ring import rings_atom_keys
from automol.graph._ring import rings_bond_keys
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
from automol.graph._func_group import Fgroup
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
import automol.create.graph
# conversions
import automol.convert.graph

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
    return automol.create.graph.from_data(
        atom_symbols=atm_sym_dct, bond_keys=bnd_keys,
        atom_implicit_hydrogen_valences=atm_imp_hyd_vlc_dct,
        atom_stereo_parities=atm_ste_par_dct, bond_orders=bnd_ord_dct,
        bond_stereo_parities=bnd_ste_par_dct
    )


def inchi(gra, remove_stereo=True):
    """ graph => inchi
    """
    return automol.convert.graph.inchi(gra, remove_stereo=remove_stereo)


def geometry(gra):
    """ graph => geometry
    """
    return automol.convert.graph.geometry(gra)


def formula(gra):
    """ graph => formula
    """
    return automol.convert.graph.formula(gra)


__all__ = [
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
    'set_atom_implicit_hydrogen_valences',
    'set_atom_stereo_parities',
    'set_bond_orders',
    'set_bond_stereo_parities',
    # I/O
    'string',
    'from_string',

    # setters
    'relabel',
    'standard_keys',
    'standard_keys_for_sequence',
    'transform_keys',
    'relabel_for_zmatrix',
    'relabel_for_geometry',
    'add_atom_implicit_hydrogen_valences',
    'without_bond_orders',
    'without_stereo_parities',
    'add_atoms',
    'add_bonds',
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
    'atom_neighbor_keys',
    'atom_neighbor_key',
    'sorted_atom_neighbor_keys',
    'atom_second_degree_neighbor_keys',
    'dummy_atom_neighbor_keys',
    'atom_bond_keys',
    'angle_keys',
    'atom_neighborhoods',
    'atom_neighborhood',
    'atom_groups',
    'radical_groups',
    'radical_group_dct',
    'radical_dissociation_prods',
    # # bond properties
    'bond_neighbor_keys',
    'bond_neighbor_bonds',
    'bond_neighborhoods',
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
    'atom_shortest_paths',
    'shortest_path_between_groups',
    'longest_chain',
    'atom_longest_chain',
    'atom_longest_chains',
    'union',
    'union_from_sequence',
    'subgraph',
    'bond_induced_subgraph',
    # # transformations
    'remove_atoms',
    'remove_bonds',
    'without_dummy_atoms',
    'add_ts_bonds',
    'ts_bond_keys',
    'without_ts_bonds',

    # implicit/explicit hydrogen functions
    # # atom properties
    'atom_explicit_hydrogen_valences',
    'atom_explicit_hydrogen_keys',
    # # other properties
    'backbone_keys',
    'explicit_hydrogen_keys',
    # # transformations
    'add_atom_explicit_hydrogen_keys',
    'implicit',
    'explicit',
    # # comparisons
    'full_subgraph_isomorphism',
    'full_isomorphism',
    'full_subgraph_isomorphism',
    'backbone_isomorphic',
    'backbone_isomorphism',
    'backbone_unique',

    # chemistry library
    # # atom properties
    'atom_element_valences',
    'atom_lone_pair_counts',
    'atom_bond_valences',
    'atom_unsaturated_valences',
    'unsaturated_atom_keys',
    # # other properties
    'maximum_spin_multiplicity',
    'possible_spin_multiplicities',
    'bond_symmetry_numbers',
    # # radical properies
    'isomorphic_radical_graphs',

    # resonance library
    # # atom properties
    'atom_hybridizations',
    'resonance_dominant_atom_hybridizations',
    'linear_atom_keys',
    'resonance_dominant_atom_centered_cumulene_keys',
    'resonance_dominant_bond_centered_cumulene_keys',
    'nonresonant_radical_atom_keys',
    'resonance_dominant_radical_atom_keys',
    'sing_res_dom_radical_atom_keys',
    # # bond properties
    'resonance_dominant_bond_orders',
    'one_resonance_dominant_bond_orders',
    'resonance_avg_bond_orders',
    # # transformations
    'resonances',
    'subresonances',
    'dominant_resonances',
    'dominant_resonance',

    # stereo graph library
    'has_stereo',
    'atom_stereo_keys',
    'bond_stereo_keys',
    'stereo_priority_vector',
    'stereogenic_atom_keys',
    'stereogenic_bond_keys',
    'sp2_bond_keys',
    'stereomers',
    'substereomers',
    'stereo_sorted_atom_neighbor_keys',

    # ring graph library
    'rings',
    'rings_atom_keys',
    'rings_bond_keys',
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
    'Fgroup',
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
