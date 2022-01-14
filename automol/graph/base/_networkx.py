""" networkx interface

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import operator
import networkx
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_orders
from automol.graph.base._core import bond_stereo_parities


def from_graph(gra):
    """ networkx graph object from a molecular graph
    """
    nxg = networkx.Graph()
    nxg.add_nodes_from(atom_keys(gra))
    nxg.add_edges_from(bond_keys(gra))
    networkx.set_node_attributes(nxg, atom_symbols(gra), 'symbol')
    networkx.set_node_attributes(nxg, atom_implicit_hydrogen_valences(gra),
                                 'implicit_hydrogen_valence')
    networkx.set_node_attributes(nxg, atom_stereo_parities(gra),
                                 'stereo_parity')
    networkx.set_edge_attributes(nxg, bond_orders(gra),
                                 'order')
    networkx.set_edge_attributes(nxg, bond_stereo_parities(gra),
                                 'stereo_parity')

    return nxg


def minimum_cycle_basis(nxg):
    """ minimum cycle basis for the graph
    """
    rng_atm_keys_lst = networkx.algorithms.cycles.minimum_cycle_basis(nxg)
    return frozenset(map(frozenset, rng_atm_keys_lst))


def connected_component_atom_keys(nxg):
    """ atom keys for the connected components in this graph
    """
    return tuple(map(frozenset, networkx.algorithms.connected_components(nxg)))


def all_pairs_shortest_path(nxg):
    """ shortest path between any two vertices in the graph
    """
    return networkx.all_pairs_shortest_path(nxg)


def all_isomorphisms(nxg1, nxg2):
    """ Find all possible isomorphisms between two graphs
    """

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq)

    return tuple(matcher.isomorphisms_iter())


def isomorphism(nxg1, nxg2):
    """ graph isomorphism
    """

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq)

    iso_dct = None
    if matcher.is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct


def subgraph_isomorphism(nxg1, nxg2):
    """ subgraph isomorphism -- a subgraph of G1 is isomorphic to G2
    """

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq)

    iso_dct = None
    if matcher.subgraph_is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct
