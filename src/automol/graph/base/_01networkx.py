""" networkx interface

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import operator

import networkx

from ... import util
from ._00core import (
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_parities,
    atom_symbols,
    bond_keys,
    bond_orders,
    bond_stereo_parities,
)


def from_graph(gra, node_attrib_dct=None, edge_attrib_dct=None):
    """networkx graph object from a molecular graph"""
    nxg = networkx.Graph()
    nxg.add_nodes_from(atom_keys(gra))
    nxg.add_edges_from(bond_keys(gra))
    networkx.set_node_attributes(nxg, atom_symbols(gra), "symbol")
    networkx.set_node_attributes(
        nxg, atom_implicit_hydrogens(gra), "implicit_hydrogens"
    )
    networkx.set_node_attributes(nxg, atom_stereo_parities(gra), "stereo_parity")
    networkx.set_edge_attributes(nxg, bond_orders(gra), "order")
    networkx.set_edge_attributes(nxg, bond_stereo_parities(gra), "stereo_parity")

    if node_attrib_dct is not None:
        for name, dct in node_attrib_dct.items():
            networkx.set_node_attributes(nxg, dct, name)

    if edge_attrib_dct is not None:
        for name, dct in edge_attrib_dct.items():
            networkx.set_edge_attributes(nxg, dct, name)

    return nxg


def minimum_cycle_basis(nxg: networkx.Graph, weight: str | None = None):
    """Cycle basis for the graph."""

    def _order_and_normalize(rkeys):
        """Order and normalize the ring keys."""
        cyc_iter = list(networkx.simple_cycles(nxg.subgraph(rkeys)))
        rkeys = next((ks for ks in cyc_iter if set(ks) == set(rkeys)), None)
        assert rkeys is not None, f"Failed to identify cycle for {rkeys} in\n{nxg}"
        return util.ring.normalize(rkeys)

    rkeys_lst = networkx.algorithms.cycles.minimum_cycle_basis(nxg, weight=weight)
    # Ensure that the ordering is correct (not guaranteed by minimum cycle basis)
    rkeys_lst = tuple(sorted(map(_order_and_normalize, rkeys_lst)))
    return rkeys_lst


def connected_component_atom_keys(nxg):
    """atom keys for the connected components in this graph"""
    return tuple(map(frozenset, networkx.algorithms.connected_components(nxg)))


def all_pairs_shortest_path(nxg):
    """shortest path between any two vertices in the graph"""
    return networkx.all_pairs_shortest_path(nxg)


def simple_paths(nxg, key1, key2, overlap: bool = False):
    """Find simple paths between any two vertices in the graph.

    :param overlap: Whether to include overlapping paths. If not, the shortest one is
        retained.
    """

    def _are_overlapping(path1, path2):
        return set(path1[1:-2]) & set(path2[1:-2])

    paths = tuple(map(tuple, networkx.all_simple_paths(nxg, source=key1, target=key2)))
    if not overlap:
        paths = tuple(
            min(ps, key=len)
            for ps in util.equivalence_partition(paths, _are_overlapping)
        )

    return paths


def all_isomorphisms(nxg1, nxg2):
    """Find all possible isomorphisms between two graphs"""

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq
    )

    return tuple(matcher.isomorphisms_iter())


def isomorphism(nxg1, nxg2):
    """graph isomorphism"""

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq
    )

    iso_dct = None
    if matcher.is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct


def subgraph_isomorphism(nxg1, nxg2):
    """subgraph isomorphism -- a subgraph of G1 is isomorphic to G2"""

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq
    )

    iso_dct = None
    if matcher.subgraph_is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct


def weighted_maximal_matching(nxg, edge_attrib_name):
    """subgraph isomorphism -- a subgraph of G1 is isomorphic to G2"""
    ret = networkx.max_weight_matching(nxg, weight=edge_attrib_name)
    bnd_keys = frozenset(map(frozenset, ret))
    return bnd_keys
