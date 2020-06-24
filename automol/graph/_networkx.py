""" networkx interface
"""
import operator
import networkx
from automol.graph._graph_base import atom_keys
from automol.graph._graph_base import bond_keys
from automol.graph._graph_base import atom_symbols
from automol.graph._graph_base import atom_implicit_hydrogen_valences
from automol.graph._graph_base import atom_stereo_parities
from automol.graph._graph_base import bond_orders
from automol.graph._graph_base import bond_stereo_parities


def from_graph(xgr):
    """ networkx graph object from a molecular graph
    """
    nxg = networkx.Graph()
    nxg.add_nodes_from(atom_keys(xgr))
    nxg.add_edges_from(bond_keys(xgr))
    networkx.set_node_attributes(nxg, atom_symbols(xgr), 'symbol')
    networkx.set_node_attributes(nxg, atom_implicit_hydrogen_valences(xgr),
                                 'implicit_hydrogen_valence')
    networkx.set_node_attributes(nxg, atom_stereo_parities(xgr),
                                 'stereo_parity')
    networkx.set_edge_attributes(nxg, bond_orders(xgr),
                                 'order')
    networkx.set_edge_attributes(nxg, bond_stereo_parities(xgr),
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


def isomorphism(nxg1, nxg2):
    """ graph isomorphism
    """

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=operator.eq, edge_match=operator.eq)

    iso_dct = None
    if matcher.is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct


if __name__ == '__main__':
    GRA = (
        {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
         6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
         frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
         frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
         frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)})

    GRA2 = (
        {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
         6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
         frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
         frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
         frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)})

    NXG = from_graph(GRA)
    NXG2 = from_graph(GRA2)
    print(isomorphism(NXG, NXG2))
    networkx.write_yaml(NXG, 'test.g.yaml')
