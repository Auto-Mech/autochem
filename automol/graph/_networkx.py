""" networkx interface
"""
import networkx


def from_graph(xgr):
    """ networkx graph object from a molecular graph
    """
    atms, bnds = xgr
    nxg = networkx.Graph()
    for atm_key, atm_vals in atms.items():
        nxg.add_node(atm_key, props=atm_vals)
    for bnd_key, bnd_val in bnds.items():
        nxg.add_edge(*bnd_key, props=bnd_val)
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

    def _same_props(dct1, dct2):
        return dct1['props'] == dct2['props']

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=_same_props, edge_match=_same_props)

    iso_dct = None
    if matcher.is_isomorphic():
        iso_dct = dict(matcher.mapping)

    return iso_dct
