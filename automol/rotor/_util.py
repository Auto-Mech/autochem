"""
 Various utility functions
"""

import automol


def graph_with_keys(zma, zrxn=None):
    """ Generate the graph
    """

    if zrxn is None:
        gra = automol.zmat.graph(zma, stereo=True, dummy=True)
        lin_keys = sorted(
            automol.graph.dummy_atoms_neighbor_atom_key(gra).values())
    else:
        gra = zrxn.forward_ts_graph
        lin_keys = sorted(
            automol.graph.dummy_atoms_neighbor_atom_key(gra).values())

    return gra, lin_keys


def sort_tors_names(tors_names):
    """ sort torsional names so that Dn where n is ascending order
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
