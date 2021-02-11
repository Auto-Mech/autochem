"""
 Various utility functions
"""

import automol


def graph_with_keys(zma, zrxn=None):
    """ Generate the graph
    """

    if zrxn is None:
        # geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
        # gra = automol.geom.graph(geo)
        # lin_keys = sorted(gdummy_key_dct.keys())
        gra = automol.zmat.graph(zma, stereo=True, dummy=True)
        lin_keys = sorted(
            automol.graph.dummy_atoms_neighbor_atom_key(gra).values())
        print('gra\n', automol.graph.string(gra))
        print('linkeys', lin_keys)
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
