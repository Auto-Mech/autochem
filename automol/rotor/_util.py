"""
 Various utility functions
"""

import automol


def graph_with_keys(zma):
    """ Generate the graph
    """

    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    gra = automol.geom.graph(geo)
    lin_keys = sorted(gdummy_key_dct.keys())

    return gra, lin_keys


def sort_tors_names(tors_names):
    """ sort torsional names so that Dn where n is ascending order
    """
    tors_names = list(tors_names)
    tors_names.sort(key=lambda x: int(x.split('D')[1]))
    return tors_names
