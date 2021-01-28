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
