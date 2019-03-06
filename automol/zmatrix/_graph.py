""" z-matrix to graph conversion
"""
from ..geom import connectivity_graph as _geom_connectivity_graph
from ._geom import geometry as _geometry


def connectivity_graph(zma):
    """ connectivity graph from a z-matrix
    """
    gra = _geom_connectivity_graph(_geometry(zma))
    return gra
