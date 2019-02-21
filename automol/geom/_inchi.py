""" functions operating on cartesian geometries (angstroms)
"""
from ._core import coordinates as _coordinates
from ._graph import connectivity_graph as _connectivity_graph
from ..graph import inchi as _inchi_from_connectivity_graph
from ..graph import (stereo_inchi_from_coordinates as
                     _stereo_inchi_from_coordinates)


def inchi(geo):
    """ connectivity InChI string of a cartesian geometry
    """
    cgr = _connectivity_graph(geo)
    return _inchi_from_connectivity_graph(cgr)


def stereo_inchi(geo):
    """ stereo-specific InChI string of a cartesian geometry
    """
    cgr = _connectivity_graph(geo)
    atm_xyz_dct = dict(enumerate(_coordinates(geo)))
    ich = _stereo_inchi_from_coordinates(cgr, atm_xyz_dct)
    return ich
