""" graph conversion helpers
"""
from itertools import combinations as _combinations
import numpy
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates
from .._cnst.graph import from_data as _graph_from_data

RQ_BOND_MAX = 3.5
RH_BOND_MAX = 2.5


def connectivity_graph(geo):
    """ connectivity graph
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo)

    def _are_bonded(idx_pair):
        xyz1, xyz2 = map(xyzs.__getitem__, idx_pair)
        sym1, sym2 = map(syms.__getitem__, idx_pair)
        dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
        return (False if 'X' in (sym1, sym2) else
                (dist < RH_BOND_MAX) if 'H' in (sym1, sym2) else
                (dist < RQ_BOND_MAX))

    idxs = range(len(xyzs))
    bnds = list(filter(_are_bonded, _combinations(idxs, r=2)))
    return _graph_from_data(atom_symbols=syms, bond_keys=bnds)
