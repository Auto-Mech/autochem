""" TS z-matrices for specific reaction classes
"""
from ..graph import ts
from ..zmat import distance_coordinate_name
from ._0core import Reaction, ts_graph


# Z-Matrix coordinate functions
def zmatrix_coordinate_names(zrxn: Reaction, zma):
    """ Get the Z-matrix coordinate names for the forming and
        breaking bonds of a reaction

        It is not always guaranteed that the bond keys will be
        present in the transition state Z-Matrix. For these cases,
        the name will be returned as None.

        :param zrxn: a Reaction object
        :rtype: str
    """

    def _zma_names(zma, bnd_keys):
        _names = ()
        for keys in bnd_keys:
            try:
                name = distance_coordinate_name(zma, *keys)
            except AssertionError:
                name = None
            _names += (name,)

        return _names

    frm_bnd_keys = ts.forming_bond_keys(ts_graph(zrxn))
    brk_bnd_keys = ts.breaking_bond_keys(ts_graph(zrxn))

    return (_zma_names(zma, frm_bnd_keys), _zma_names(zma, brk_bnd_keys))
