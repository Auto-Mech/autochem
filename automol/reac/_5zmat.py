""" TS z-matrices for specific reaction classes
"""
import automol.geom
import automol.graph
from automol.graph import ts
from automol.zmat import distance_coordinate_name
from automol.reac._0core import Reaction
from automol.reac._0core import ts_graph


def ts_zmatrix(rxn: Reaction, ts_geo):
    """ reaction-class-specific embedding info

    :param rxn: a hydrogen migration Reaction object
    :param ts_geo: the TS geometry
    :returns: the TS z-matrix, the row keys, and the dummy index dictionary
    """
    tsg = ts_graph(rxn)
    ts_zma, dc_ = automol.geom.zmatrix_with_conversion_info(ts_geo, gra=tsg)
    return ts_zma, dc_


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

    frm_bnd_keys = ts.ts_forming_bond_keys(ts_graph(zrxn))
    brk_bnd_keys = ts.ts_breaking_bond_keys(ts_graph(zrxn))

    return (_zma_names(zma, frm_bnd_keys), _zma_names(zma, brk_bnd_keys))
