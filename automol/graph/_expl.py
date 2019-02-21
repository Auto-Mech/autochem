""" backbone and explicit hydrogen library
"""
from itertools import chain as _chain
from future.moves.itertools import filterfalse as _filterfalse
import numpy
from ._networkx import from_graph as _nxg_from_graph
from ._networkx import isomorphism as _nxg_isomorphism
from ._math import unique as _unique
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict import keys_by_value as _keys_by_value
from ._dict import transform_values as _transform_values
from ._dict import transform_items_to_values as _transform_items_to_values
from ._core import add_atoms as _add_atoms
from ._core import add_bonds as _add_bonds
from ._core import atom_keys as _atom_keys
from ._core import atom_symbols as _atom_symbols
from ._core import (atom_implicit_hydrogen_valences as
                    _atom_implicit_hydrogen_valences)
from ._core import (set_atom_implicit_hydrogen_valences as
                    _set_atom_implicit_hydrogen_valences)
from ._graph import atom_neighbor_keys as _atom_neighbor_keys
from ._graph import atom_neighborhoods as _atom_neighborhoods
from ._graph import delete_atoms as _delete_atoms


# atom properties
def atom_explicit_hydrogen_valences(xgr):
    """ explicit hydrogen valences, by atom
    """
    return _transform_values(atom_explicit_hydrogen_keys(xgr), len)


def atom_explicit_hydrogen_keys(xgr):
    """ explicit hydrogen valences, by atom
    """
    def _explicit_hydrogen_keys(atm_key, atm_nbh):
        return frozenset(explicit_hydrogen_keys(atm_nbh) - {atm_key})

    atm_exp_hyd_keys_dct = _transform_items_to_values(_atom_neighborhoods(xgr),
                                                      _explicit_hydrogen_keys)
    return atm_exp_hyd_keys_dct


# other properties
def backbone_keys(xgr):
    """ backbone atom keys
    """
    bbn_keys = _atom_keys(xgr) - explicit_hydrogen_keys(xgr)
    return bbn_keys


def explicit_hydrogen_keys(xgr):
    """ explicit hydrogen keys (H types: explicit, implicit, backbone)
    """
    hyd_keys = _keys_by_value(_atom_symbols(xgr), lambda x: x == 'H')
    atm_ngb_keys_dct = _atom_neighbor_keys(xgr)

    def _is_backbone(hyd_key):
        return all(ngb_key in hyd_keys and hyd_key < ngb_key
                   for ngb_key in atm_ngb_keys_dct[hyd_key])

    exp_hyd_keys = frozenset(_filterfalse(_is_backbone, hyd_keys))
    return exp_hyd_keys


# transformations
def add_explicit_hydrogens(xgr, atm_exp_hyd_vlc_dct):
    """ add explicit hydrogens by atom
    """
    assert set(atm_exp_hyd_vlc_dct.keys()) <= _atom_keys(xgr)
    for atm_key, atm_exp_hyd_vlc in atm_exp_hyd_vlc_dct.items():
        next_atm_key = max(_atom_keys(xgr)) + 1
        atm_exp_hyd_keys = set(range(next_atm_key,
                                     next_atm_key + atm_exp_hyd_vlc))
        atm_exp_hyd_bnd_keys = {frozenset({atm_key, atm_exp_hyd_key})
                                for atm_exp_hyd_key in atm_exp_hyd_keys}
        atm_exp_hyd_sym_dct = _by_key({}, atm_exp_hyd_keys, fill_val='H')
        xgr = _add_atoms(xgr, atm_exp_hyd_sym_dct)
        xgr = _add_bonds(xgr, atm_exp_hyd_bnd_keys)
    return xgr


def implicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms implicit
    """
    atm_keys = backbone_keys(xgr) if atm_keys is None else atm_keys
    atm_keys = list(atm_keys)
    atm_imp_hyd_vlcs = _values_by_key(
        _atom_implicit_hydrogen_valences(xgr), atm_keys)

    atm_exp_hyd_keys = _values_by_key(
        atom_explicit_hydrogen_keys(xgr), atm_keys)
    atm_exp_hyd_vlcs = tuple(map(len, atm_exp_hyd_keys))
    atm_tot_hyd_vlcs = numpy.add(atm_imp_hyd_vlcs, atm_exp_hyd_vlcs)

    exp_hyd_keys = tuple(_chain(*atm_exp_hyd_keys))

    xgr = _set_atom_implicit_hydrogen_valences(
        xgr, dict(zip(atm_keys, atm_tot_hyd_vlcs)))
    xgr = _delete_atoms(xgr, exp_hyd_keys)
    return xgr


def explicit(xgr, atm_keys=None):
    """ make the hydrogens at these atoms explicit
    """
    atm_keys = backbone_keys(xgr) if atm_keys is None else atm_keys
    atm_keys = list(atm_keys)
    atm_imp_hyd_vlcs = _values_by_key(
        _atom_implicit_hydrogen_valences(xgr), atm_keys)

    xgr = _set_atom_implicit_hydrogen_valences(
        xgr, _by_key({}, atm_keys, fill_val=0))
    xgr = add_explicit_hydrogens(
        xgr, dict(zip(atm_keys, atm_imp_hyd_vlcs)))
    return xgr


# comparisons
def backbone_isomorphic(xgr1, xgr2):
    """ are these molecular graphs backbone isomorphic?
    """
    return backbone_isomorphism(xgr1, xgr2) is not None


def backbone_isomorphism(xgr1, xgr2):
    """ graph backbone isomorphism

    for implicit graphs, this is the relabeling of `xgr1` to produce `xgr2`
    for other graphs, it gives the correspondences between backbone atoms
    """
    xgr1 = implicit(xgr1)
    xgr2 = implicit(xgr2)
    nxg1 = _nxg_from_graph(xgr1)
    nxg2 = _nxg_from_graph(xgr2)
    iso_dct = _nxg_isomorphism(nxg1, nxg2)
    return iso_dct


def backbone_unique(xgrs):
    """ unique non-isomorphic graphs from a series
    """
    xgrs = _unique(xgrs, equiv=backbone_isomorphic)
    return xgrs
