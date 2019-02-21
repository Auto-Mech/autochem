""" specific resonance graph functions
"""
from functools import partial as _partial
from itertools import product as _product
import numpy
from ..atom import valence as _atom_valence
from ._dict import values_by_key as _values_by_key
from ._dict import transform_values as _transform_values
from ._core import frozen as _frozen
from ._core import atom_keys as _atom_keys
from ._core import bond_keys as _bond_keys
from ._core import atom_symbols as _atom_symbols
from ._core import bond_orders as _bond_orders
from ._core import set_bond_orders as _set_bond_orders
from ._core import (atom_implicit_hydrogen_valences as
                    _atom_implicit_hydrogen_valences)
from ._core import without_bond_orders as _without_bond_orders
from ._graph import atom_bond_keys as _atom_bond_keys
from ._graph import atom_neighborhoods as _atom_neighborhoods


# atom properties
def atom_bond_valences(rgr):
    """ bond valences, by atom
    """
    atm_keys = list(_atom_keys(rgr))
    atm_nbhs = _values_by_key(_atom_neighborhoods(rgr), atm_keys)
    atm_exp_bnd_vlcs = [sum(_bond_orders(nbh).values()) for nbh in atm_nbhs]
    atm_imp_hyd_vlcs = _values_by_key(
        _atom_implicit_hydrogen_valences(rgr), atm_keys)
    atm_bnd_vlcs = numpy.add(atm_exp_bnd_vlcs, atm_imp_hyd_vlcs)
    atm_bnd_vlc_dct = dict(zip(atm_keys, atm_bnd_vlcs))
    return atm_bnd_vlc_dct


def atom_radical_valences(rgr):
    """ radical valences, by atom
    """
    atm_keys = list(_atom_keys(rgr))
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_tot_vlcs = _values_by_key(_atom_total_valences(rgr), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    return dict(zip(atm_keys, atm_rad_vlcs))


def _atom_total_valences(xgr):
    atm_sym_dct = _atom_symbols(xgr)
    atm_tot_vlc_dct = _transform_values(atm_sym_dct, func=_atom_valence)
    return atm_tot_vlc_dct


# bond properties
def resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(_bond_keys(rgr))
    bnd_ords_by_res = [_values_by_key(_bond_orders(dom_rgr), bnd_keys)
                       for dom_rgr in dominant_resonances(rgr)]
    bnd_ords_lst = list(map(frozenset, zip(*bnd_ords_by_res)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


# other properties
def maximum_spin_multiplicity(rgr):
    """ the highest possible spin multiplicity for this molecular graph
    """
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), _atom_keys(rgr))
    return sum(atm_rad_vlcs) + 1


def possible_spin_multiplicities(rgr):
    """ possible spin multiplicities for this molecular graph
    """
    mult_max = maximum_spin_multiplicity(rgr)
    mult_min = 2 if mult_max % 2 == 0 else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


# transformations
def dominant_resonance(rgr):
    """ *a* dominant (minimum spin/maximum pi) resonance graph
    """
    return next(iter(dominant_resonances(rgr)))


def dominant_resonances(rgr):
    """ all dominant (minimum spin/maximum pi) resonance graphs
    """
    rgrs = resonances(rgr)
    mult_min = min(map(maximum_spin_multiplicity, rgrs))
    dom_rgrs = tuple(
        rgr for rgr in rgrs if maximum_spin_multiplicity(rgr) == mult_min)
    return dom_rgrs


def resonances(rgr):
    """ all resonance graphs with this connectivity
    """
    return subresonances(_without_bond_orders(rgr))


def subresonances(rgr):
    """ this graph and its lower-spin (more pi-bonded) resonances
    """
    def _inc_range(bnd_cap):
        return tuple(range(0, bnd_cap+1))

    add_pi_bonds_ = _partial(_add_pi_bonds, rgr)

    atm_keys = list(_atom_keys(rgr))
    bnd_keys = list(_bond_keys(rgr))
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), atm_keys)
    atm_bnd_keys_lst = _values_by_key(_atom_bond_keys(rgr), atm_keys)
    bnd_caps = _values_by_key(_bond_capacities(rgr), bnd_keys)
    bnd_ord_dct = _bond_orders(rgr)

    def _is_valid(bnd_ord_inc_dct):
        # check if radical decrements are more than radical valences
        def __tally(atm_bnd_keys):
            return sum(_values_by_key(bnd_ord_inc_dct, atm_bnd_keys))
        atm_rad_vlc_decs = tuple(map(__tally, atm_bnd_keys_lst))
        enough_elecs = numpy.all(
            numpy.less_equal(atm_rad_vlc_decs, atm_rad_vlcs))
        # check if all bond orders are less than 4 (should only affect C2)
        bnd_inc_keys = bnd_ord_inc_dct.keys()
        bnd_incs = _values_by_key(bnd_ord_inc_dct, bnd_inc_keys)
        bnd_ords = _values_by_key(bnd_ord_dct, bnd_inc_keys)
        new_bnd_ords = numpy.add(bnd_ords, bnd_incs)
        not_too_many = numpy.all(numpy.less(new_bnd_ords, 4))
        return enough_elecs and not_too_many

    def _bond_value_dictionary(bnd_vals):
        return dict(zip(bnd_keys, bnd_vals))

    bnd_ord_incs_itr = _product(*map(_inc_range, bnd_caps))
    bnd_ord_inc_dct_itr = map(_bond_value_dictionary, bnd_ord_incs_itr)
    bnd_ord_inc_dct_itr = filter(_is_valid, bnd_ord_inc_dct_itr)
    rgrs = tuple(sorted(map(add_pi_bonds_, bnd_ord_inc_dct_itr), key=_frozen))
    return rgrs


def _bond_capacities(rgr):
    """ the number of electron pairs available for further pi-bonding, by bond
    """
    atm_rad_vlc_dct = atom_radical_valences(rgr)

    def _pi_capacities(bnd_key):
        return min(map(atm_rad_vlc_dct.__getitem__, bnd_key))

    bnd_keys = list(_bond_keys(rgr))
    bnd_caps = tuple(map(_pi_capacities, bnd_keys))
    bnd_cap_dct = dict(zip(bnd_keys, bnd_caps))
    return bnd_cap_dct


def _add_pi_bonds(rgr, bnd_ord_inc_dct):
    """ add pi bonds to this graph
    """
    bnd_keys = _bond_keys(rgr)
    assert set(bnd_ord_inc_dct.keys()) <= bnd_keys

    bnd_keys = list(bnd_keys)
    bnd_ords = _values_by_key(_bond_orders(rgr), bnd_keys)
    bnd_ord_incs = _values_by_key(bnd_ord_inc_dct, bnd_keys, fill_val=0)
    new_bnd_ords = numpy.add(bnd_ords, bnd_ord_incs)
    bnd_ord_dct = dict(zip(bnd_keys, new_bnd_ords))
    rgr = _set_bond_orders(rgr, bnd_ord_dct)
    return rgr
