""" resonance graph library
"""
import itertools
import functools
import numpy
from automol import dict_
from automol.graph._graph import frozen as _frozen
from automol.graph._graph import atom_keys as _atom_keys
from automol.graph._graph import bond_keys as _bond_keys
from automol.graph._graph import bond_orders as _bond_orders
from automol.graph._graph import set_bond_orders as _set_bond_orders
from automol.graph._graph import without_bond_orders as _without_bond_orders
from automol.graph._graph import atom_bond_keys as _atom_bond_keys
from automol.graph._graph import atom_neighbor_keys as _atom_neighbor_keys
from automol.graph._graph import (atom_unsaturated_valences as
                                  _atom_unsaturated_valences)
from automol.graph._graph import atom_bond_valences as _atom_bond_valences
from automol.graph._graph import (atom_lone_pair_counts as
                                  _atom_lone_pair_counts)
from automol.graph._graph import (maximum_spin_multiplicity as
                                  _maximum_spin_multiplicity)
from automol.graph._graph import explicit as _explicit
from automol.graph._graph import (atom_explicit_hydrogen_valences as
                                  _atom_explicit_hydrogen_valences)


# atom properties
def atom_hybridizations(rgr):
    """ atom hybridizations, by atom
    """
    atm_keys = list(_atom_keys(rgr))
    atm_unsat_vlc_dct = _atom_unsaturated_valences(rgr, bond_order=True)
    atm_bnd_vlc_dct = _atom_bond_valences(rgr, bond_order=False)     # note!!
    atm_unsat_vlcs = numpy.array(
        dict_.values_by_key(atm_unsat_vlc_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(
        dict_.values_by_key(_atom_lone_pair_counts(rgr), atm_keys))
    atm_hybs = atm_unsat_vlcs + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def resonance_dominant_atom_hybridizations(rgr):
    """ resonance-dominant atom hybridizations, by atom
    """
    atm_keys = list(_atom_keys(rgr))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def resonance_dominant_atom_centered_cumulene_keys(rgr):
    """ resonance dominant keys for atom-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}), cent_atm_key)
    where the first pair contains the sp2 atoms at the cumulene ends and
    `cent_atm_key` is the key of the central atom
    """
    cum_chains = _cumulene_chains(rgr)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 1:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}),
                 cum_chain[size // 2])
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


def resonance_dominant_bond_centered_cumulene_keys(rgr):
    """ resonance dominant keys for bond-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}),
         frozenset({cent_atm_key1, cent_atm_key2}))
    where the first pair contains the sp2 atoms at the cumulene ends and the
    second pair is the bond key for the central bond
    """
    cum_chains = _cumulene_chains(rgr)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 0:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}),
                 frozenset({cum_chain[size // 2 - 1], cum_chain[size // 2]}))
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


def _cumulene_chains(rgr):
    atm_hyb_dct = resonance_dominant_atom_hybridizations(rgr)
    sp1_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)

    atm_ngb_keys_dct = _atom_neighbor_keys(rgr)

    def _cumulene_chain(chain):
        ret = None
        atm_key = chain[-1]
        next_atm_keys = atm_ngb_keys_dct[atm_key] - {chain[-2]}
        if next_atm_keys:
            assert len(next_atm_keys) == 1
            next_atm_key, = next_atm_keys
            if next_atm_key in sp1_atm_keys:
                chain.append(next_atm_key)
                ret = _cumulene_chain(chain)
            elif next_atm_key in sp2_atm_keys:
                chain.append(next_atm_key)
                ret = chain
        return ret

    cum_chains = []
    for atm_key in sp2_atm_keys:
        sp1_atm_ngb_keys = atm_ngb_keys_dct[atm_key] & sp1_atm_keys
        chains = [[atm_key, atm_ngb_key] for atm_ngb_key in sp1_atm_ngb_keys]
        for chain in chains:
            cum_chain = _cumulene_chain(chain)
            if cum_chain is not None:
                cum_chains.append(cum_chain)

    cum_chains = tuple(map(tuple, cum_chains))
    return cum_chains


def resonance_dominant_radical_atom_keys(rgr):
    """ resonance-dominant radical atom keys

    (keys of resonance-dominant radical sites)
    """
    atm_keys = list(_atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(_atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def sing_res_dom_radical_atom_keys(rgr):
    """ resonance-dominant radical atom keys,for one resonance
    """
    atm_keys = list(_atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(_atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    first_atm_rad_val = [atm_rad_vlcs_by_res[0]]
    atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*first_atm_rad_val)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


# bond properties
def resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(_bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(_bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    bnd_ords_lst = list(map(frozenset, zip(*bnd_ords_by_res)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


def one_resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(_bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(_bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    first_bnd_ords = [bnd_ords_by_res[0]]
    bnd_ords_lst = list(map(frozenset, zip(*first_bnd_ords)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


def resonance_avg_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(_bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(_bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    nres = len(bnd_ords_by_res)
    bnd_ords_lst = zip(*bnd_ords_by_res)
    avg_bnd_ord_lst = [sum(bnd_ords)/nres for bnd_ords in bnd_ords_lst]
    avg_bnd_ord_dct = dict(zip(bnd_keys, avg_bnd_ord_lst))
    return avg_bnd_ord_dct


# transformations
def dominant_resonance(rgr):
    """ *a* dominant (minimum spin/maximum pi) resonance graph
    """
    return next(iter(dominant_resonances(rgr)))


def dominant_resonances(rgr):
    """ all dominant (minimum spin/maximum pi) resonance graphs
    """
    rgrs = resonances(rgr)
    mult_min = min(map(_maximum_spin_multiplicity, rgrs))
    dom_rgrs = tuple(
        rgr for rgr in rgrs if _maximum_spin_multiplicity(rgr) == mult_min)
    return dom_rgrs


def resonances(rgr):
    """ all resonance graphs with this connectivity
    """
    return subresonances(_without_bond_orders(rgr))


def subresonances(rgr):
    """ this connected graph and its lower-spin (more pi-bonded) resonances
    """
    def _inc_range(bnd_cap):
        return tuple(range(0, bnd_cap+1))

    add_pi_bonds_ = functools.partial(_add_pi_bonds, rgr)

    atm_keys = list(_atom_keys(rgr))
    bnd_keys = list(_bond_keys(rgr))
    atm_unsat_vlcs = dict_.values_by_key(
        _atom_unsaturated_valences(rgr), atm_keys)
    atm_bnd_keys_lst = dict_.values_by_key(_atom_bond_keys(rgr), atm_keys)
    bnd_caps = dict_.values_by_key(_bond_capacities(rgr), bnd_keys)
    bnd_ord_dct = _bond_orders(rgr)

    def _is_valid(bnd_ord_inc_dct):
        # check if pi bonds exceed unsaturated valences
        def __tally(atm_bnd_keys):
            return sum(dict_.values_by_key(bnd_ord_inc_dct, atm_bnd_keys))
        atm_unsat_vlc_decs = tuple(map(__tally, atm_bnd_keys_lst))
        enough_elecs = numpy.all(
            numpy.less_equal(atm_unsat_vlc_decs, atm_unsat_vlcs))
        # check if all bond orders are less than 4 (should only affect C2)
        bnd_inc_keys = bnd_ord_inc_dct.keys()
        bnd_incs = dict_.values_by_key(bnd_ord_inc_dct, bnd_inc_keys)
        bnd_ords = dict_.values_by_key(bnd_ord_dct, bnd_inc_keys)
        new_bnd_ords = numpy.add(bnd_ords, bnd_incs)
        not_too_many = numpy.all(numpy.less(new_bnd_ords, 4))
        return enough_elecs and not_too_many

    def _bond_value_dictionary(bnd_vals):
        return dict(zip(bnd_keys, bnd_vals))

    bnd_ord_incs_itr = itertools.product(*map(_inc_range, bnd_caps))
    bnd_ord_inc_dct_itr = map(_bond_value_dictionary, bnd_ord_incs_itr)
    bnd_ord_inc_dct_itr = filter(_is_valid, bnd_ord_inc_dct_itr)
    rgrs = tuple(sorted(map(add_pi_bonds_, bnd_ord_inc_dct_itr), key=_frozen))
    return rgrs


def _bond_capacities(rgr):
    """ the number of electron pairs available for further pi-bonding, by bond
    """
    atm_unsat_vlc_dct = _atom_unsaturated_valences(rgr)

    def _pi_capacities(bnd_key):
        return min(map(atm_unsat_vlc_dct.__getitem__, bnd_key))

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
    bnd_ords = dict_.values_by_key(_bond_orders(rgr), bnd_keys)
    bnd_ord_incs = dict_.values_by_key(bnd_ord_inc_dct, bnd_keys, fill_val=0)
    new_bnd_ords = numpy.add(bnd_ords, bnd_ord_incs)
    bnd_ord_dct = dict(zip(bnd_keys, new_bnd_ords))
    rgr = _set_bond_orders(rgr, bnd_ord_dct)
    return rgr


# other utilities
def rotational_bond_keys(xgr, with_h_rotors=True):
    """ determine rotational bonds in this molecular graph
    """
    xgr = _explicit(xgr)
    atm_bnd_vlc_dct = _atom_bond_valences(xgr, bond_order=False)
    atm_exp_hyd_vlc_dct = _atom_explicit_hydrogen_valences(xgr)
    res_dom_bnd_ords_dct = resonance_dominant_bond_orders(xgr)

    bnd_keys = []
    for bnd_key, bnd_ords in res_dom_bnd_ords_dct.items():
        if all(bnd_ord <= 1 for bnd_ord in bnd_ords):
            atm_keys = list(bnd_key)
            bnd_ord = min(bnd_ords)
            rot_vlcs = numpy.array(
                list(map(atm_bnd_vlc_dct.__getitem__, atm_keys)))
            rot_vlcs -= bnd_ord
            if not with_h_rotors:
                atm_exp_hyd_vlcs = numpy.array(list(
                    map(atm_exp_hyd_vlc_dct.__getitem__, atm_keys)))
                rot_vlcs -= atm_exp_hyd_vlcs
            if all(rot_vlcs):
                bnd_keys.append(bnd_key)

    return frozenset(bnd_keys)
