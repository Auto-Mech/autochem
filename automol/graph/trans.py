""" molecular graph transformations (representing transformations)
"""
import numbers
import itertools
from automol import formula
from automol import dict_
import automol.convert.graph
from automol.graph._graph import atom_keys as _atom_keys
from automol.graph._graph import explicit as _explicit
from automol.graph._graph import union as _union
from automol.graph._graph import relabel as _relabel
from automol.graph._graph import connected_components as _connected_components
from automol.graph._graph import full_isomorphism as _full_isomorphism
from automol.graph._graph import add_bonds as _add_bonds
from automol.graph._graph import remove_bonds as _remove_bonds
from automol.graph._graph import atom_neighbor_keys as _atom_neighbor_keys
from automol.graph._graph import (without_stereo_parities as
                                  _without_stereo_parities)
from automol.graph._graph import (add_atom_explicit_hydrogen_keys as
                                  _add_atom_explicit_hydrogen_keys)
from automol.graph._graph import (unsaturated_atom_keys as
                                  _unsaturated_atom_keys)
from automol.graph._graph import atom_stereo_parities as _atom_stereo_parities
from automol.graph._graph import bond_stereo_parities as _bond_stereo_parities
from automol.graph._res import (resonance_dominant_radical_atom_keys as
                                _resonance_dominant_radical_atom_keys)
from automol.graph._stereo import atom_stereo_keys as _atom_stereo_keys
from automol.graph._stereo import bond_stereo_keys as _bond_stereo_keys
from automol.graph._stereo import (stereo_sorted_atom_neighbor_keys as
                                   _stereo_sorted_atom_neighbor_keys)


def from_data(frm_bnd_keys, brk_bnd_keys):
    """ define a transformation from data
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))
    assert all(map(_is_bond_key, frm_bnd_keys))
    assert all(map(_is_bond_key, brk_bnd_keys))
    return (frm_bnd_keys, brk_bnd_keys)


def formed_bond_keys(tra):
    """ keys for bonds that are formed in the transformation
    """
    frm_bnd_keys, _ = tra
    return frm_bnd_keys


def broken_bond_keys(tra):
    """ keys for bonds that are broken in the transformation
    """
    _, brk_bnd_keys = tra
    return brk_bnd_keys


def apply(tra, xgr):
    """ apply this transformation to a graph
    """
    brk_bnd_keys = broken_bond_keys(tra)
    frm_bnd_keys = formed_bond_keys(tra)
    # in case some bonds are broken *and* formed, we subtract the other set
    xgr = _remove_bonds(xgr, brk_bnd_keys - frm_bnd_keys)
    xgr = _add_bonds(xgr, frm_bnd_keys - brk_bnd_keys)
    return xgr


def form_dummy_bonds(tra, xgr):
    """ for each bond formed in this transformation, add a bond with order 0
    """
    bnd_keys = formed_bond_keys(tra) - broken_bond_keys(tra)
    bnd_ord_dct = {bnd_key: 0 for bnd_key in bnd_keys}
    xgr = _add_bonds(xgr, bnd_keys, bnd_ord_dct)
    return xgr


def is_stereo_compatible(tra, sgr1, sgr2):
    """ is this transformation compatible with the applyant/product stereo
    assignments?
    """
    cgr1 = _without_stereo_parities(sgr1)
    cgr2 = _without_stereo_parities(sgr2)
    atm_key_dct = _full_isomorphism(apply(tra, cgr1), cgr2)

    # determine the stereo centers which are preserved in the transformation
    sgr1 = _relabel(sgr1, atm_key_dct)
    atm_keys = sorted(_atom_stereo_keys(sgr1) & _atom_stereo_keys(sgr2))
    bnd_keys = sorted(_bond_stereo_keys(sgr1) & _bond_stereo_keys(sgr2))

    atm_pars1 = dict_.values_by_key(_atom_stereo_parities(sgr1), atm_keys)
    atm_pars2 = dict_.values_by_key(_atom_stereo_parities(sgr2), atm_keys)
    bnd_pars1 = dict_.values_by_key(_bond_stereo_parities(sgr1), bnd_keys)
    bnd_pars2 = dict_.values_by_key(_bond_stereo_parities(sgr2), bnd_keys)

    atm_ngb_keys_dct1 = _atom_neighbor_keys(sgr1)
    atm_ngb_keys_dct2 = _atom_neighbor_keys(sgr2)

    ret = True

    for atm_key, par1, par2 in zip(atm_keys, atm_pars1, atm_pars2):
        atm_ngb_keys1 = _stereo_sorted_atom_neighbor_keys(
            sgr1, atm_key, atm_ngb_keys_dct1[atm_key])
        atm_ngb_keys2 = _stereo_sorted_atom_neighbor_keys(
            sgr2, atm_key, atm_ngb_keys_dct2[atm_key])

        if _permutation_parity(atm_ngb_keys1, atm_ngb_keys2):
            ret &= (par1 == par2)
        else:
            ret &= (par1 != par2)

    for bnd_key, par1, par2 in zip(bnd_keys, bnd_pars1, bnd_pars2):
        atm1_key, atm2_key = bnd_key

        atm1_ngb_key1 = _stereo_sorted_atom_neighbor_keys(
            sgr1, atm1_key, atm_ngb_keys_dct1[atm1_key] - {atm2_key})[0]
        atm2_ngb_key1 = _stereo_sorted_atom_neighbor_keys(
            sgr1, atm2_key, atm_ngb_keys_dct1[atm2_key] - {atm1_key})[0]
        atm1_ngb_key2 = _stereo_sorted_atom_neighbor_keys(
            sgr2, atm1_key, atm_ngb_keys_dct2[atm1_key] - {atm2_key})[0]
        atm2_ngb_key2 = _stereo_sorted_atom_neighbor_keys(
            sgr2, atm2_key, atm_ngb_keys_dct2[atm2_key] - {atm1_key})[0]

        if not ((atm1_ngb_key1 != atm1_ngb_key2) ^
                (atm2_ngb_key1 != atm2_ngb_key2)):
            ret &= (par1 == par2)
        else:
            ret &= (par1 != par2)

    return ret


def _permutation_parity(seq, ref_seq):
    size = len(seq)
    assert sorted(seq) == sorted(ref_seq) and len(set(seq)) == size
    perm = [ref_seq.index(val) for val in seq]

    sgn = 1
    for idx in range(size):
        if perm[idx] != idx:
            sgn *= -1
            swap_idx = perm.index(idx)
            perm[idx], perm[swap_idx] = perm[swap_idx], perm[idx]

    par = (sgn == 1)

    return par


def hydrogen_atom_migration(xgr1, xgr2):
    """ find a hydrogen migration transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    if len(xgrs1) == 1 and len(xgrs2) == 1:
        xgr1, = xgrs1
        xgr2, = xgrs2

        h_atm_key1 = max(_atom_keys(xgr1)) + 1
        h_atm_key2 = max(_atom_keys(xgr2)) + 1

        rad_atm_keys1 = _resonance_dominant_radical_atom_keys(xgr1)
        rad_atm_keys2 = _resonance_dominant_radical_atom_keys(xgr2)
        for atm_key1, atm_key2 in itertools.product(rad_atm_keys1,
                                                    rad_atm_keys2):
            xgr1_h = _add_atom_explicit_hydrogen_keys(
                xgr1, {atm_key1: [h_atm_key1]})
            xgr2_h = _add_atom_explicit_hydrogen_keys(
                xgr2, {atm_key2: [h_atm_key2]})

            inv_atm_key_dct = _full_isomorphism(xgr2_h, xgr1_h)
            if inv_atm_key_dct:
                tra = from_data(
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}])

    return tra


def proton_migration(xgr1, xgr2):
    """ find a proton migration transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)
        
    h_atm_key1 = max(_atom_keys(xgr1)) + 1
    h_atm_key2 = max(_atom_keys(xgr2)) + 1

    if len(xgrs1) == 1 and len(xgrs2) == 1:
        xgr1, = xgrs1
        xgr2, = xgrs2
        atm_keys1 = _unsaturated_atom_keys(xgr1)
        atm_keys2 = _unsaturated_atom_keys(xgr2)
        for atm_key1, atm_key2 in itertools.product(atm_keys1, atm_keys2):
            xgr1_h = _add_atom_explicit_hydrogen_keys(
                xgr1, {atm_key1: [h_atm_key1]})
            xgr2_h = _add_atom_explicit_hydrogen_keys(
                xgr2, {atm_key2: [h_atm_key2]})

            inv_atm_key_dct = _full_isomorphism(xgr2_h, xgr1_h)
            if inv_atm_key_dct:
                tra = from_data(
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}])

    return tra


def beta_scission(xgr1, xgr2):
    """ find a beta scission transformation
    """
    tra = None

    rev_tra = addition(xgr2, xgr1)
    if rev_tra:
        tra = _reverse(rev_tra, xgr2, xgr1)

    return tra


def addition(xgr1, xgr2):
    """ find an addition transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    if len(xgrs1) == 2 and len(xgrs2) == 1:
        x_xgr, y_xgr = xgrs1
        xgr2, = xgrs2
        x_atm_keys = _unsaturated_atom_keys(x_xgr)
        y_atm_keys = _unsaturated_atom_keys(y_xgr)
        for x_atm_key, y_atm_key in itertools.product(x_atm_keys, y_atm_keys):
            xy_xgr = _add_bonds(
                _union(x_xgr, y_xgr), [{x_atm_key, y_atm_key}])

            atm_key_dct = _full_isomorphism(xy_xgr, xgr2)
            if atm_key_dct:
                tra = from_data(frm_bnd_keys=[{x_atm_key, y_atm_key}],
                                brk_bnd_keys=[])

    return tra


def hydrogen_abstraction(xgr1, xgr2):
    """ find an addition transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    ret = formula.reac.argsort_hydrogen_abstraction(
        list(map(automol.convert.graph.formula, xgrs1)),
        list(map(automol.convert.graph.formula, xgrs2)))
    if ret is not None:
        idxs1, idxs2 = ret
        q1h_xgr, q2_xgr = list(map(xgrs1.__getitem__, idxs1))
        q1_xgr, q2h_xgr = list(map(xgrs2.__getitem__, idxs2))
        q1_tra = _partial_hydrogen_abstraction(q1h_xgr, q1_xgr)
        q2_rev_tra = _partial_hydrogen_abstraction(q2h_xgr, q2_xgr)
        if q1_tra and q2_rev_tra:
            xgr1_ = _union(apply(q1_tra, q1h_xgr), q2_xgr)
            xgr2_ = _union(q1_xgr, q2h_xgr)

            q2_tra = _reverse(q2_rev_tra, xgr2_, xgr1_)
            tra = from_data(
                frm_bnd_keys=formed_bond_keys(q2_tra),
                brk_bnd_keys=broken_bond_keys(q1_tra))

    return tra


def _partial_hydrogen_abstraction(qh_xgr, q_xgr):
    tra = None
    h_atm_key = max(_atom_keys(q_xgr)) + 1
    #rad_atm_keys = _resonance_dominant_radical_atom_keys(q_xgr)
    uns_atm_keys = automol.graph.unsaturated_atom_keys(q_xgr)
    for atm_key in uns_atm_keys:
    #for atm_key in rad_atm_keys:
        q_xgr_h = _add_atom_explicit_hydrogen_keys(
            q_xgr, {atm_key: [h_atm_key]})
        inv_atm_key_dct = _full_isomorphism(q_xgr_h, qh_xgr)
        if inv_atm_key_dct:
            brk_bnd_keys = [frozenset(
                {inv_atm_key_dct[atm_key], inv_atm_key_dct[h_atm_key]})]
            tra = from_data(frm_bnd_keys=[], brk_bnd_keys=brk_bnd_keys)

    return tra


def _reverse(tra, xgr1, xgr2):
    frm_bnd_keys = formed_bond_keys(tra)
    brk_bnd_keys = broken_bond_keys(tra)
    atm_key_dct = _full_isomorphism(apply(tra, xgr1), xgr2)
    rev_frm_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in brk_bnd_keys]
    rev_brk_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in frm_bnd_keys]
    rev_tra = from_data(frm_bnd_keys=rev_frm_bnd_keys,
                        brk_bnd_keys=rev_brk_bnd_keys)
    return rev_tra


def _is_atom_key(obj):
    return isinstance(obj, numbers.Integral)


def _is_bond_key(obj):
    return (isinstance(obj, frozenset) and len(obj) == 2 and
            all(map(_is_atom_key, obj)))
