""" molecular graph reactions
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
    """ define a reaction from data
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))
    assert all(map(_is_bond_key, frm_bnd_keys))
    assert all(map(_is_bond_key, brk_bnd_keys))
    return (frm_bnd_keys, brk_bnd_keys)


def formed_bond_keys(rxn):
    """ keys for bonds that are formed in the reaction
    """
    frm_bnd_keys, _ = rxn
    return frm_bnd_keys


def broken_bond_keys(rxn):
    """ keys for bonds that are broken in the reaction
    """
    _, brk_bnd_keys = rxn
    return brk_bnd_keys


def react(rxn, xgr):
    """ apply this transformation to a graph
    """
    brk_bnd_keys = broken_bond_keys(rxn)
    frm_bnd_keys = formed_bond_keys(rxn)
    # in case some bonds are broken *and* formed, we subtract the other set
    xgr = _remove_bonds(xgr, brk_bnd_keys - frm_bnd_keys)
    xgr = _add_bonds(xgr, frm_bnd_keys - brk_bnd_keys)
    return xgr


def is_stereo_compatible(rxn, sgr1, sgr2):
    """ is this reaction compatible with the reactant/product stereo assignments?
    """
    cgr1 = _without_stereo_parities(sgr1)
    cgr2 = _without_stereo_parities(sgr2)
    atm_key_dct = _full_isomorphism(react(rxn, cgr1), cgr2)

    # determine the stereo centers which are preserved in the reaction
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


def hydrogen_migration(xgr1, xgr2):
    """ find a hydrogen migration reaction
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    rxn = None
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
                rxn = from_data(
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}])

    return rxn


def beta_scission(xgr1, xgr2):
    """ find a beta scission reaction
    """
    rxn = None

    rev_rxn = addition(xgr2, xgr1)
    if rev_rxn:
        rxn = _reverse(rev_rxn, xgr2, xgr1)

    return rxn


def addition(xgr1, xgr2):
    """ find an addition reaction
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    rxn = None
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
                rxn = from_data(frm_bnd_keys=[{x_atm_key, y_atm_key}],
                                brk_bnd_keys=[])

    return rxn


def hydrogen_abstraction(xgr1, xgr2):
    """ find an addition reaction
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    rxn = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    if len(xgrs1) == 2 and len(xgrs2) == 2:
        for xgrs1_ in itertools.permutations(xgrs1):
            for xgrs2_ in itertools.permutations(xgrs2):
                fmls1 = list(map(automol.convert.graph.formula, xgrs1_))
                fmls2 = list(map(automol.convert.graph.formula, xgrs2_))
                if (formula.add_hydrogen(fmls2[0]) == fmls1[0] and
                        formula.add_hydrogen(fmls1[1]) == fmls2[1]):
                    q1h_xgr, q2_xgr = xgrs1_
                    q1_xgr, q2h_xgr = xgrs2_
                    q1_rxn = _partial_hydrogen_abstraction(q1h_xgr, q1_xgr)
                    q2_rev_rxn = _partial_hydrogen_abstraction(q2h_xgr, q2_xgr)
                    if q1_rxn and q2_rev_rxn:
                        xgr1_ = _union(react(q1_rxn, q1h_xgr), q2_xgr)
                        xgr2_ = _union(q1_xgr, q2h_xgr)

                        q2_rxn = _reverse(q2_rev_rxn, xgr2_, xgr1_)
                        rxn = from_data(
                            frm_bnd_keys=formed_bond_keys(q2_rxn),
                            brk_bnd_keys=broken_bond_keys(q1_rxn))

    return rxn


def _partial_hydrogen_abstraction(qh_xgr, q_xgr):
    rxn = None
    h_atm_key = max(_atom_keys(q_xgr)) + 1
    rad_atm_keys = _resonance_dominant_radical_atom_keys(q_xgr)
    for atm_key in rad_atm_keys:
        q_xgr_h = _add_atom_explicit_hydrogen_keys(
            q_xgr, {atm_key: [h_atm_key]})
        inv_atm_key_dct = _full_isomorphism(q_xgr_h, qh_xgr)
        if inv_atm_key_dct:
            brk_bnd_keys = [frozenset(
                {inv_atm_key_dct[atm_key], inv_atm_key_dct[h_atm_key]})]
            rxn = from_data(frm_bnd_keys=[], brk_bnd_keys=brk_bnd_keys)

    return rxn


def _reverse(rxn, xgr1, xgr2):
    frm_bnd_keys = formed_bond_keys(rxn)
    brk_bnd_keys = broken_bond_keys(rxn)
    atm_key_dct = _full_isomorphism(react(rxn, xgr1), xgr2)
    rev_frm_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in brk_bnd_keys]
    rev_brk_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in frm_bnd_keys]
    rev_rxn = from_data(frm_bnd_keys=rev_frm_bnd_keys,
                        brk_bnd_keys=rev_brk_bnd_keys)
    return rev_rxn


def _is_atom_key(obj):
    return isinstance(obj, numbers.Integral)


def _is_bond_key(obj):
    return (isinstance(obj, frozenset) and len(obj) == 2 and
            all(map(_is_atom_key, obj)))
