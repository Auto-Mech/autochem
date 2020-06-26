""" molecular graph transformations (representing reactions)
"""
import numbers
from automol import dict_
from automol.graph._graph import relabel as _relabel
from automol.graph._graph import full_isomorphism
from automol.graph._graph import atom_neighbor_keys
from automol.graph._graph import add_bonds
from automol.graph._graph import remove_bonds
from automol.graph._stereo import atom_stereo_keys
from automol.graph._stereo import bond_stereo_keys
from automol.graph._graph import atom_stereo_parities
from automol.graph._graph import bond_stereo_parities
from automol.graph._graph import without_stereo_parities
from automol.graph._stereo import stereo_sorted_atom_neighbor_keys


def from_data(frm_bnd_keys, brk_bnd_keys):
    """ define a transformation from data
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))
    assert all(map(_is_bond_key, frm_bnd_keys))
    assert all(map(_is_bond_key, brk_bnd_keys))
    assert not frm_bnd_keys & brk_bnd_keys
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


def relabel(tra, atm_key_dct):
    """ relabel the atom keys in the transformation
    """
    def _relabel_bond_key(bnd_key):
        return frozenset(map(atm_key_dct.__getitem__, bnd_key))

    frm_bnd_keys = list(map(_relabel_bond_key, formed_bond_keys(tra)))
    brk_bnd_keys = list(map(_relabel_bond_key, broken_bond_keys(tra)))

    return from_data(frm_bnd_keys, brk_bnd_keys)


def apply(tra, xgr):
    """ apply this transformation to a graph
    """
    brk_bnd_keys = broken_bond_keys(tra)
    frm_bnd_keys = formed_bond_keys(tra)
    # in case some bonds are broken *and* formed, we subtract the other set
    xgr = remove_bonds(xgr, brk_bnd_keys - frm_bnd_keys)
    xgr = add_bonds(xgr, frm_bnd_keys - brk_bnd_keys)
    return xgr


def reverse(tra, xgr1, xgr2):
    """ reverse a transformation to get the one taking products into reactants
    """
    frm_bnd_keys = formed_bond_keys(tra)
    brk_bnd_keys = broken_bond_keys(tra)
    atm_key_dct = full_isomorphism(apply(tra, xgr1), xgr2)
    rev_frm_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in brk_bnd_keys]
    rev_brk_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
                        for bnd_key in frm_bnd_keys]
    rev_tra = from_data(frm_bnd_keys=rev_frm_bnd_keys,
                        brk_bnd_keys=rev_brk_bnd_keys)
    return rev_tra


def is_stereo_compatible(tra, sgr1, sgr2):
    """ is this transformation compatible with the applyant/product stereo
    assignments?
    """
    cgr1 = without_stereo_parities(sgr1)
    cgr2 = without_stereo_parities(sgr2)
    atm_key_dct = full_isomorphism(apply(tra, cgr1), cgr2)

    # determine the stereo centers which are preserved in the transformation
    sgr1 = _relabel(sgr1, atm_key_dct)
    atm_keys = sorted(atom_stereo_keys(sgr1) & atom_stereo_keys(sgr2))
    bnd_keys = sorted(bond_stereo_keys(sgr1) & bond_stereo_keys(sgr2))

    atm_pars1 = dict_.values_by_key(atom_stereo_parities(sgr1), atm_keys)
    atm_pars2 = dict_.values_by_key(atom_stereo_parities(sgr2), atm_keys)
    bnd_pars1 = dict_.values_by_key(bond_stereo_parities(sgr1), bnd_keys)
    bnd_pars2 = dict_.values_by_key(bond_stereo_parities(sgr2), bnd_keys)

    atm_ngb_keys_dct1 = atom_neighbor_keys(sgr1)
    atm_ngb_keys_dct2 = atom_neighbor_keys(sgr2)

    ret = True

    for atm_key, par1, par2 in zip(atm_keys, atm_pars1, atm_pars2):
        atm_ngb_keys1 = stereo_sorted_atom_neighbor_keys(
            sgr1, atm_key, atm_ngb_keys_dct1[atm_key])
        atm_ngb_keys2 = stereo_sorted_atom_neighbor_keys(
            sgr2, atm_key, atm_ngb_keys_dct2[atm_key])

        if _permutation_parity(atm_ngb_keys1, atm_ngb_keys2):
            ret &= (par1 == par2)
        else:
            ret &= (par1 != par2)

    for bnd_key, par1, par2 in zip(bnd_keys, bnd_pars1, bnd_pars2):
        atm1_key, atm2_key = bnd_key

        atm1_ngb_key1 = stereo_sorted_atom_neighbor_keys(
            sgr1, atm1_key, atm_ngb_keys_dct1[atm1_key] - {atm2_key})[0]
        atm2_ngb_key1 = stereo_sorted_atom_neighbor_keys(
            sgr1, atm2_key, atm_ngb_keys_dct1[atm2_key] - {atm1_key})[0]
        atm1_ngb_key2 = stereo_sorted_atom_neighbor_keys(
            sgr2, atm1_key, atm_ngb_keys_dct2[atm1_key] - {atm2_key})[0]
        atm2_ngb_key2 = stereo_sorted_atom_neighbor_keys(
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


def _is_bond_key(obj):
    return (isinstance(obj, frozenset) and len(obj) == 2 and
            all(map(_is_atom_key, obj)))


def _is_atom_key(obj):
    return isinstance(obj, numbers.Integral)
