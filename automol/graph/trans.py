""" molecular graph transformations (representing reactions)
"""
import itertools
import numbers
import yaml
from automol.util import dict_
from automol import formula
from automol import convert
from automol.graph._graph_dep import relabel as _relabel
from automol.graph._graph_dep import atom_neighbor_keys
from automol.graph._graph_dep import atoms
from automol.graph._graph_dep import bond_keys
from automol.graph._graph_dep import atom_lone_pair_counts
from automol.graph._graph_dep import add_bonds
from automol.graph._graph_dep import remove_bonds
from automol.graph._graph import (
    connected_components as _connected_components)
from automol.graph._graph_dep import explicit as _explicit
from automol.graph._graph_dep import union as _union
from automol.graph._graph_dep import (
    unsaturated_atom_keys as _unsaturated_atom_keys)
from automol.graph._graph_dep import atom_keys as _atom_keys
from automol.graph._graph_dep import reverse as _reverse
from automol.graph._graph_dep import (
    resonance_dominant_radical_atom_key as
    _resonance_dominant_radical_atom_keys)
from automol.graph._graph_dep import (
    add_atom_explicit_hydrogen_keys as
    _add_atom_explicit_hydrogen_keys)
from automol.graph._graph_dep import atom_stereo_keys
from automol.graph._graph_dep import bond_stereo_keys
from automol.graph._graph_dep import atom_stereo_parities
from automol.graph._graph_dep import bond_stereo_parities
from automol.graph._graph_dep import without_stereo_parities
from automol.graph._graph_dep import stereo_sorted_atom_neighbor_keys
from automol.graph._graph import full_isomorophism as _full_isomorphism
from automol.reac._find import _partial_hydrogen_abstraction


# old
def old_from_data(frm_bnd_keys, brk_bnd_keys):
    """ define a transformation from data
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))
    assert all(map(_is_bond_key, frm_bnd_keys))
    assert all(map(_is_bond_key, brk_bnd_keys))
    assert not frm_bnd_keys & brk_bnd_keys
    return (frm_bnd_keys, brk_bnd_keys)


def old_formed_bond_keys(tra):
    """ keys for bonds that are formed in the transformation
    """
    frm_bnd_keys, _ = tra
    return frm_bnd_keys


def old_broken_bond_keys(tra):
    """ keys for bonds that are broken in the transformation
    """
    _, brk_bnd_keys = tra
    return brk_bnd_keys


def old_string(tra):
    """ write the transformation to a string
    """
    def _encode_bond(bnd_key):
        atm1_key, atm2_key = bnd_key
        bnd_str = '{}-{}'.format(atm1_key+1, atm2_key+1)
        return bnd_str

    frm_bnd_keys = sorted(map(sorted, old_formed_bond_keys(tra)))
    brk_bnd_keys = sorted(map(sorted, old_broken_bond_keys(tra)))

    if any(frm_bnd_keys):
        frm_bnd_strs = list(map(_encode_bond, frm_bnd_keys))
    else:
        frm_bnd_strs = None

    if any(brk_bnd_keys):
        brk_bnd_strs = list(map(_encode_bond, brk_bnd_keys))
    else:
        brk_bnd_strs = None

    tra_dct = {'bonds formed': frm_bnd_strs,
               'bonds broken': brk_bnd_strs}

    tra_str = yaml.dump(tra_dct, default_flow_style=None, sort_keys=False)
    return tra_str


def old_from_string(tra_str):
    """ read the transformation from a string
    """
    def _decode_bond(bnd_str):
        atm1_key, atm2_key = map(int, bnd_str.split('-'))
        bnd_key = frozenset({atm1_key-1, atm2_key-1})
        return bnd_key

    tra_dct = yaml.load(tra_str, Loader=yaml.FullLoader)
    frm_bnd_strs = tra_dct['bonds formed']
    brk_bnd_strs = tra_dct['bonds broken']

    if frm_bnd_strs is not None:
        frm_bnd_keys = frozenset(map(_decode_bond, frm_bnd_strs))
    else:
        frm_bnd_keys = frozenset({})

    if brk_bnd_strs is not None:
        brk_bnd_keys = frozenset(map(_decode_bond, brk_bnd_strs))
    else:
        brk_bnd_keys = frozenset({})

    tra = old_from_data(frm_bnd_keys, brk_bnd_keys)

    return tra


# New
def from_data(rxn_class, frm_bnd_keys, brk_bnd_keys):
    """ define a transformation from data
    """
    frm_bnd_keys = frozenset(map(frozenset, frm_bnd_keys))
    brk_bnd_keys = frozenset(map(frozenset, brk_bnd_keys))
    assert all(map(_is_bond_key, frm_bnd_keys))
    assert all(map(_is_bond_key, brk_bnd_keys))
    assert not frm_bnd_keys & brk_bnd_keys
    return (rxn_class, frm_bnd_keys, brk_bnd_keys)


def reaction_class(tra):
    """ string describing the reaction class
    """
    rxn_class, _, _ = tra
    # assert par.is_reaction_class(rxn_class), (
    #     '{} is not an allowed reaction class'.format(rxn_class)
    # )
    return rxn_class


def formed_bond_keys(tra):
    """ keys for bonds that are formed in the transformation
    """
    _, frm_bnd_keys, _ = tra
    return frm_bnd_keys


def broken_bond_keys(tra):
    """ keys for bonds that are broken in the transformation
    """
    _, _, brk_bnd_keys = tra
    return brk_bnd_keys


def string(tra):
    """ write the transformation to a string
    """
    def _encode_bond(bnd_key):
        atm1_key, atm2_key = bnd_key
        bnd_str = '{}-{}'.format(atm1_key+1, atm2_key+1)
        return bnd_str

    rxn_class = reaction_class(tra)
    frm_bnd_keys = sorted(map(sorted, formed_bond_keys(tra)))
    brk_bnd_keys = sorted(map(sorted, broken_bond_keys(tra)))

    if any(frm_bnd_keys):
        frm_bnd_strs = list(map(_encode_bond, frm_bnd_keys))
    else:
        frm_bnd_strs = None

    if any(brk_bnd_keys):
        brk_bnd_strs = list(map(_encode_bond, brk_bnd_keys))
    else:
        brk_bnd_strs = None

    tra_dct = {'reaction class': rxn_class,
               'bonds formed': frm_bnd_strs,
               'bonds broken': brk_bnd_strs}

    tra_str = yaml.dump(tra_dct, sort_keys=False)
    return tra_str


def from_string(tra_str):
    """ read the transformation from a string
    """
    def _decode_bond(bnd_str):
        atm1_key, atm2_key = map(int, bnd_str.split('-'))
        bnd_key = frozenset({atm1_key-1, atm2_key-1})
        return bnd_key

    tra_dct = yaml.load(tra_str, Loader=yaml.FullLoader)
    rxn_class = tra_dct['reaction class']
    frm_bnd_strs = tra_dct['bonds formed']
    brk_bnd_strs = tra_dct['bonds broken']

    if frm_bnd_strs is not None:
        frm_bnd_keys = frozenset(map(_decode_bond, frm_bnd_strs))
    else:
        frm_bnd_keys = frozenset({})

    if brk_bnd_strs is not None:
        brk_bnd_keys = frozenset(map(_decode_bond, brk_bnd_strs))
    else:
        brk_bnd_keys = frozenset({})

    tra = from_data(rxn_class, frm_bnd_keys, brk_bnd_keys)

    return tra


def relabel(tra, atm_key_dct):
    """ relabel the atom keys in the transformation
    """
    def _relabel_bond_key(bnd_key):
        return frozenset(map(atm_key_dct.__getitem__, bnd_key))

    rxn_class = reaction_class(tra)
    frm_bnd_keys = list(map(_relabel_bond_key, formed_bond_keys(tra)))
    brk_bnd_keys = list(map(_relabel_bond_key, broken_bond_keys(tra)))

    return from_data(rxn_class, frm_bnd_keys, brk_bnd_keys)


def apply(tra, xgr):
    """ apply this transformation to a graph
    """
    brk_bnd_keys = broken_bond_keys(tra)
    frm_bnd_keys = formed_bond_keys(tra)
    # in case some bonds are broken *and* formed, we subtract the other set
    xgr = remove_bonds(xgr, brk_bnd_keys - frm_bnd_keys)
    xgr = add_bonds(xgr, frm_bnd_keys - brk_bnd_keys)
    return xgr


# NEED TO FIX THIS FUNCTION (NEED REVERSE REaCTION ID CODE)
# def reverse(tra, xgr1, xgr2):
#     """ reverse a transformation to get the one taking products intoreactants
#     """
#     rxn_class = reaction_class(tra)
#     frm_bnd_keys = formed_bond_keys(tra)
#     brk_bnd_keys = broken_bond_keys(tra)
#     atm_key_dct = full_isomorphism(apply(tra, xgr1), xgr2)
#
#     rev_rxn_class = reverse_class(rxn_class)
#     rev_frm_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
#                         for bnd_key in brk_bnd_keys]
#     rev_brk_bnd_keys = [frozenset(map(atm_key_dct.__getitem__, bnd_key))
#                         for bnd_key in frm_bnd_keys]
#     rev_tra = from_data(rxn_class=rev_rxn_class,
#                         frm_bnd_keys=rev_frm_bnd_keys,
#                         brk_bnd_keys=rev_brk_bnd_keys)
#     return rev_tra
#
#
def is_stereo_compatible(tra, sgr1, sgr2):
    """ is this transformation compatible with the reactant/product stereo
    assignments?
    """
    cgr1 = without_stereo_parities(sgr1)
    cgr2 = without_stereo_parities(sgr2)
    atm_key_dct = _full_isomorphism(apply(tra, cgr1), cgr2)

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


def hydrogen_atom_migration(xgr1, xgr2):
    """ find a hydrogen migration transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tras = []
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
                tras.append(from_data(
                    'hydrogen migration',
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}]))

    if len(tras) < 1:
        tras = None
    return tras


def proton_migration(xgr1, xgr2):
    """ find a proton migration transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tras = []
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
                tras.append(from_data(
                    'proton migration',
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}]))
    if len(tras) < 1:
        tras = None
    return tras


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
        # xgeo = geometry(xgr2)
        for x_atm_key, y_atm_key in itertools.product(x_atm_keys, y_atm_keys):
            xy_xgr = add_bonds(
                _union(x_xgr, y_xgr), [{x_atm_key, y_atm_key}])

            # xgeo = geometry(xy_xgr)
            atm_key_dct = _full_isomorphism(xy_xgr, xgr2)
            if atm_key_dct:
                tra = from_data(
                    'addition',
                    frm_bnd_keys=[{x_atm_key, y_atm_key}],
                    brk_bnd_keys=[])

    return tra


def elimination(xgr1, xgr2):
    """identifies elimination reactions
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)
    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)
    tras = []
    if len(xgrs1) == 1 and len(xgrs2) == 2:
        atms = atoms(xgr1)
        neighs = atom_neighbor_keys(xgr1)
        bnds = bond_keys(xgr1)
        radicals = _resonance_dominant_radical_atom_keys(xgr1)
        lonepairs = atom_lone_pair_counts(xgr1)
        for atmi in atms:
            i_neighs = neighs[atmi]
            for atmj in i_neighs:
                bnd_break_key_ij = _get_bnd_key(atmi, atmj, bnds)
                new_xgr = remove_bonds(xgr1, [bnd_break_key_ij])
                new_xgrs = _connected_components(new_xgr)
                if len(new_xgrs) == 2:
                    xgra, xgrb = new_xgrs
                    atmsa = atoms(xgra)
                    if atmi not in atmsa.keys():
                        xgrb, xgra = xgra, xgrb
                        atmsa = atoms(xgra)
                    neighsa = atom_neighbor_keys(xgra)
                    atmsb = atoms(xgrb)
                    neighs_i = neighsa[atmi]
                    for atmk in atmsb:
                        if atmk in radicals:
                            for atml in neighs_i:
                                neighs_l = neighsa[atml]
                                if atml != atmj:
                                    bnd_break_key_il = _get_bnd_key(
                                        atmi, atml, bnds)
                                    bnd_form_key_kl = frozenset({atmk, atml})
                                    newnew_xgr = remove_bonds(
                                        new_xgr, [bnd_break_key_il])
                                    newnew_xgr = add_bonds(
                                        newnew_xgr, [bnd_form_key_kl])
                                    atm_key_dct = _full_isomorphism(
                                        newnew_xgr, xgr2)
                                    if atm_key_dct:
                                        tra = [
                                            [bnd_form_key_kl],
                                            [bnd_break_key_ij,
                                             bnd_break_key_il]]
                                        return tra
                                for atmm in neighs_l:
                                    if atmm != atmi:
                                        bnd_break_key_lm = _get_bnd_key(
                                            atml, atmm, bnds)
                                        bnd_form_key_km = frozenset(
                                            {atmk, atmm})
                                        newnew_xgr = remove_bonds(
                                            new_xgr, [bnd_break_key_lm])
                                        newnew_xgr = add_bonds(
                                            newnew_xgr, [bnd_form_key_km])
                                        atm_key_dct = _full_isomorphism(
                                            newnew_xgr, xgr2)
                                        if atm_key_dct:
                                            tras.append([
                                                [bnd_form_key_km],
                                                [bnd_break_key_ij,
                                                 bnd_break_key_lm]])
        for atmi in atms:
            i_neighs = neighs[atmi]
            print('atmi test:', atmi)
            print('i_neighs test:', i_neighs)
            for atmj in i_neighs:
                bnd_break_key_ij = _get_bnd_key(atmi, atmj, bnds)
                new_xgr = remove_bonds(xgr1, [bnd_break_key_ij])
                new_xgrs = _connected_components(new_xgr)
                if len(new_xgrs) == 2:
                    xgra, xgrb = new_xgrs
                    atmsa = atoms(xgra)
                    if atmi not in atmsa.keys():
                        xgrb, xgra = xgra, xgrb
                        atmsa = atoms(xgra)
                    neighsa = atom_neighbor_keys(xgra)
                    atmsb = atoms(xgrb)
                    neighs_i = neighsa[atmi]
                    print('len atmsb test:', len(atmsb))
                    for atmk in atmsb:
                        if lonepairs[atmk] > 0 or len(atmsb) == 1:
                            # if lonepairs[atmk] > 0:
                            for atml in neighs_i:
                                neighs_l = neighsa[atml]
                                if atml != atmj:
                                    bnd_break_key_il = _get_bnd_key(
                                        atmi, atml, bnds)
                                    bnd_form_key_kl = frozenset({atmk, atml})
                                    newnew_xgr = remove_bonds(
                                        new_xgr, [bnd_break_key_il])
                                    newnew_xgr = add_bonds(
                                        newnew_xgr, [bnd_form_key_kl])
                                    atm_key_dct = _full_isomorphism(
                                        newnew_xgr, xgr2)
                                    if atm_key_dct:
                                        tra = [
                                            [bnd_form_key_kl],
                                            [bnd_break_key_ij,
                                             bnd_break_key_il]]
                                        return tra
                                for atmm in neighs_l:
                                    if atmm != atmi:
                                        bnd_break_key_lm = _get_bnd_key(
                                            atml, atmm, bnds)
                                        bnd_form_key_km = frozenset(
                                            {atmk, atmm})
                                        newnew_xgr = remove_bonds(
                                            new_xgr, [bnd_break_key_lm])
                                        newnew_xgr = add_bonds(
                                            newnew_xgr, [bnd_form_key_km])
                                        atm_key_dct = _full_isomorphism(
                                            newnew_xgr, xgr2)
                                        if atm_key_dct:
                                            tras.append([
                                                [bnd_form_key_km],
                                                [bnd_break_key_ij,
                                                 bnd_break_key_lm]])
    if len(tras) < 1:
        tras = None
    return tras


def substitution(xgr1, xgr2):
    """identifies substitution reactions
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    idxs = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    print('len xgrs test:', len(xgrs1), len(xgrs2))
    print('xgrs test:', xgrs1, xgrs2)

    if len(xgrs1) == 2 and len(xgrs2) == 2:
        xgra, xgrb = xgrs1
        xgrc, xgrd = xgrs2
        atmsa = atoms(xgra)
        neighsa = atom_neighbor_keys(xgra)
        bndsa = bond_keys(xgra)
        atmsb = atoms(xgrb)
        neighsb = atom_neighbor_keys(xgrb)
        bndsb = bond_keys(xgrb)
        atmsc = atoms(xgrc)
        neighsc = atom_neighbor_keys(xgrc)
        bndsc = bond_keys(xgrc)
        atmsd = atoms(xgrd)
        neighsd = atom_neighbor_keys(xgrd)
        bndsd = bond_keys(xgrd)
        tra = _substitution(atmsa, neighsa, bndsa, atmsb, xgr1, xgr2)
        # tra = _substitution(
        # atmsa, neighsa, bndsa, atmsb, neighsb, xgr1, xgr2)
        idxs = [[0, 1], [0, 1]]
        if not tra:
            tra = _substitution(
                atmsb, neighsb, bndsb, atmsa, xgr1, xgr2)
            # atmsb, neighsb, bndsb, atmsa, neighsa, xgr1, xgr2)
            idxs = [[0, 1], [1, 0]]
            if not tra:
                tra = _substitution(
                    atmsc, neighsc, bndsc, atmsd, xgr2, xgr1)
                # atmsc, neighsc, bndsc, atmsd, neighsd, xgr2, xgr1)
                idxs = [[1, 0], [0, 1]]
                if not tra:
                    tra = _substitution(
                        atmsd, neighsd, bndsd, atmsc, xgr2, xgr1)
                    # atmsd, neighsd, bndsd, atmsc, neighsc, xgr2, xgr1)
                    idxs = [[1, 0], [1, 0]]
                    if not tra:
                        idxs = None

        # return not substitution for radical + unsaturated reactions
        unsat_atm_keys = _unsaturated_atom_keys(xgra)
        # print('unsat test:', tra[0][0], unsat_atm_keys)
        tra_list = list(tra[0])
        for key in unsat_atm_keys:
            if key in tra_list[0]:
                pass
                # tra = None
                # commented out the tra = None and added pass to run CO + HO2

    return tra, idxs


def _substitution(atmsa, neighsa, bndsa, atmsb, xgr1, xgr2):
    ret_tra = None
    for atmi in atmsa:
        if _is_heavy(atmi, atmsa):
            i_neighs = neighsa[atmi]
            for atmj in i_neighs:
                bnd_break_key_ij = _get_bnd_key(atmi, atmj, bndsa)
                new_xgr = remove_bonds(xgr1, [bnd_break_key_ij])
                for atmk in atmsb:
                    if atmk not in (atmi, atmj):
                        bnd_form_key_ik = frozenset({atmi, atmk})
                        newnew_xgr = add_bonds(new_xgr, [bnd_form_key_ik])
                        atm_key_dct = _full_isomorphism(newnew_xgr, xgr2)
                        if atm_key_dct:
                            tra = [[bnd_form_key_ik], [bnd_break_key_ij]]
                            ret_tra = tra
    return ret_tra


def insertion(xgr1, xgr2):
    """ find an insertion transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    idxs = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)
    if len(xgrs1) == 2 and len(xgrs2) == 1:
        xgra, xgrb = xgrs1
        atmsa = atoms(xgra)
        atmsb = atoms(xgrb)
        neighsa = atom_neighbor_keys(xgra)
        neighsb = atom_neighbor_keys(xgrb)
        bndsa = bond_keys(xgra)
        bndsb = bond_keys(xgrb)
        tra = _insertion(atmsa, neighsa, bndsa, atmsb, neighsb, xgr1, xgr2)
        idxs = [0, 1]
        if not tra:
            tra = _insertion(atmsb, neighsb, bndsb, atmsa, neighsa, xgr1, xgr2)
            if tra:
                idxs = [1, 0]
            else:
                idxs = None
    elif len(xgrs1) == 1 and len(xgrs2) == 1:
        xgra = xgr1
        idxs = [0]
        atmsa = atoms(xgra)
        neighsa = atom_neighbor_keys(xgra)
        bndsa = bond_keys(xgra)
        tra = _insertion(atmsa, neighsa, bndsa, atmsa, neighsa, xgr1, xgr2)
    return tra, idxs


def _insertion(atmsa, neighsa, bndsa, atmsb, neighsb, xgr1, xgr2):
    """Do the insertion for an order of reactants
    """
    ret_tra = None
    for i in atmsa:
        if _is_heavy(i, atmsa):
            i_neighs = neighsa[i]
            for j in i_neighs:
                bnd_break_key_ij = _get_bnd_key(i, j, bndsa)
                new_xgr = remove_bonds(xgr1, [bnd_break_key_ij])
                for k in atmsb:
                    if _is_heavy(k, atmsb) and (
                            k != i and k != j and i not in neighsb[k] and
                            j not in neighsb[k]):
                        bnd_form_key_ik = {i, k}
                        bnd_form_key_jk = {j, k}
                        newnew_xgr = add_bonds(
                            new_xgr, [bnd_form_key_ik, bnd_form_key_jk])
                        atm_key_dct = _full_isomorphism(newnew_xgr, xgr2)
                        if atm_key_dct:
                            tra = [
                                [bnd_form_key_ik, bnd_form_key_jk],
                                [bnd_break_key_ij]]
                            ret_tra = tra
    return ret_tra


def _is_heavy(idx, atms):
    """Determine if atom index points to a heavy atom
    """
    heavy = True
    if atms[idx][0] == 'H':
        heavy = False
    return heavy


def _get_bnd_key(idx1, idx2, _bond_keys):
    """get bond key for two atoms
    """
    ret = None
    for bnd in _bond_keys:
        if idx1 in bnd and idx2 in bnd:
            ret = bnd
            break
    return ret


def hydrogen_abstraction(xgr1, xgr2):
    """ find an addition transformation
    """
    assert xgr1 == _explicit(xgr1) and xgr2 == _explicit(xgr2)

    tra = None
    xgrs1 = _connected_components(xgr1)
    xgrs2 = _connected_components(xgr2)

    ret = formula.reac.argsort_hydrogen_abstraction(
        list(map(convert.graph.formula, xgrs1)),
        list(map(convert.graph.formula, xgrs2)))
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
                'hydrogen abstraction',
                frm_bnd_keys=formed_bond_keys(q2_tra),
                brk_bnd_keys=broken_bond_keys(q1_tra))

    return tra


def _is_bond_key(obj):
    return (isinstance(obj, frozenset) and len(obj) == 2 and
            all(map(_is_atom_key, obj)))


def _is_atom_key(obj):
    return isinstance(obj, numbers.Integral)
