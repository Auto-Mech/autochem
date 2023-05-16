""" functions for working with rotational bonds and groups

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
from automol.util import dict_
from automol.graph.base._kekule import linear_segments_atom_keys
from automol.graph.base._kekule import kekules_bond_orders_collated
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import explicit
from automol.graph.base._core import implicit
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atom_neighbor_atom_key
from automol.graph.base._algo import branch_atom_keys
from automol.graph.base._algo import rings_bond_keys


def rotational_bond_keys(
        gra, lin_keys=None, with_h_rotors=True, with_chx_rotors=True):
    """ get all rotational bonds for a graph

    :param gra: the graph
    :param lin_keys: keys to linear atoms in the graph
    """
    gra = explicit(gra)
    sym_dct = atom_symbols(gra)
    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    bnd_ord_dct = kekules_bond_orders_collated(gra)
    rng_bnd_keys = list(itertools.chain(*rings_bond_keys(gra)))

    def _is_rotational_bond(bnd_key):
        ngb_keys_lst = [ngb_keys_dct[k] - bnd_key for k in bnd_key]

        is_single = max(bnd_ord_dct[bnd_key]) <= 1
        has_neighbors = all(ngb_keys_lst)
        not_in_ring = bnd_key not in rng_bnd_keys

        is_h_rotor = any(set(map(sym_dct.__getitem__, ks)) == {'H'}
                         for ks in ngb_keys_lst)
        is_chx_rotor = is_h_rotor and any(
            sym_dct[k] == 'C' for k in bnd_key)

        return is_single and has_neighbors and not_in_ring and (
            not is_h_rotor or with_h_rotors) and (
            not is_chx_rotor or with_chx_rotors)

    rot_bnd_keys = frozenset(filter(_is_rotational_bond, bond_keys(gra)))
    lin_keys_lst = linear_segments_atom_keys(gra, lin_keys=lin_keys)
    dum_keys = tuple(atom_keys(gra, symb='X'))
    for keys in lin_keys_lst:
        bnd_keys = sorted((k for k in rot_bnd_keys if k & set(keys)),
                          key=sorted)

        # Check whether there are neighboring atoms on either side of the
        # linear segment
        excl_keys = set(keys) | set(dum_keys)

        end_key1 = atom_neighbor_atom_key(
            gra, keys[0], excl_atm_keys=excl_keys)

        excl_keys |= {end_key1}

        end_key2 = atom_neighbor_atom_key(
            gra, keys[-1], excl_atm_keys=excl_keys)

        if end_key1 is None or end_key2 is None:
            has_neighbors = False
        else:
            end_keys = {end_key1, end_key2}
            ngb_keys_lst = [ngb_keys_dct[k] - excl_keys for k in end_keys]
            has_neighbors = all(ngb_keys_lst)

        if not has_neighbors:
            rot_bnd_keys -= set(bnd_keys)
        else:
            rot_bnd_keys -= set(bnd_keys[:-1])

    return rot_bnd_keys


def rotational_groups(gra, key1, key2, dummy=False):
    """ get the rotational groups for a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """

    if not dummy:
        gra = without_dummy_atoms(gra)

    bnd_key = frozenset({key1, key2})
    assert bnd_key in bond_keys(gra)
    grp1 = branch_atom_keys(gra, key2, bnd_key) - {key1}
    grp2 = branch_atom_keys(gra, key1, bnd_key) - {key2}
    grp1 = tuple(sorted(grp1))
    grp2 = tuple(sorted(grp2))
    return grp1, grp2


def rotational_symmetry_number(gra, key1, key2, lin_keys=None):
    """ get the rotational symmetry number along a given rotational axis

    :param gra: the graph
    :param key1: the first atom key
    :param key2: the second atom key
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))
    imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(implicit(gra))

    axis_keys = {key1, key2}
    # If the keys are part of a linear chain, use the ends of that for the
    # symmetry number calculation
    lin_keys_lst = linear_segments_atom_keys(gra, lin_keys=lin_keys)
    for keys in lin_keys_lst:
        if key1 in keys or key2 in keys:
            if len(keys) == 1:
                key1, key2 = sorted(ngb_keys_dct[keys[0]])
            else:
                key1, = ngb_keys_dct[keys[0]] - {keys[1]}
                key2, = ngb_keys_dct[keys[-1]] - {keys[-2]}
                axis_keys |= set(keys)
                break

    sym_num = 1
    for key in (key1, key2):
        if key in imp_hyd_vlc_dct:
            ngb_keys = ngb_keys_dct[key] - axis_keys
            if len(ngb_keys) == imp_hyd_vlc_dct[key] == 3:
                sym_num = 3
                break
    return sym_num


def bond_symmetry_numbers(gra, frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers, by bond

    TODO: DEPRECATE -- I think this function can be replaced with
    rotational_symmetry_number(). Passing in formed and broken keys is
    unnecessary if one passes in a TS graph, which is stored in the reaction
    object.

    the (approximate) symmetry number of the torsional potential for this bond,
    based on the hydrogen counts for each atom
    It is reduced to 1 if one of the H atoms in the torsional bond is a
    neighbor to the special bonding atom (the atom that is being transferred)
    """
    imp_gra = implicit(gra)
    atm_imp_hyd_vlc_dct = atom_implicit_hydrogen_valences(imp_gra)

    bnd_keys = bond_keys(imp_gra)

    tfr_atm = None
    if frm_bnd_key and brk_bnd_key:
        for atm_f in list(frm_bnd_key):
            for atm_b in list(brk_bnd_key):
                if atm_f == atm_b:
                    tfr_atm = atm_f

        if tfr_atm:
            neighbor_dct = atoms_neighbor_atom_keys(gra)
            nei_tfr = neighbor_dct[tfr_atm]

            atms = gra[0]
            all_hyds = []
            for atm in atms:
                if atms[atm][0] == 'H':
                    all_hyds.append(atm)
        else:
            nei_tfr = {}

    bnd_symb_num_dct = {}
    bnd_symb_nums = []
    for bnd_key in bnd_keys:
        bnd_sym = 1
        vlc = max(map(atm_imp_hyd_vlc_dct.__getitem__, bnd_key))
        if vlc == 3:
            bnd_sym = 3
            if tfr_atm:
                for atm in nei_tfr:
                    nei_s = neighbor_dct[atm]
                    h_nei = 0
                    for nei in nei_s:
                        if nei in all_hyds:
                            h_nei += 1
                    if h_nei == 3:
                        bnd_sym = 1
        bnd_symb_nums.append(bnd_sym)

    bnd_symb_num_dct = dict(zip(bnd_keys, bnd_symb_nums))

    # fill in the rest of the bonds for completeness
    bnd_symb_num_dct = dict_.by_key(
        bnd_symb_num_dct, bond_keys(gra), fill_val=1)

    return bnd_symb_num_dct
