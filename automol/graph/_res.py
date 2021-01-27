""" resonance graph library
"""

import itertools
import functools
import numpy
from automol.util import dict_
from automol.graph._graph import atom_keys
from automol.graph._graph import bond_keys
from automol.graph._graph import remove_bonds
from automol.graph._graph import remove_atoms
from automol.graph._graph import bond_orders
from automol.graph._graph import set_bond_orders
from automol.graph._graph import without_bond_orders
from automol.graph._graph import atoms_neighbor_atom_keys
from automol.graph._graph import atom_unsaturated_valences
from automol.graph._graph import atom_bond_valences
from automol.graph._graph import atom_lone_pair_counts
from automol.graph._graph import maximum_spin_multiplicity
from automol.graph._graph import explicit
from automol.graph._graph import implicit
from automol.graph._graph import atom_explicit_hydrogen_valences
from automol.graph._graph import atoms
from automol.graph._graph import atom_groups
from automol.graph._graph import full_isomorphism


# atom properties
def atom_hybridizations(rgr):
    """ atom hybridizations, by atom
    """
    atm_keys = list(atom_keys(rgr))
    atm_unsat_vlc_dct = atom_unsaturated_valences(rgr, bond_order=True)
    atm_bnd_vlc_dct = atom_bond_valences(rgr, bond_order=False)     # note!!
    atm_unsat_vlcs = numpy.array(
        dict_.values_by_key(atm_unsat_vlc_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(
        dict_.values_by_key(atom_lone_pair_counts(rgr), atm_keys))
    atm_hybs = atm_unsat_vlcs + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(
        dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def resonance_dominant_atom_hybridizations(rgr):
    """ resonance-dominant atom hybridizations, by atom
    """
    atm_keys = list(atom_keys(rgr))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def linear_atom_keys(rgr):
    """ atoms forming linear bonds, based on their hybridization
    """
    atm_hyb_dct = resonance_dominant_atom_hybridizations(implicit(rgr))
    atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
    return atm_keys


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

    atm_ngb_keys_dct = atoms_neighbor_atom_keys(rgr)

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


def nonresonant_radical_atom_keys(rgr):
    """ keys for radical atoms that are not in resonance
    """
    atm_keys = list(atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_rad_vlcs = [min(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def resonance_dominant_radical_atom_keys(rgr):
    """ resonance-dominant radical atom keys

    (keys of resonance-dominant radical sites)
    """
    atm_keys = list(atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def sing_res_dom_radical_atom_keys(rgr):
    """ resonance-dominant radical atom keys,for one resonance
    """
    atm_keys = list(atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    first_atm_rad_val = [atm_rad_vlcs_by_res[0]]
    atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*first_atm_rad_val)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def radical_groups(gra):
    """ returns a list of lists of groups attached each radical
    """

    groups = []
    rads = sing_res_dom_radical_atom_keys(gra)
    for rad in rads:
        groups.append(atom_groups(gra, rad))
    return groups


def radical_group_dct(gra):
    """ return a dictionary of lists of groups attached each radical
    """
    groups = {}
    rads = list(sing_res_dom_radical_atom_keys(gra))
    atms = atoms(gra)
    for rad in rads:
        key = atms[rad][0]
        if key in groups:
            groups[atms[rad][0]] += (atom_groups(gra, rad),)
        else:
            groups[atms[rad][0]] = atom_groups(gra, rad)

    return groups


def radical_dissociation_prods(gra, pgra1):
    """ given a dissociation product, determine the other product
    """

    pgra2 = None
    rads = sing_res_dom_radical_atom_keys(gra)
    adj_atms = atoms_neighbor_atom_keys(gra)
    # adj_idxs = tuple(adj_atms[rad] for rad in rads)
    for rad in rads:
        for adj in adj_atms[rad]:
            for group in atom_groups(gra, adj):
                if full_isomorphism(explicit(group), explicit(pgra1)):
                    pgra2 = remove_atoms(gra, atom_keys(group))
                    # pgra2 = remove_bonds(pgra2, bond_keys(group))
                    if bond_keys(group) in pgra2:
                        pgra2 = remove_bonds(pgra2, bond_keys(group))
    return (pgra1, pgra2)


# bond properties
def resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    bnd_ords_lst = list(map(frozenset, zip(*bnd_ords_by_res)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


def one_resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    first_bnd_ords = [bnd_ords_by_res[0]]
    bnd_ords_lst = list(map(frozenset, zip(*first_bnd_ords)))
    bnd_dom_res_ords_dct = dict(zip(bnd_keys, bnd_ords_lst))
    return bnd_dom_res_ords_dct


def resonance_avg_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    bnd_keys = list(bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(bond_orders(dom_rgr), bnd_keys)
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
    mult_min = min(map(maximum_spin_multiplicity, rgrs))
    dom_rgrs = tuple(
        rgr for rgr in rgrs if maximum_spin_multiplicity(rgr) == mult_min)
    return dom_rgrs


def resonances(rgr):
    """ all resonance graphs with this connectivity
    """
    return subresonances(without_bond_orders(rgr))


def subresonances(rgr):
    """ this connected graph and its lower-spin (more pi-bonded) resonances
    """
    # get the bond capacities (room for increasing bond order), filtering out
    # the negative ones to avoid complications with hypervalent atoms in TSs
    bnd_cap_dct = dict_.by_value(_bond_capacities(rgr), lambda x: x > 0)

    ret_rgrs = []
    if bnd_cap_dct:
        bnd_keys, bnd_caps = zip(*bnd_cap_dct.items())
        atm_keys = list(functools.reduce(frozenset.union, bnd_keys))

        # Loop over all possible combinations of bond order increments (amounts
        # by which to increase the bond order), filtering out combinations that
        # exceed the valences of the atoms involved.
        # (Note that we are only testing the bonds with available pi electrons,
        # so this is compatible with having hypervalent atoms elsewhere in the
        # molecule)
        bnd_ord_inc_ranges = [range(bnd_cap+1) for bnd_cap in bnd_caps]
        for bnd_ord_incs in itertools.product(*bnd_ord_inc_ranges):
            bnd_ord_inc_dct = dict(zip(bnd_keys, bnd_ord_incs))
            ret_rgr = _add_pi_bonds(rgr, bnd_ord_inc_dct)

            atm_unsat_vlcs = dict_.values_by_key(
                atom_unsaturated_valences(ret_rgr), atm_keys)

            if not any(atm_unsat_vlc < 0 for atm_unsat_vlc in atm_unsat_vlcs):
                ret_rgrs.append(ret_rgr)

    if not ret_rgrs:
        ret_rgrs = (rgr,)
    else:
        ret_rgrs = tuple(ret_rgrs)

    return ret_rgrs


def _bond_capacities(rgr):
    """ the number of electron pairs available for further pi-bonding, by bond
    """
    atm_unsat_vlc_dct = atom_unsaturated_valences(rgr)

    def _pi_capacities(bnd_key):
        return min(map(atm_unsat_vlc_dct.__getitem__, bnd_key))

    bnd_keys = list(bond_keys(rgr))
    bnd_caps = tuple(map(_pi_capacities, bnd_keys))
    bnd_cap_dct = dict(zip(bnd_keys, bnd_caps))
    return bnd_cap_dct


def _add_pi_bonds(rgr, bnd_ord_inc_dct):
    """ add pi bonds to this graph
    """
    bnd_keys = bond_keys(rgr)
    assert set(bnd_ord_inc_dct.keys()) <= bnd_keys

    bnd_keys = list(bnd_keys)
    bnd_ords = dict_.values_by_key(bond_orders(rgr), bnd_keys)
    bnd_ord_incs = dict_.values_by_key(bnd_ord_inc_dct, bnd_keys, fill_val=0)
    new_bnd_ords = numpy.add(bnd_ords, bnd_ord_incs)
    bnd_ord_dct = dict(zip(bnd_keys, new_bnd_ords))
    rgr = set_bond_orders(rgr, bnd_ord_dct)
    return rgr


# other utilities
def rotational_bond_keys(gra, with_h_rotors=True):
    """ determine rotational bonds in this molecular graph
    """
    gra = explicit(gra)
    atm_bnd_vlc_dct = atom_bond_valences(gra, bond_order=False)
    atm_exp_hyd_vlc_dct = atom_explicit_hydrogen_valences(gra)
    res_dom_bnd_ords_dct = resonance_dominant_bond_orders(gra)

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
