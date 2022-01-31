""" graph functions that depend on resonance structure

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import functools
import numpy
from automol.util import dict_
from automol.graph.base._algo import atom_groups
from automol.graph.base._algo import isomorphism
from automol.graph.base._algo import connected_components
from automol.graph.base._core import atoms
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bond_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import set_bond_orders
from automol.graph.base._core import atom_hybridizations
from automol.graph.base._core import atom_unsaturated_valences
from automol.graph.base._core import remove_atoms
from automol.graph.base._core import remove_bonds
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_dummy_bonds
from automol.graph.base._core import without_fractional_bonds
from automol.graph.base._core import maximum_spin_multiplicity
from automol.graph.base._core import implicit
from automol.graph.base._core import subgraph
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import dummy_atoms_neighbor_atom_key


# # core functions
def dominant_resonance(rgr):
    """ *a* dominant (minimum spin/maximum pi) resonance graph
    """
    return next(iter(dominant_resonances(rgr)))


def dominant_resonances(rgr):
    """ all dominant (minimum spin/maximum pi) resonance graphs
    """
    rgr = without_fractional_bonds(rgr)
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
    rgr = without_fractional_bonds(rgr)
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

            max_bnd_ord = max(bond_orders(ret_rgr).values())

            atm_unsat_vlcs = dict_.values_by_key(
                atom_unsaturated_valences(ret_rgr), atm_keys)

            if not any(atm_unsat_vlc < 0 for atm_unsat_vlc in atm_unsat_vlcs):
                if max_bnd_ord < 4:
                    ret_rgrs.append(ret_rgr)

    if not ret_rgrs:
        ret_rgrs = (rgr,)
    else:
        ret_rgrs = tuple(ret_rgrs)

    return ret_rgrs


def resonance_dominant_bond_orders(rgr):
    """ resonance-dominant bond orders, by bond
    """
    rgr = without_fractional_bonds(rgr)
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
    rgr = without_fractional_bonds(rgr)
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
    rgr = without_fractional_bonds(rgr)
    bnd_keys = list(bond_keys(rgr))
    bnd_ords_by_res = [
        dict_.values_by_key(bond_orders(dom_rgr), bnd_keys)
        for dom_rgr in dominant_resonances(rgr)]
    nres = len(bnd_ords_by_res)
    bnd_ords_lst = zip(*bnd_ords_by_res)
    avg_bnd_ord_lst = [sum(bnd_ords)/nres for bnd_ords in bnd_ords_lst]
    avg_bnd_ord_dct = dict(zip(bnd_keys, avg_bnd_ord_lst))
    return avg_bnd_ord_dct


# # derived properties
def linear_atom_keys(rgr, dummy=True):
    """ atoms forming linear bonds, based on their hybridization

    :param rgr: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    rgr = without_fractional_bonds(rgr)
    atm_hyb_dct = resonance_dominant_atom_hybridizations(implicit(rgr))
    lin_atm_keys = set(dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1))

    if dummy:
        dum_ngb_key_dct = dummy_atoms_neighbor_atom_key(rgr)
        lin_atm_keys |= set(dum_ngb_key_dct.values())

    lin_atm_keys = tuple(sorted(lin_atm_keys))
    return lin_atm_keys


def linear_segments_atom_keys(gra, lin_keys=None):
    """ atom keys for linear segments in the graph
    """
    ngb_keys_dct = atoms_neighbor_atom_keys(without_dummy_atoms(gra))

    lin_keys = (linear_atom_keys(gra, dummy=True)
                if lin_keys is None else lin_keys)

    lin_keys = [k for k in lin_keys if len(ngb_keys_dct[k]) <= 2]

    lin_segs = connected_components(subgraph(gra, lin_keys))

    lin_keys_lst = []
    for lin_seg in lin_segs:
        lin_seg_keys = atom_keys(lin_seg)
        if len(lin_seg_keys) == 1:
            key, = lin_seg_keys
            lin_keys_lst.append([key])
        else:
            end_key1, end_key2 = sorted([
                key for key, ngb_keys in
                atoms_neighbor_atom_keys(lin_seg).items()
                if len(ngb_keys) == 1])
            ngb_keys_dct = atoms_neighbor_atom_keys(lin_seg)

            key = None
            keys = [end_key1]
            while key != end_key2:
                key, = ngb_keys_dct[keys[-1]] - set(keys)
                keys.append(key)
            lin_keys_lst.append(keys)

    lin_keys_lst = tuple(map(tuple, lin_keys_lst))
    return lin_keys_lst


def radical_atom_keys(gra, single_res=False, min_valence=1.):
    """ Radical atom keys for this molecular graph

    Radical atoms are based on the lowest-spin resonance structures for this
    graph. If the `single_res` flag is set, a single low-spin resonance
    structure will be chosen when there are multiple such structures.

    This function should eventually replace both
    `resonance_dominant_radical_atom_keys` and
    `sing_res_dom_radical_atom_keys`.

    If you wish to identify radical atom keys based on the bond orders in
    `gra`, this can be done by using the `radical_atom_keys_from_resonance`
    function.

    :param gra: the molecular graph
    :param single_res: only include radical keys for a single (arbitrary)
        resonance structure, or include all atoms that are radicals in any of
        the low-spin resonance structures?
    :type single_res: bool
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    gra = without_fractional_bonds(gra)
    atm_keys = list(atom_keys(gra))

    if single_res:
        atm_rad_vlcs = dict_.values_by_key(
            atom_unsaturated_valences(dominant_resonance(gra)), atm_keys)
    else:
        atm_rad_vlcs_by_res = [
            dict_.values_by_key(atom_unsaturated_valences(dom_gra), atm_keys)
            for dom_gra in dominant_resonances(gra)]
        atm_rad_vlcs = [
            max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]

    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs)
                             if atm_rad_vlc >= min_valence)
    return atm_rad_keys


def radical_atom_keys_from_resonance(rgr, min_valence=1.):
    """ Radical atom keys for a resonance molecular graph

    Assumes the graph has already been assinged to a resonance structure.

    :param rgr: a resonance-structure molecular graph
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    rgr = without_fractional_bonds(rgr)
    atm_keys = list(atom_keys(rgr))

    atm_rad_vlcs = dict_.values_by_key(
        atom_unsaturated_valences(rgr), atm_keys)

    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs)
                             if atm_rad_vlc >= min_valence)
    return atm_rad_keys


def has_separated_radical_sites(gra):
    """ does this radical have two or more separated radical sites?

    The identification is performed based on one of its lowest-spin resonance
    structures. It shouldn't matter which of the low-spin resonance structures
    is used -- if one of them has separated radical sites, they all should.

    This identifies polyradical molecules, but excludes things like carbenes.

    :param gra: the graph
    :returns: True if it has, False if not
    :rtype: bool
    """
    rad_atm_keys = radical_atom_keys(gra, single_res=True)
    return len(rad_atm_keys) > 1


def nonresonant_radical_atom_keys(rgr):
    """ keys for radical atoms that are not in resonance
    """
    rgr = without_fractional_bonds(rgr)
    atm_keys = list(atom_keys(rgr))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unsaturated_valences(dom_rgr), atm_keys)
        for dom_rgr in dominant_resonances(rgr)]
    atm_rad_vlcs = [min(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(atm_key for atm_key, atm_rad_vlc
                             in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc)
    return atm_rad_keys


def vinyl_radical_atom_keys(rgr):
    """ Vinyl radical atom keys for this molecular graph

    :param rgr: the molecular graph
    :returns: the vinyl radical atom keys
    :rtype: frozenset[int]
    """
    rgr = without_fractional_bonds(rgr)
    atm_rad_keys = nonresonant_radical_atom_keys(rgr)
    bnd_ords_dct = resonance_dominant_bond_orders(rgr)
    atm_bnd_keys_dct = atoms_bond_keys(rgr)
    atm_vin_keys = []
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 2 in bnd_ords_dct[bnd_key]:
                atm_vin_keys.append(atm_key)
                break
    atm_vin_keys = frozenset(atm_vin_keys)
    return atm_vin_keys


def sigma_radical_atom_keys(rgr):
    """ keys for sigma radical atoms

    :param rgr: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    rgr = without_fractional_bonds(rgr)
    atm_rad_keys = nonresonant_radical_atom_keys(rgr)
    bnd_ords_dct = resonance_dominant_bond_orders(rgr)
    atm_bnd_keys_dct = atoms_bond_keys(rgr)
    atm_sig_keys = []
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 3 in bnd_ords_dct[bnd_key]:
                atm_sig_keys.append(atm_key)
                break
    atm_sig_keys = frozenset(atm_sig_keys)
    return atm_sig_keys


def resonance_dominant_radical_atom_keys(rgr):
    """ resonance-dominant radical atom keys

    TODO: DEPRECATE

    (keys of resonance-dominant radical sites)
    """
    rgr = without_fractional_bonds(rgr)
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

    TODO: DEPRECATE
    """
    rgr = without_fractional_bonds(rgr)
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
    gra = without_fractional_bonds(gra)

    groups = []
    rads = sing_res_dom_radical_atom_keys(gra)
    for rad in rads:
        groups.append(atom_groups(gra, rad))
    return groups


def radical_group_dct(gra):
    """ return a dictionary of lists of groups attached each radical
    """
    gra = without_fractional_bonds(gra)

    groups = {}
    rads = list(sing_res_dom_radical_atom_keys(gra))
    atms = atoms(gra)
    for rad in rads:
        key = atms[rad][0]
        if key in groups:
            groups[atms[rad][0]] += atom_groups(gra, rad)
        else:
            groups[atms[rad][0]] = atom_groups(gra, rad)

    return groups


def radical_dissociation_prods(gra, pgra1):
    """ given a dissociation product, determine the other product

        stereo has been turned on atom_groups call to see if we
        can return the graphs with stereo, previously it was off
    """
    gra = without_fractional_bonds(gra)

    pgra2 = None
    rads = sing_res_dom_radical_atom_keys(gra)
    adj_atms = atoms_neighbor_atom_keys(gra)
    # adj_idxs = tuple(adj_atms[rad] for rad in rads)
    for rad in rads:
        for i, adj in enumerate(adj_atms[rad]):
            for j, group in enumerate(atom_groups(gra, adj, stereo=False)):  
                # print('\n', i, j)
                # print(group)
                # print(pgra1)
                if isomorphism(group, pgra1, backbone_only=True):
                    pgra2 = remove_atoms(gra, atom_keys(group))
                    # pgra2 = remove_bonds(pgra2, bond_keys(group))
                    if bond_keys(group) in pgra2:
                        pgra2 = remove_bonds(pgra2, bond_keys(group))
    return (pgra1, pgra2)


def sp2_bond_keys(gra):
    """ determine the sp2 bonds in this graph
    """
    gra = without_bond_orders(gra)
    bnd_keys = dict_.keys_by_value(
        resonance_dominant_bond_orders(gra), lambda x: 2 in x)

    # make sure both ends are sp^2 (excludes cumulenes)
    atm_hyb_dct = resonance_dominant_atom_hybridizations(gra)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)
    bnd_keys = frozenset({bnd_key for bnd_key in bnd_keys
                          if bnd_key <= sp2_atm_keys})
    return bnd_keys


def resonance_dominant_atom_hybridizations(rgr):
    """ resonance-dominant atom hybridizations, by atom
    """
    rgr = without_fractional_bonds(rgr)
    atm_keys = list(atom_keys(rgr))
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
    rgr = without_fractional_bonds(rgr)
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
    rgr = without_fractional_bonds(rgr)
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


# # helpers
def _bond_capacities(rgr):
    """ the number of electron pairs available for further pi-bonding, by bond
    """
    rgr = without_dummy_bonds(rgr)
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
