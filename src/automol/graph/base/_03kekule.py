"""graph functions associated with kekule (resonance) structures

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import numbers
from typing import Any, Dict, List, Optional, Tuple

import numpy

from ...util import dict_
from ...util.ring import edges
from ._00core import (
    AtomKey,
    AtomKeys,
    BondKey,
    atom_bond_counts,
    atom_implicit_hydrogens,
    atom_keys,
    atom_lone_pairs,
    atom_unpaired_electrons,
    atoms,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    bond_keys,
    bond_neighbor_atom_keys,
    bond_orders,
    bond_stereo_keys,
    bond_unpaired_electrons,
    bonds_neighbor_bond_keys,
    dummy_source_dict,
    has_atom_stereo,
    implicit,
    is_ts_graph,
    local_stereo_priorities,
    set_bond_orders,
    subgraph,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reagents_graphs_without_stereo,
    ts_transferring_atoms,
    vinyl_radical_bond_candidates,
    without_dummy_atoms,
    without_pi_bonds,
)
from ._02algo import (
    branches,
    connected_components,
    connected_components_atom_keys,
    rings_atom_keys,
    shortest_path_between_atoms,
)


# # core functions
def kekule(gra, max_stereo_overlap=True):
    """One low-spin kekule graph, ignoring current bond orders

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param max_stereo_overlap: optionally, request as many stereo bonds as
        possible to have a bond order of 2
    :type max_stereo_overlap: bool
    :returns: a kekule graph
    """
    ste_bkeys = bond_stereo_keys(gra)
    nbkeys_dct = bonds_neighbor_bond_keys(gra, group=False)

    def _count_stereo_double_bonds(gra):
        good_bkeys = good_stereo_bond_keys_from_kekule(
            gra, ste_bkeys=ste_bkeys, nbkeys_dct=nbkeys_dct
        )
        return len(good_bkeys)

    def _count_ring_double_bonds(gra):
        len_rng = 0
        # extract rings and count alternate double bonds
        # maximize the total N of alternate double bonds in the molecule
        rngs_atm_keys = rings_atom_keys(gra)
        for rng_atm_keys in rngs_atm_keys:
            rng_bnd_keys = frozenset(map(frozenset, edges(rng_atm_keys)))
            good_bkeys = good_stereo_bond_keys_from_kekule(
                gra, ste_bkeys=rng_bnd_keys, nbkeys_dct=nbkeys_dct
            )
            len_rng += len(good_bkeys)
        return len_rng

    gras = kekules(gra)

    if not max_stereo_overlap:
        gra = next(iter(gras))
    else:
        # prioritize aromaticity over stereo
        gra = max(
            gras,
            key=lambda x: (_count_ring_double_bonds(x), _count_stereo_double_bonds(x)),
        )

    return gra


def kekules(gra):
    """All possible low-spin kekule graphs, ignoring current bond orders

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: all possible low-spin kekule graphs
    """
    orig_gra = gra
    bnd_ord_dcts = kekules_bond_orders(gra)
    gras = []
    for bnd_ord_dct in bnd_ord_dcts:
        gra = set_bond_orders(orig_gra, bnd_ord_dct)
        gras.append(gra)
    return tuple(gras)


def kekule_bond_orders(gra, max_stereo_overlap=True):
    """Bond orders for one low-spin kekule graph

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param max_stereo_overlap: optionally, request as many stereo bonds as
        possible to have a bond order of 2
    :type max_stereo_overlap: bool
    :returns: a bond order dictionary
    """
    return bond_orders(kekule(gra, max_stereo_overlap=max_stereo_overlap))


def kekules_bond_orders(gra):
    """Bond orders for all possible low-spin kekule graphs

    Low-spin kekule graphs have double and triple bonds assigned to
    minimize the number of unpaired electrons.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    gra = without_pi_bonds(gra)
    bord_dct0 = bond_orders(gra)

    gra = implicit(gra)

    # identify all of the independent pi systems and assign kekules to each
    pi_keys_lst = pi_system_atom_keys(gra)
    pi_bord_dcts_lst = [
        pi_system_kekules_bond_orders_brute_force(gra, pi_keys)
        for pi_keys in pi_keys_lst
    ]

    bnd_ord_dcts = []
    # combine the kekules from each pi system together in all possible ways
    for bord_dcts in itertools.product(*pi_bord_dcts_lst):
        bord_dct = bord_dct0.copy()
        for dct in bord_dcts:
            bord_dct.update(dct)
        bnd_ord_dcts.append(bord_dct)

    if bnd_ord_dcts:
        bnd_ord_dcts = tuple(bnd_ord_dcts)
    else:
        bnd_ord_dcts = (bord_dct0,)

    return bnd_ord_dcts


def kekules_bond_orders_collated(gra):
    """Bond orders for all possible low-spin kekule graphs, collated into a
    single dictionary

    For TS graphs, collates possible bond orders from both reactants and
    products.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    gras = ts_reagents_graphs_without_stereo(gra) if is_ts_graph(gra) else [gra]

    bnd_keys = list(bond_keys(gra))
    bnd_ords_lst = list(
        dict_.values_by_key(d, bnd_keys, fill_val=0)
        for g in gras
        for d in kekules_bond_orders(g)
    )
    bnd_ords_dct = dict(zip(bnd_keys, zip(*bnd_ords_lst)))
    return bnd_ords_dct


def kekules_bond_orders_averaged(gra):
    """Bond orders for all possible low-spin kekule graphs, averaged

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: bond orders for all possible low-spin kekule graphs
    :rtype: tuple[dict]
    """
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    avg_bnd_ord_dct = {k: sum(v) / len(v) for k, v in bnd_ords_dct.items()}
    return avg_bnd_ord_dct


# # derived properties
def ts_linear_reacting_atoms(tsg, ring: bool = False) -> dict[int, tuple[int, int]]:
    """Identify linear reacting atoms in a TS graph, along with their in-line neighbors.

    :param tsg: TS graph
    :param ring: Include atoms in rings?; default `False`
    :return: Mapping of linear reacting atom keys onto their in-line neighbors
    """
    rcts_gra = ts_reactants_graph_without_stereo(tsg)

    # Start with transferring atoms
    lin_dct = ts_transferring_atoms(tsg)

    # Add bond-forming sigma radicals
    sig_dct = sigma_radical_atom_bond_keys(rcts_gra)
    for sig_key, sig_bkey in sig_dct.items():
        frm_bkey = next((bk for bk in ts_forming_bond_keys(tsg) if sig_key in bk), None)
        if frm_bkey is not None:
            (att_key,) = frm_bkey - sig_bkey
            (sig_nkey,) = sig_bkey - frm_bkey
            lin_dct[sig_key] = tuple(sorted([att_key, sig_nkey]))

    # If requested, remove ring atoms
    if not ring:
        rng_keys = set(itertools.chain(*rings_atom_keys(tsg)))
        lin_dct = {k: v for k, v in lin_dct.items() if k not in rng_keys}

    return lin_dct


def ts_linear_reacting_atom_keys(tsg, ring: bool = False) -> list[int]:
    """Identify linear reacting atoms in a TS graph.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param ring: Include atoms in rings?; default `False`
    :return: Linear reacting atom keys
    """
    return frozenset(ts_linear_reacting_atoms(tsg, ring=ring).keys())


def linear_atom_keys(gra, dummy=True):
    """Find linear atoms, based on their hybridization.

    For TS graphs, includes atoms that are linear for *either* the reactants
    *or* the products. This both simplifies the way Reaction objects can be
    handled and anticipates cases where the TS structure is close to either
    reactants or products.

    :param gra: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    return frozenset(linear_atoms_neighbor_atom_keys(gra, dummy=dummy))


def linear_atoms_neighbor_atom_keys(gra, dummy=True):
    """Find linear atoms and their in-line neighbors, based on hybridization.

    For TS graphs, includes atoms that are linear for *either* the reactants
    *or* the products. This both simplifies the way Reaction objects can be
    handled and anticipates cases where the TS structure is close to either
    reactants or products.

    :param gra: the graph
    :param dummy: whether or not to consider atoms connected to dummy atoms as
        linear, if different from what would be predicted based on their
        hybridization
    :returns: the linear atom keys
    :rtype: tuple[int]
    """
    gra0 = without_dummy_atoms(gra)
    nkeys_dct0 = atoms_neighbor_atom_keys(gra0)
    ts_ = is_ts_graph(gra)
    gras = ts_reagents_graphs_without_stereo(gra) if ts_ else [gra]

    lin_nkeys_dct = {}
    for gra_ in gras:
        gra_ = without_dummy_atoms(gra_)
        nkeys_dct = atoms_neighbor_atom_keys(gra_)
        # To be linear, an atom must be both (a.) sp1 hybridized and (b.) have more than
        # 1 neighbor (exactly 2)
        atm_hyb_dct = atom_hybridizations(gra_)
        sp1_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
        lin_nkeys_dct.update(
            {k: nkeys_dct[k] for k in sp1_atm_keys if len(nkeys_dct[k]) > 1}
        )

    # If requested, include all keys associated with dummy atoms
    if dummy:
        dum_src_keys = list(dummy_source_dict(gra, dir_=False).values())
        lin_nkeys_dct.update(
            {k: nkeys_dct0[k] for k in dum_src_keys if k not in lin_nkeys_dct}
        )

    if ts_:
        lin_nkeys_dct.update(
            dict_.transform_values(ts_linear_reacting_atoms(gra, ring=False), frozenset)
        )

    return lin_nkeys_dct


def linear_segment_cap_keys(
    gra, lin_keys: Optional[List[int]] = None, extend: bool = False
) -> Dict[List[int], Tuple[Optional[int], Optional[int]]]:
    """Find linear segments in the graph, along with the keys of their "capping" atoms,
    that is the final in-line atom on either side.

    For sigma radicals and linear nitrogens, the capping atom is the same as the last
    atom in the segment.

    Otherwise, the capping atom is the neighbor of the last atom in the segment, which
    falls within the line of the linear segment.

    Examples:
        H3C-C#C-C#C-C#C-CH3
          ^(* * * * * *)^

        .C#C-C#C-C#C-CH3
        (* * * * * *)^
         ^

         N#C-H
        (* *)^
         ^

        Cl....CH3--X--Br
        ^    (*)      ^

         ^  = cap
        (*) = linear segment

    :param gra: A graph
    :type gra: automol graph data structure
    :param lin_keys: Specify the keys of linear atoms, instead of determining from graph
    :type lin_keys: Optional[List[int]]
    :param extend: Extend each segment, to include in-line neighbors on either side?
    :type extend: bool, optional
    :returns: A dictionary mapping linear segments onto their in-line neighbors
    :rtype: Dict[List[int], Tuple[Optional[int], Optional[int]]]
    """
    nkeys_dct = atoms_neighbor_atom_keys(gra)
    lin_nkeys_dct = linear_atoms_neighbor_atom_keys(gra, dummy=True)
    lin_keys = frozenset(lin_nkeys_dct.keys()) if lin_keys is None else lin_keys

    # 1. Get graphs for each linear segment
    segs = connected_components(subgraph(gra, lin_keys))

    # 2. Build sorted lists of keys for each linear segment graph
    keys_lst = []
    for seg in segs:
        seg_keys = atom_keys(seg)
        if len(seg_keys) == 1:
            keys_lst.append(seg_keys)
        else:
            # a. Find the segment ends
            seg_nkeys_dct = atoms_neighbor_atom_keys(seg)
            end_key1, end_key2 = sorted(
                [k for k, ns in seg_nkeys_dct.items() if len(ns) == 1]
            )
            # b. Sort the segment keys in order, end_key1-...-end_key2
            keys = shortest_path_between_atoms(seg, end_key1, end_key2)
            keys_lst.append(keys)

    # 3. If requested, extend the ends
    keys_lst = list(map(tuple, keys_lst))
    lin_seg_dct = {}
    for keys in keys_lst:
        end_key1 = keys[0]
        end_key2 = keys[-1]

        # We fill in the remaining atoms, in case we are determining caps from
        # pre-defined linear keys
        ext_nkeys_dct = {**nkeys_dct, **lin_nkeys_dct}

        # Identify in-line neighbors
        excl_keys = set(keys)
        ext_nkey1s = ext_nkeys_dct[end_key1] - excl_keys

        excl_keys |= ext_nkey1s
        ext_nkey2s = ext_nkeys_dct[end_key2] - excl_keys

        # Split up neighbor keys if this is a single-atom segment
        if end_key1 == end_key2 and len(ext_nkey1s) > 1:
            ext_nkeys = sorted(ext_nkey1s)
            ext_nkey1s = ext_nkeys[:1]
            ext_nkey2s = ext_nkeys[1:]

        # Add in-line neighbors to the extended keys list, if not None
        # Bug fix: Require len(ext_nkey1s) == 1 -> only do the extension if unambiguous
        ext_keys = keys
        if ext_nkey1s and len(ext_nkey1s) == 1:
            ext_keys = (*ext_nkey1s, *ext_keys)

        if ext_nkey2s and len(ext_nkey2s) == 1:
            ext_keys = (*ext_keys, *ext_nkey2s)

        if extend:
            keys = ext_keys

        # The cap keys are the ends of the extended keys lists
        lin_seg_dct[keys] = (ext_keys[0], ext_keys[-1])

    return lin_seg_dct


def linear_segments_atom_keys(
    gra, lin_keys: Optional[List[int]] = None, extend: bool = False
) -> List[List[int]]:
    """Atom keys for linear segments in the graph

    :param gra: A graph
    :type gra: automol graph data structure
    :param lin_keys: Specify the keys of linear atoms, instead of determining from graph
    :type lin_keys: Optional[List[int]]
    :param extend: Extend each segment, to include in-line neighbors on either side?
    :type extend: bool, optional
    :returns: A list of lists of keys for each segment
    :rtype: List[List[int]]
    """
    lin_seg_dct = linear_segment_cap_keys(gra, extend=extend, lin_keys=lin_keys)
    return tuple(lin_seg_dct.keys())


def unneeded_dummy_atom_keys(gra) -> frozenset:
    """Get the keys of dummy atoms which are not connected to linear atoms

    :param gra: A graph
    :type gra: automol graph data structure
    :returns: Keys for the unneeded dummy atoms
    :rtype: frozenset[int]
    """
    dum_lin_key_dct = dummy_source_dict(gra, dir_=False)
    lin_keys = linear_atom_keys(gra, dummy=False)
    return frozenset(dk for dk, lk in dum_lin_key_dct.items() if lk not in lin_keys)


def atom_hybridizations(gra, kgrs: Optional[List[Any]] = None):
    """resonance-dominant atom hybridizations, by atom

    :param gra: A graph
    :param kgrs: Optionally, pass in the kekule graphs, to avoid recalculating them
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    kgrs = kekules(gra) if kgrs is None else kgrs

    atm_keys = list(atom_keys(gra))
    atm_hybs_by_res = [
        dict_.values_by_key(atom_hybridizations_from_kekule(g), atm_keys) for g in kgrs
    ]
    atm_hybs = [min(hybs) for hybs in zip(*atm_hybs_by_res)]
    atm_hyb_dct = dict(zip(atm_keys, atm_hybs))
    return atm_hyb_dct


def atom_hybridizations_from_kekule(gra):
    """atom hybridizations, by atom"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    atm_unsat_dct = atom_unpaired_electrons(gra, bond_order=True)
    atm_bnd_vlc_dct = atom_bond_counts(gra, bond_order=False)  # note!!
    atm_unsats = numpy.array(dict_.values_by_key(atm_unsat_dct, atm_keys))
    atm_bnd_vlcs = numpy.array(dict_.values_by_key(atm_bnd_vlc_dct, atm_keys))
    atm_lpcs = numpy.array(dict_.values_by_key(atom_lone_pairs(gra), atm_keys))
    atm_hybs = atm_unsats + atm_bnd_vlcs + atm_lpcs - 1
    atm_hyb_dct = dict_.transform_values(dict(zip(atm_keys, atm_hybs)), int)
    return atm_hyb_dct


def bad_stereo_bond_keys_from_kekule(
    gra,
    ste_bkeys: Optional[frozenset[BondKey]] = None,
    nbkeys_dct: Optional[Dict[BondKey, frozenset[BondKey]]] = None,
) -> frozenset[BondKey]:
    """Identify stereo bonds which are *not* well-represented by a Kekule graph

    A stereo bond is well-represented by the Kekule graph if it is doubly-bonded and
    surrounded by single bonds

    :param gra: A kekule graph
    :param ste_bkeys: Stereo bond keys, to avoid recalculating
    :param nbkeys_dct: Neighboring bond keys, by bond key, to avoid recalculating
    :returns: The bond keys
    """
    ste_bkeys = bond_stereo_keys(gra) if ste_bkeys is None else ste_bkeys
    good_bkeys = good_stereo_bond_keys_from_kekule(
        gra, ste_bkeys=ste_bkeys, nbkeys_dct=nbkeys_dct
    )
    return ste_bkeys - good_bkeys


def good_stereo_bond_keys_from_kekule(
    gra,
    ste_bkeys: Optional[frozenset[BondKey]] = None,
    nbkeys_dct: Optional[Dict[BondKey, frozenset[BondKey]]] = None,
) -> frozenset[BondKey]:
    """Identify stereo bonds which are well-represented by a Kekule graph

    A stereo bond is well-represented by the Kekule graph if it is doubly-bonded and
    surrounded by single bonds

    :param gra: A kekule graph
    :param ste_bkeys: Stereo bond keys, to avoid recalculating
    :param nbkeys_dct: Neighboring bond keys, by bond key, to avoid recalculating
    :returns: The bond keys
    """
    ord_dct = bond_orders(gra)
    ste_bkeys = bond_stereo_keys(gra) if ste_bkeys is None else ste_bkeys
    nbkeys_dct = (
        bonds_neighbor_bond_keys(gra, group=False) if nbkeys_dct is None else nbkeys_dct
    )

    o1_bkeys = set(dict_.keys_by_value(ord_dct, lambda x: x == 1))
    o2_bkeys = set(dict_.keys_by_value(ord_dct, lambda x: x == 2))
    # Look for stereo bonds to be doubly-bonded and surrounded by single bonds
    good_bkeys = frozenset(
        bk for bk in ste_bkeys & o2_bkeys if nbkeys_dct[bk] <= o1_bkeys
    )
    return good_bkeys


def radical_atom_keys(gra, sing_res=False, min_valence=1.0):
    """Radical atom keys for this molecular graph

    Radical atoms are based on the lowest-spin resonance structures for this
    graph. If the `sing_res` flag is set, a single low-spin resonance
    structure will be chosen when there are multiple such structures.

    For TS graphs, this will only return radical atoms that are not involved in the
    reaction.

    If you wish to identify radical atom keys based on the bond orders already
    in `gra`, this can be done by using the `radical_atom_keys_from_kekule`
    function.

    :param gra: the molecular graph
    :param sing_res: only include radical keys for a single (arbitrary)
        resonance structure, or include all atoms that are radicals in any of
        the low-spin resonance structures?
    :type sing_res: bool
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]
    """
    # If this is not a TS graph, call the main function directly
    if not is_ts_graph(gra):
        return _radical_atom_keys(gra, sing_res=sing_res, min_valence=min_valence)

    # Otherwise, determine the non-reacting radical atoms
    rgra, pgra = ts_reagents_graphs_without_stereo(gra)
    rkeys = _radical_atom_keys(rgra)
    pkeys = _radical_atom_keys(pgra)
    return frozenset(rkeys & pkeys)


def _radical_atom_keys(gra, sing_res=False, min_valence=1.0):
    """Radical atom keys for this molecular graph

    Radical atoms are based on the lowest-spin resonance structures for this
    graph. If the `sing_res` flag is set, a single low-spin resonance
    structure will be chosen when there are multiple such structures.

    For TS graphs, this will only return radical atoms that are not involved in the
    reaction.

    If you wish to identify radical atom keys based on the bond orders already
    in `gra`, this can be done by using the `radical_atom_keys_from_kekule`
    function.

    :param gra: the molecular graph
    :param sing_res: only include radical keys for a single (arbitrary)
        resonance structure, or include all atoms that are radicals in any of
        the low-spin resonance structures?
    :type sing_res: bool
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))

    if sing_res:
        atm_rad_vlcs = dict_.values_by_key(
            atom_unpaired_electrons(kekule(gra)), atm_keys
        )
    else:
        atm_rad_vlcs_by_res = [
            dict_.values_by_key(atom_unpaired_electrons(dom_gra), atm_keys)
            for dom_gra in kekules(gra)
        ]
        atm_rad_vlcs = [max(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]

    atm_rad_keys = frozenset(
        atm_key
        for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
        if atm_rad_vlc >= min_valence
    )
    return atm_rad_keys


def radical_atom_keys_from_kekule(gra, min_valence=1.0):
    """Radical atom keys for a particular kekule graph

    Assumes the graph already has assigned bond orders

    :param gra: a resonance-structure molecular graph
    :param min_valence: optionally, specify that only sites with at least a
        certain number of radical electrons be included
    :type min_valence: int
    :returns: the radical atom keys
    :rtype: frozenset[int]

    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))

    atm_rad_vlcs = dict_.values_by_key(atom_unpaired_electrons(gra), atm_keys)

    atm_rad_keys = frozenset(
        atm_key
        for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
        if atm_rad_vlc >= min_valence
    )
    return atm_rad_keys


def nonresonant_radical_atom_keys(gra):
    """keys for radical atoms that are not in resonance"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_keys = list(atom_keys(gra))
    atm_rad_vlcs_by_res = [
        dict_.values_by_key(atom_unpaired_electrons(g), atm_keys) for g in kekules(gra)
    ]
    atm_rad_vlcs = [min(rad_vlcs) for rad_vlcs in zip(*atm_rad_vlcs_by_res)]
    atm_rad_keys = frozenset(
        atm_key for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs) if atm_rad_vlc
    )
    return atm_rad_keys


def vinyl_radical_atom_bond_keys(gra):
    """Get a mapping of vinyl radical bonds onto their radical atoms

    :param gra: A graph
    :return: A mappping of vinyl radical bond keys onto vinyl atom keys
    :rtype: Dict[int, frozenset[int]]
    """
    cand_dct = vinyl_radical_bond_candidates(gra)
    hyb_dct = atom_hybridizations(gra)
    vin_dct = {k: bk for bk, k in cand_dct.items() if hyb_dct[k] == 2}
    return vin_dct


def sigma_radical_atom_bond_keys(gra):
    """keys for sigma radical atoms

    :param gra: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    atm_rad_keys = nonresonant_radical_atom_keys(gra)
    bnd_ords_dct = kekules_bond_orders_collated(gra)
    atm_bnd_keys_dct = atoms_bond_keys(gra)
    sig_dct = {}
    for atm_key in atm_rad_keys:
        for bnd_key in atm_bnd_keys_dct[atm_key]:
            if 3 in bnd_ords_dct[bnd_key]:
                sig_dct[atm_key] = bnd_key
                break
    return sig_dct


def vinyl_radical_atom_keys(gra):
    """Vinyl radical atom keys for this molecular graph

    :param gra: the molecular graph
    :returns: the vinyl radical atom keys
    :rtype: frozenset[int]
    """
    return frozenset(vinyl_radical_atom_bond_keys(gra))


def sigma_radical_atom_keys(gra):
    """keys for sigma radical atoms

    :param gra: the molecular graph
    :returns: the sigma radical atom keys
    :rtype: frozenset[int]
    """
    return frozenset(sigma_radical_atom_bond_keys(gra))


def has_separated_radical_sites(gra):
    """does this radical have two or more separated radical sites?

    The identification is performed based on one of its lowest-spin resonance
    structures. It shouldn't matter which of the low-spin resonance structures
    is used -- if one of them has separated radical sites, they all should.

    This identifies polyradical molecules, but excludes things like carbenes.

    :param gra: the graph
    :returns: True if it has, False if not
    :rtype: bool
    """
    rad_atm_keys = radical_atom_keys(gra, sing_res=True)
    return len(rad_atm_keys) > 1


def addition_atom_keys(gra):
    """Determine the atoms where addition is likely to occur

    Radical sites, if available, otherwise pi-bonds.

    :param gra: a molecular graph
    """
    # 1. First, check for radical sites
    add_keys = radical_atom_keys(gra, sing_res=False)
    # 2. If none are found, check for other unsaturated sites (pi bonds)
    if not add_keys:
        unp_dct = atom_unpaired_electrons(gra, bond_order=False)
        add_keys = dict_.keys_by_value(unp_dct)
    return frozenset(add_keys)


def beta_scission_bond_keys(gra):
    """Determine the bonds where beta scissions can occur

    :param gra: a molecular graph
    :returns: the beta scission bond keys
    """
    kgrs = kekules(gra)
    bkeys = frozenset(itertools.chain(*map(beta_scission_bond_keys_from_kekule, kgrs)))
    return bkeys


def beta_scission_bond_keys_from_kekule(gra):
    """Determine beta scission bonds for a particular kekule graph

    Assumes the graph already has assigned bond orders

    :param gra: a resonance-structure molecular graph
    :returns: the beta scission bond keys
    """
    nkeys_dct = atoms_neighbor_atom_keys(gra)

    # Determine radical sites and their neighboring atoms
    rad_keys = radical_atom_keys_from_kekule(gra)
    rad_nkeys = set(itertools.chain(*map(nkeys_dct.get, rad_keys)))

    # Grab all single bonds
    bkeys = dict_.keys_by_value(bond_orders(gra), lambda x: x == 1)
    # Remove bonds with radicals
    bkeys = {bk for bk in bkeys if not bk & rad_keys}
    # Require bonds to be adjacent to radicals
    bkeys = {bk for bk in bkeys if bk & rad_nkeys}

    return bkeys


def resonance_bond_stereo_keys(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    res_bnd_ords_dct = kekules_bond_orders_collated(gra)

    res_bnd_keys = []
    for bnd_key in ste_bnd_keys:
        if bnd_key in res_bnd_ords_dct and 1 in res_bnd_ords_dct[bnd_key]:
            res_bnd_keys.append(bnd_key)

    res_bnd_keys = tuple(res_bnd_keys)
    return res_bnd_keys


def vinyl_bond_stereo_keys(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    vin_atm_keys = vinyl_radical_atom_keys(gra)

    vin_bnd_keys = []
    for bnd_key in ste_bnd_keys:
        if bnd_key & vin_atm_keys:
            vin_bnd_keys.append(bnd_key)
    vin_bnd_keys = tuple(vin_bnd_keys)
    return vin_bnd_keys


def has_resonance_bond_stereo(gra):
    """does this graph have stereo at a resonance bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    return bool(resonance_bond_stereo_keys(gra))


def has_vinyl_bond_stereo(gra):
    """does this graph have stereo at a vinyl bond?

    :param gra: the molecular graph
    :rtype: bool
    """
    return bool(vinyl_bond_stereo_keys(gra))


def has_nonkekule_bond_stereo(gra):
    """does this graph have stereo at a resonance bond that cannot be
    represented as a double bond in a kekule structure?

    :param gra: the molecular graph
    :rtype: bool
    """
    ste_bnd_keys = bond_stereo_keys(gra)
    bnd_ord_dcts = kekules_bond_orders(gra)

    nsingles_lst = []
    for bnd_ord_dct in bnd_ord_dcts:
        nsingles = len([k for k in ste_bnd_keys if bnd_ord_dct[k] == 1])
        nsingles_lst.append(nsingles)

    return min(nsingles_lst) > 0


def has_noninchi_stereo(gra):
    """does this graph have stereo that cannot be captured by an InChI string?

    :param gra: the molecular graph
    :rtype: bool
    """
    return (
        has_resonance_bond_stereo(gra)
        or has_vinyl_bond_stereo(gra)
        or has_atom_stereo(gra, symb="N")
    )


def radical_groups(gra):
    """returns a list of lists of groups attached each radical"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"

    groups = []
    rads = radical_atom_keys(gra, sing_res=True)
    for rad in rads:
        groups.append(branches(gra, rad))
    return groups


def radical_group_dct(gra):
    """return a dictionary of lists of groups attached each radical"""
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"

    groups = {}
    rads = list(radical_atom_keys(gra, sing_res=True))
    atms = atoms(gra)
    for rad in rads:
        key = atms[rad][0]
        if key in groups:
            groups[atms[rad][0]] += branches(gra, rad)
        else:
            groups[atms[rad][0]] = branches(gra, rad)

    return groups


def strict_rigid_planar_bond_keys(gra) -> frozenset[BondKey]:
    """Get bonds which are guaranteed to be rigid and planar

    :param gra: A graph
    :type gra: automol graph data structure
    :returns: The bond keys
    """
    kgrs = kekules(gra)
    hyb_dct = atom_hybridizations(gra, kgrs=kgrs)

    bkeys = set()
    for kgr in kgrs:
        ord_dct = bond_orders(kgr)
        sp2_keys = {k for k, h in hyb_dct.items() if h == 2}
        bkeys |= {bk for bk, o in ord_dct.items() if o == 2 and bk <= sp2_keys}
    return frozenset(bkeys)


def possible_rigid_planar_bond_keys(gra) -> frozenset[BondKey]:
    """Get bonds which are could be rigid and planar

    Includes all bonds bewteen sp2 atoms, whether rigid or not

    :param gra: A graph
    :type gra: automol graph data structure
    :returns: The bond keys
    """
    bnd_unp_dct = bond_unpaired_electrons(gra, bond_order=False)
    bkeys = dict_.keys_by_value(bnd_unp_dct, lambda x: x == 1)
    return frozenset(bkeys)


def rigid_planar_bonds(
    gra,
    min_ncount: int = 1,
    min_ring_size: int = 8,
    strict: bool = True,
    excl_keys: frozenset[int] = frozenset(),
) -> Dict[BondKey, Tuple[AtomKeys, AtomKeys]]:
    """Get a mapping of rigid, planary bond keys onto their neighbor keys

    The neighbor keys are sorted by local priority

    :param gra: A graph
    :type gra: automol graph data structure
    :param min_ncount: Minimum # neighbors on either side for inclusion, defaults to 1
        (If min_ncount = 0, this will still require at least one neighbor on one side)
    :param min_ring_size: Minimum ring size for inclusion, defaults to 8
    :param strict: Only include bonds that are guaranteed to be rigid?
    :param excl_keys: Atom keys whose bonds should be excluded (for linear atoms)
    :returns: A mapping of rigid, planar bond keys onto their neighbor keys; The pair of
        neighbor key lists is sorted by the atom key that they are neighbors to
    :rtype: Dict[BondKey, Tuple[AtomKeys, AtomKeys]]
    """
    gra = without_dummy_atoms(gra)  # remove dummy atoms
    excl_keys = frozenset(excl_keys)

    gras = ts_reagents_graphs_without_stereo(gra) if is_ts_graph(gra) else [gra]
    nhyd_dct = atom_implicit_hydrogens(gra)
    pri_dct = local_stereo_priorities(gra, with_none=True)

    # 1. Find rigid, planar bonds in the graph, along with their neighbors
    # (Initially, the neighbor keys will be stored as sets)
    rp_dct = {}
    for gra_ in gras:
        rp_bkeys = (
            strict_rigid_planar_bond_keys(gra_)
            if strict
            else possible_rigid_planar_bond_keys(gra_)
        )

        rp_bkeys = {bk for bk in rp_bkeys if not bk & excl_keys}

        for bkey in rp_bkeys:
            key1, key2 = sorted(bkey)
            # Get previous neighbor keys, if any
            nkeys1, nkeys2 = rp_dct[bkey] if bkey in rp_dct else (set(), set())
            # Get current neighbor keys
            nkeys1_, nkeys2_ = bond_neighbor_atom_keys(gra_, key1, key2)
            # Combine to give full sets of neighbor keys for TS graphs
            rp_dct[bkey] = (nkeys1 | nkeys1_, nkeys2 | nkeys2_)

    # Convert to tuples
    rp_dct = dict_.transform_values(rp_dct, lambda bnks: tuple(map(tuple, bnks)))

    # 2. Enforce the minimum neighbor count
    # (Neighbor keys will be converted to tuples, with Nones for implicit hydrogens)
    for bkey in list(rp_dct.keys()):
        # Convert to tuple
        bnkeys = rp_dct[bkey]
        # Create placeholders for implicit hydrogens
        bhkeys = [(None,) * nhyd_dct[k] for k in sorted(bkey)]
        # Append the implicit hydrogen placeholders
        bnkeys = [n + h for n, h in zip(bnkeys, bhkeys)]
        # Get the total numbers of neighbors on either side
        ncounts = list(map(len, bnkeys))
        # Remove the bond if it doesn't have enough neighbors
        # Even if min_ncount = 0 we need at least one neighbor on one side for planarity
        if any(ncount < min_ncount for ncount in ncounts) or not any(ncounts):
            rp_dct.pop(bkey)
        else:
            rp_dct[bkey] = tuple(tuple(sorted(n, key=pri_dct.get)) for n in bnkeys)

    # 3. Enforce the minimum ring size
    rp_rng_const_dct = rigid_planar_bonds_with_ring_constraints(
        gra, rp_dct, min_ring_size=min_ring_size
    )
    rp_dct = dict_.filter_by_key(rp_dct, lambda k: k not in rp_rng_const_dct)
    return rp_dct


def rigid_planar_bonds_with_ring_constraints(
    gra, rp_dct: Dict[frozenset[int], tuple[tuple, tuple]], min_ring_size: int = 8
) -> Dict[BondKey, tuple[AtomKey, AtomKey]]:
    """Given a rigid, planar bond dictionary, identify bonds which are constrained by
    small rings, along with the pairs of neighbors that are constrained to be cis

    :param gra: A graph
    :type gra: automol graph data structure
    :returns: A mapping of rigid, planar bond keys onto their neighbor keys
    :type rp_dct: Dict[frozenset[int], tuple[tuple, tuple]]
    :param min_ring_size: Minimum ring size for ignoring the constraint, defaults to 8
    :type min_ring_size: int, optional
    :return: A mapping of constrained rigid, planar bond keys onto the pair of neighbors
        which are forced by the ring to be cis
    :rtype: Dict[frozenset[int], tuple[int, int]]
    """
    # Get the pool of ring keys
    rng_keys_pool = list(map(set, rings_atom_keys(gra, ts_=True)))

    # Search for constrained rigid, planar bonds
    rp_rng_const_dct = {}
    for bkey, bnkeys in rp_dct.items():
        if not isinstance(bkey, numbers.Number):  # Allow dict to contain atoms
            nkeys1, nkeys2 = bnkeys
            # Identify sufficiently small rings containing the bond
            rng_keys_lst = [
                rks for rks in rng_keys_pool if bkey <= rks and len(rks) < min_ring_size
            ]
            # Check for constraints in smaller rings first
            for rng_keys in sorted(rng_keys_lst, key=len):
                nkey1 = next((k for k in nkeys1 if k in rng_keys), None)
                nkey2 = next((k for k in nkeys2 if k in rng_keys), None)
                if nkey1 is not None and nkey2 is not None:
                    rp_rng_const_dct[bkey] = (nkey1, nkey2)
                    break

    return rp_rng_const_dct


def rigid_planar_bond_keys(
    gra, min_ncount: int = 1, min_ring_size: int = 8
) -> frozenset[BondKey]:
    """Get the keys to bonds which are rigid and planar

    This can be used to find candidates for bond stereochemistry

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param min_ncount: Minimum # neighbors on either side for inclusion, defaults to 1
        (If min_ncount = 0, this will still require at least one neighbor on one side)
    :type min_ncount: int, optional
    :param min_ring_size: Minimum ring size for inclusion, defaults to 8
    :type min_ring_size: int, optional
    """
    return frozenset(
        rigid_planar_bonds(gra, min_ncount=min_ncount, min_ring_size=min_ring_size)
    )


def atom_centered_cumulene_keys(gra):
    """resonance dominant keys for atom-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}), cent_atm_key)
    where the first pair contains the sp2 atoms at the cumulene ends and
    `cent_atm_key` is the key of the central atom
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 1:
            cum_keys.add(
                (frozenset({cum_chain[0], cum_chain[-1]}), cum_chain[size // 2])
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


def bond_centered_cumulene_keys(gra):
    """resonance dominant keys for bond-centered cumulenes

    the bond-centered cumulenes are described by
        (frozenset({end_atm_key1, end_atm_key2}),
         frozenset({cent_atm_key1, cent_atm_key2}))
    where the first pair contains the sp2 atoms at the cumulene ends and the
    second pair is the bond key for the central bond
    """
    assert not is_ts_graph(gra), f"This doesn't work for TS graphs:\n{gra}"
    cum_chains = _cumulene_chains(gra)
    cum_keys = set()
    for cum_chain in cum_chains:
        size = len(cum_chain)
        if size % 2 == 0:
            cum_keys.add(
                (
                    frozenset({cum_chain[0], cum_chain[-1]}),
                    frozenset({cum_chain[size // 2 - 1], cum_chain[size // 2]}),
                )
            )
    cum_keys = frozenset(cum_keys)
    return cum_keys


# # helpers
def pi_system_atom_keys(gra):
    """Extract keys for each closed, connected pi-system of a molecule

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: keys for each closed, connected pi system
    """
    atm_unsat_dct = atom_unpaired_electrons(gra, bond_order=False)
    all_pi_keys = dict_.keys_by_value(atm_unsat_dct)
    pi_keys_lst = tuple(
        ks
        for ks in connected_components_atom_keys(subgraph(gra, all_pi_keys))
        if len(ks) > 1
    )
    return pi_keys_lst


def pi_system_kekules_bond_orders_brute_force(gra, pi_keys, log=False):
    """Determine kekules for a closed, connected pi-system

    A general brute-force algorithm

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pi_keys: keys of an closed, connected pi system
    :type pi_keys: frozenset[int]
    :param log: print a debugging log of the number of recursive calls?
    :type log: bool
    """
    pi_sy = subgraph(gra, pi_keys)
    atm_keys = list(atom_keys(pi_sy))
    bnd_keys = list(bond_keys(pi_sy))

    aus_dct = dict_.by_key(atom_unpaired_electrons(gra, bond_order=False), atm_keys)
    bus_dct = dict_.by_key(bond_unpaired_electrons(gra, bond_order=False), bnd_keys)

    spin_max = spin_min = sum(aus_dct.values())

    # Allows no more than 2 copies of a given key in the pool, to prohibit
    # quadruple bonding for [C][C]
    bnd_pool = list(
        itertools.chain(*(itertools.repeat(k, min(bus_dct[k], 2)) for k in bnd_keys))
    )

    bnd_ord_dcts = []
    niter = 0
    for inc in range(spin_max // 2, 0, -1):
        spin = spin_max - inc * 2

        if spin > spin_min:
            break

        for bkeys in itertools.combinations(bnd_pool, inc):
            niter += 1

            akeys = list(itertools.chain(*bkeys))
            excess = next((k for k in set(akeys) if aus_dct[k] < akeys.count(k)), None)
            if excess is None:
                spin_min = spin
                bnd_ord_dct = {k: 1 + bkeys.count(k) for k in set(bkeys)}
                if bnd_ord_dct not in bnd_ord_dcts:
                    bnd_ord_dcts.append(bnd_ord_dct)

    if log:
        print(f"niter: {niter}")

    return bnd_ord_dcts


def _cumulene_chains(gra):
    atm_hyb_dct = atom_hybridizations(gra)
    sp1_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 1)
    sp2_atm_keys = dict_.keys_by_value(atm_hyb_dct, lambda x: x == 2)

    atm_ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    def _cumulene_chain(chain):
        ret = None
        atm_key = chain[-1]
        next_atm_keys = atm_ngb_keys_dct[atm_key] - {chain[-2]}
        if next_atm_keys:
            assert len(next_atm_keys) == 1
            (next_atm_key,) = next_atm_keys
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


if __name__ == "__main__":
    # # C=C[CH2]
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=CC=C[CH2]
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({1, 2}): (1, None)})

    # # C=C=C=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 0, None), 2: ('C', 0, None),
    #         3: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None)})

    # # C=CC=CC=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 2, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None)})

    # # C=C-C(-[CH2])=C
    # GRA = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 2, None), 6: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({2, 6}): (1, None)})

    # # C1=CC=CC=C1 (benzene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)})

    # # C12=CC=C1C=C2  _  _
    # #              ||_||_||
    # GRA = ({0: ('C', 0, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({0, 5}): (1, None)})

    # # C1=CC=C2C=CC=CC2=C1 (naphthalene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({8, 3}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({0, 9}): (1, None)})

    # # C1=CC=C2C=C3C=CC=CC3=CC2=C1 (anthracene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 1, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({0, 13}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({3, 12}): (1, None), frozenset({10, 5}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC=C2C(=C1)C=CC3=CC=CC=C32 (phenanthrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 1, None), 13: ('C', 0, None)},
    #        {frozenset({3, 4}): (1, None), frozenset({4, 6}): (1, None),
    #         frozenset({2, 3}): (1, None), frozenset({11, 12}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({4, 5}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
    #         frozenset({3, 13}): (1, None), frozenset({0, 5}): (1, None),
    #         frozenset({8, 7}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({12, 13}): (1, None), frozenset({10, 11}): (1, None)})

    # # C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2 (pyrene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 1, None), 10: ('C', 1, None), 11: ('C', 1, None),
    #         12: ('C', 0, None), 13: ('C', 0, None), 14: ('C', 1, None),
    #         15: ('C', 1, None)},
    #        {frozenset({4, 6}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({14, 15}): (1, None), frozenset({10, 11}): (1, None),
    #         frozenset({3, 4}): (1, None), frozenset({3, 13}): (1, None),
    #         frozenset({0, 1}): (1, None), frozenset({6, 7}): (1, None),
    #         frozenset({8, 13}): (1, None), frozenset({12, 14}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({2, 15}): (1, None),
    #         frozenset({12, 13}): (1, None)})

    # # C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7 (coronene)
    # GRA = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
    #         3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
    #         6: ('C', 1, None), 7: ('C', 1, None), 8: ('C', 0, None),
    #         9: ('C', 0, None), 10: ('C', 0, None), 11: ('C', 0, None),
    #         12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 1, None),
    #         15: ('C', 1, None), 16: ('C', 0, None), 17: ('C', 0, None),
    #         18: ('C', 0, None), 19: ('C', 0, None), 20: ('C', 1, None),
    #         21: ('C', 1, None), 22: ('C', 1, None), 23: ('C', 1, None)},
    #        {frozenset({17, 10}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({11, 12}): (1, None), frozenset({4, 5}): (1, None),
    #         frozenset({0, 5}): (1, None), frozenset({8, 7}): (1, None),
    #         frozenset({9, 4}): (1, None), frozenset({14, 15}): (1, None),
    #         frozenset({10, 11}): (1, None), frozenset({3, 4}): (1, None),
    #         frozenset({16, 17}): (1, None), frozenset({22, 23}): (1, None),
    #         frozenset({16, 15}): (1, None), frozenset({5, 6}): (1, None),
    #         frozenset({18, 19}): (1, None), frozenset({18, 3}): (1, None),
    #         frozenset({11, 14}): (1, None), frozenset({0, 1}): (1, None),
    #         frozenset({6, 7}): (1, None), frozenset({2, 21}): (1, None),
    #         frozenset({17, 18}): (1, None), frozenset({8, 13}): (1, None),
    #         frozenset({20, 21}): (1, None), frozenset({19, 20}): (1, None),
    #         frozenset({9, 10}): (1, None), frozenset({1, 2}): (1, None),
    #         frozenset({8, 9}): (1, None), frozenset({16, 23}): (1, None),
    #         frozenset({19, 22}): (1, None), frozenset({12, 13}): (1, None)})

    # # [C][C]
    # GRA = ({0: ('C', 0, None), 1: ('C', 0, None)},
    #        {frozenset({0, 1}): (1, None)})

    # [CH]=CCC#CC1C=CC=CC=1
    GRA = (
        {
            0: ("C", 1, None),
            1: ("C", 1, None),
            2: ("C", 2, None),
            3: ("C", 1, None),
            4: ("C", 0, None),
            5: ("C", 1, None),
            6: ("C", 1, None),
            7: ("C", 0, None),
            8: ("C", 1, None),
            9: ("C", 1, None),
            10: ("C", 0, None),
        },
        {
            frozenset({9, 6}): (1, None),
            frozenset({9, 10}): (1, None),
            frozenset({10, 7}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({0, 1}): (1, True),
            frozenset({3, 6}): (1, None),
            frozenset({8, 10}): (1, None),
            frozenset({2, 4}): (1, None),
            frozenset({3, 5}): (1, None),
            frozenset({8, 5}): (1, None),
            frozenset({4, 7}): (1, None),
        },
    )

    # # C#CC=C
    # GRA = ({0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 1, None),
    #         3: ('C', 2, None)},
    #        {frozenset({0, 1}): (1, None), frozenset({2, 3}): (1, None),
    #         frozenset({1, 2}): (1, None)})

    # # C#CC1C=CC=CC=1
    # GRA = ({3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
    #         6: ('C', 1, None), 7: ('C', 0, None), 8: ('C', 1, None),
    #         9: ('C', 1, None), 10: ('C', 0, None)},
    #        {frozenset({9, 6}): (1, None), frozenset({9, 10}): (1, None),
    #         frozenset({10, 7}): (1, None), frozenset({3, 6}): (1, None),
    #         frozenset({8, 10}): (1, None), frozenset({3, 5}): (1, None),
    #         frozenset({8, 5}): (1, None), frozenset({4, 7}): (1, None)})

    print(len(kekules(GRA)))
