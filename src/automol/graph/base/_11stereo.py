""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import numbers
from typing import Dict, Optional, Tuple

import more_itertools as mit
import numpy
from phydat import phycon

from ... import util
from ...geom import base as geom_base
from ...util import dict_
from ._00core import (
    align_with_geometry,
    atom_neighbor_atom_keys,
    atom_stereo_keys,
    atomic_numbers,
    bond_stereo_keys,
    equivalent_without_dummy_atoms,
    frozen,
    has_atom_stereo,
    has_stereo,
    invert_atom_stereo_parities,
    is_ts_graph,
    relabel,
    set_stereo_parities,
    stereo_keys,
    stereo_parities,
    ts_reacting_atom_keys,
    ts_reagents_graph_without_stereo,
    ts_reverse,
    without_dummy_atoms,
    without_stereo,
)
from ._02algo import (
    branch_atom_keys,
    branch_dict,
    ring_systems_atom_keys,
)
from ._03kekule import linear_atom_keys
from ._05stereo import (
    geometry_atom_parity,
    geometry_bond_parity,
    parity_evaluator_measure_from_geometry_,
    parity_evaluator_reactants_from_local_ts_graph_,
    stereoatom_bridgehead_pairs,
    stereocenter_candidates,
    stereocenter_candidates_grouped,
    unassigned_stereocenter_keys_from_candidates,
)
from ._07geom import (
    geometry_correct_linear_vinyls,
    geometry_correct_nonplanar_pi_bonds,
    geometry_rotate_bond,
)
from ._08canon import (
    calculate_stereo,
    canonical_priorities,
    from_local_stereo,
    is_canonical_direction,
    is_canonical_enantiomer,
    refine_priorities,
    stereo_assignment_representation,
    to_local_stereo,
)


# # core functions
def expand_stereo(
    gra,
    symeq: bool = False,
    enant: bool = True,
    strained: bool = False,
    rcts_gra: object | None = None,
    prds_gra: object | None = None,
):
    """Obtain all possible stereoisomers of a graph, ignoring its assignments.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symeq: Include symmetrically equivalent stereoisomers?
    :param enant: Include all enantiomers, or only canonical ones?
    :param strained: Include stereoisomers which are too strained to exist?
    :param rcts_gra: For TS graphs, optionally specify a reactants graph to match
    :param prds_gra: For TS graphs, optionally specify a products graph to match
    :returns: a series of molecular graphs for the stereoisomers
    """
    cand_dct = stereocenter_candidates(gra)
    rcand_dct = stereocenter_candidates(ts_reverse(gra))

    # 1. Run the core stereo expansion algorithm
    gps = _expand_stereo_core(
        gra, cand_dct, rcand_dct=rcand_dct, rcts_gra=rcts_gra, prds_gra=prds_gra
    )

    # 2. If requested, filter out strained stereoisomers
    # (If not, filter out only strained ring H-migrations, which are guaranteed to cause
    # problems)
    gps = _remove_strained_stereoisomers_from_expansion(
        gps, cand_dct=cand_dct, migration_only=strained
    )

    # 3. If requested, filter out non-canonical enantiomers
    if not enant:
        gps = _remove_noncanonical_enantiomers_from_expansion(gps, cand_dct=cand_dct)

    # 4. If requested, filter out symmetry equivalents
    if not symeq:
        gps = _remove_symmetry_equivalents_from_expansion(gps, cand_dct=cand_dct)

    sgras = [sgra for sgra, _ in gps]
    sgras = tuple(sorted(sgras, key=frozen))
    return sgras


def _expand_stereo_core(
    gra,
    cand_dct,
    rcand_dct: dict | None = None,
    rcts_gra: object | None = None,
    prds_gra: object | None = None,
):
    """Obtain all possible stereoisomers of a graph, ignoring its assignments.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param rcts_gra: For TS graphs, optionally specify a reactants graph to match
    :param prds_gra: For TS graphs, optionally specify a products graph to match
    :returns: a series of molecular graphs for the stereoisomers
    """
    rcts_gra = None if rcts_gra is None or not has_stereo(rcts_gra) else rcts_gra
    prds_gra = None if prds_gra is None or not has_stereo(prds_gra) else prds_gra

    ts_ = is_ts_graph(gra)

    bools = (False, True)

    gra0 = without_stereo(gra)
    gprs0 = None
    gprs = [(gra0, None, None)]

    # 1. Expand all possible stereoisomers, along with their priority mappings
    while gprs0 != gprs:
        gprs0 = gprs
        gprs = []

        for gra0, pri_dct0, rpri_dct0 in gprs0:
            # a. Refine priorities based on current assignments
            pri_dct = refine_priorities(gra0, pri_dct=pri_dct0)
            rpri_dct = (
                refine_priorities(ts_reverse(gra0), pri_dct=rpri_dct0) if ts_ else None
            )

            # c. Find stereogenic atoms and bonds based on current priorities
            keys = unassigned_stereocenter_keys_from_candidates(gra0, cand_dct, pri_dct)

            # d. Assign True/False parities in all possible ways
            for pars in itertools.product(bools, repeat=len(keys)):
                gra = set_stereo_parities(gra0, dict(zip(keys, pars)))
                gprs.append((gra, pri_dct, rpri_dct))

    if ts_:
        gps = _select_ts_canonical_direction_priorities(
            gprs, fcand_dct=cand_dct, rcand_dct=rcand_dct
        )
    else:
        gps = [(g, p) for g, p, _ in gprs]

    if rcts_gra is not None:
        gps = [
            (g, p)
            for g, p in gps
            if equivalent_without_dummy_atoms(ts_reactants_graph(g), rcts_gra)
        ]

    if prds_gra is not None:
        gps = [
            (g, p)
            for g, p in gps
            if equivalent_without_dummy_atoms(ts_products_graph(g), prds_gra)
        ]

    return gps


def _select_ts_canonical_direction_priorities(gprs, fcand_dct: dict, rcand_dct: dict):
    """Select select priorities for the canonical directions of each TS"""
    gps = []
    for ts_gra, pri_dct, rpri_dct in gprs:
        ts_rgra = ts_reverse(ts_gra)
        is_can_dir = is_canonical_direction(
            ts_gra, pri_dct, ts_rgra, rpri_dct, fcand_dct=fcand_dct, rcand_dct=rcand_dct
        )
        gps.append((ts_gra, pri_dct if is_can_dir else rpri_dct))
    return gps


def _remove_strained_stereoisomers_from_expansion(
    gps, cand_dct, migration_only: bool = False
):
    """Remove strained stereoisomers from an expansion.

    :param migration_only: Only remove strained H-migration bridgeheads?
    """
    if not gps:
        return gps

    gps0 = list(gps)
    gra = without_stereo(gps0[0][0])
    bhp_dct = stereoatom_bridgehead_pairs(gra, cand_dct, migration_only=migration_only)

    gps = gps0.copy()
    for gra, pri_dct in gps0:
        par_dct = stereo_parities(gra)
        can_nkeys_dct, _ = stereocenter_candidates_grouped(cand_dct, pri_dct=pri_dct)
        for (key1, key2), (conn_nkeys1, conn_nkeys2) in bhp_dct.items():
            free_nkeys1 = tuple(k for k in cand_dct[key1] if k not in conn_nkeys1)
            free_nkeys2 = tuple(k for k in cand_dct[key2] if k not in conn_nkeys2)

            srt_nkeys1 = free_nkeys1 + conn_nkeys1
            srt_nkeys2 = free_nkeys2 + conn_nkeys2
            can_nkeys1 = can_nkeys_dct[key1]
            can_nkeys2 = can_nkeys_dct[key2]

            par1, par2 = map(par_dct.get, (key1, key2))

            # If the parities relative to the above ordering are not opposite, then the
            # configuration of the bridgehead pair is strained
            # Exception: Adjacent bridgehead pairs will have groups on the same side
            if par1 is not None and par2 is not None:
                sgn1 = util.is_odd_permutation(srt_nkeys1, can_nkeys1)
                sgn2 = util.is_odd_permutation(srt_nkeys2, can_nkeys2)
                is_strained = not par1 ^ par2 ^ sgn1 ^ sgn2
                if is_strained:
                    gps.remove((gra, pri_dct))
                    break

    return gps


def _remove_noncanonical_enantiomers_from_expansion(gps, cand_dct: dict):
    """Remove non-canonical enantiomers from an expansion."""
    gps = list(gps)

    # a. Augment the list of graphs and priorities with local stereo graphs
    gpls = [(g, p, to_local_stereo(g, p)) for g, p in gps]

    # b. Find pairs of enantiomers and remove the non-canonical ones
    for ugpl, rgpl in itertools.combinations(gpls, r=2):
        ugra, upri_dct, uloc_gra = ugpl
        rgra, rpri_dct, rloc_gra = rgpl
        if rloc_gra == invert_atom_stereo_parities(uloc_gra):
            is_can = is_canonical_enantiomer(
                ugra, upri_dct, rgra, rpri_dct, cand_dct=cand_dct
            )

            if is_can is True:
                gps.remove((rgra, rpri_dct))
            elif is_can is False:
                gps.remove((ugra, upri_dct))

    return gps


def _remove_symmetry_equivalents_from_expansion(gps, cand_dct: dict):
    """Remove symmetry-equivalent stereoisomers from an expansion."""
    gps0 = gps
    gps = []
    seen_reps = []
    for gra, pri_dct in gps0:
        rep = stereo_assignment_representation(gra, pri_dct, cand_dct=cand_dct)
        if rep not in seen_reps:
            gps.append((gra, pri_dct))
            seen_reps.append(rep)
    return gps


# # TS functions
def expand_reaction_stereo(ts_gra, flat: bool = False):
    """Obtain all possible stereoisomeric pathways of a reaction, grouped by reactants
    and products (unless requesting `flat` expansion)

    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
    :param flat: Return a flat list, instead of grouping by reagents?, default False
    :type flat: bool, optional
    :returns: A list of triples: a reactants graph, a products graph, and a list of TS
        graphs associated with them
    :rtype: List[Tuple[graph, graph, List[graph]]]
    """

    def _from_groupby_result(groupby_result):
        """For formatting groupby results"""
        _, group = groupby_result
        rcts_gras, prds_gras, ts_gras = zip(*group)
        return (rcts_gras[0], prds_gras[0], ts_gras)

    # Allow *all* possibilities for the TS, including symmetry equivalent ones
    ts_sgras = expand_stereo(ts_gra, symeq=True, enant=True)

    rxn_sgras_lst = []
    for ts_sgra in ts_sgras:
        rcts_sgra = ts_reactants_graph(ts_sgra)
        prds_sgra = ts_products_graph(ts_sgra)
        rxn_sgras_lst.append((rcts_sgra, prds_sgra, ts_sgra))

    rxn_sgras_lst = sorted(rxn_sgras_lst, key=lambda gs: tuple(map(frozen, gs)))

    if not flat:
        rxn_sgras_groupby_iter = itertools.groupby(
            rxn_sgras_lst, key=lambda gs: tuple(map(frozen, gs[:-1]))
        )
        rxn_sgras_lst = list(map(_from_groupby_result, rxn_sgras_groupby_iter))
    return rxn_sgras_lst


def ts_reactants_graph(tsg, stereo=True, dummy=True):
    """Get the reactants from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default True
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return ts_reagents_graph(tsg, prod=False, stereo=stereo, dummy=dummy)


def ts_products_graph(tsg, stereo=True, dummy=True):
    """Get the products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default True
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    return ts_reagents_graph(tsg, prod=True, stereo=stereo, dummy=dummy)


def ts_reagents_graph(tsg, prod=False, stereo=True, dummy=True):
    """Get the reactants or products from a TS graph

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Do this for the products, instead of the reactants?
    :type prod: bool
    :param stereo: Keep stereo, even though it is invalid?
    :type stereo: bool, optional
    :param dummy: Keep dummy atoms? default True
    :type dummy: bool, optional
    :returns: The TS graph, without reacting bond orders
    :rtype: automol graph data structure
    """
    gra = ts_reagents_graph_without_stereo(tsg, prod=prod, dummy=dummy)
    if stereo:
        tsg = ts_reverse(tsg) if prod else tsg
        loc_tsg = to_local_stereo(tsg)
        par_eval_ = parity_evaluator_reactants_from_local_ts_graph_(loc_tsg)
        gra, *_ = calculate_stereo(gra, par_eval_=par_eval_)
    return gra


def ts_fleeting_stereocenter_keys(tsg, strict: bool = True):
    """Find keys to fleeting atom stereocenters in a TS graph

    A fleeting stereocenter is a stereocenter which is not present in the reactants or
    products

    A non-strict fleeting stereocenter is implicitly present in the reactants or
    products, whose configuration determines one of the reacant or product stereocenters
    (Generally, a TS atom that determines a reactant or product bond)

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param strict: Only include strict fleeting stereocenters?
    :type strict: bool, optional
    :returns: The keys of fleeting atom stereocenters
    :rtype: List[int]
    """
    rcts_gra = ts_reactants_graph(tsg)
    prds_gra = ts_products_graph(tsg)
    r_akeys = atom_stereo_keys(rcts_gra)
    p_akeys = atom_stereo_keys(prds_gra)
    r_bkeys = bond_stereo_keys(rcts_gra)
    p_bkeys = bond_stereo_keys(prds_gra)
    rp_keys = r_akeys | r_bkeys | p_akeys | p_bkeys

    keys = {k for k in stereo_keys(tsg) if k not in rp_keys}

    if strict:
        # Missing case captured here: In some cases, the configuration of an adjacent
        # pair of atom stereocenters is implied by the combination of one of the atoms
        # being a stereocenter for the reactant or product and the bond being a
        # stereocenter for the reactant or product (both must be present)
        # In this case, the stereochemsitry of the reaction site is fully specified by
        # the reactants and products and neither atom is truly a fleeting stereocenter
        for bkey in r_bkeys | p_bkeys:
            if any(akey in bkey for akey in r_akeys | p_akeys):
                keys -= bkey

    return frozenset(keys)


def has_fleeting_atom_or_bond_stereo(tsg, strict: bool = True) -> Tuple[bool, bool]:
    """Does this graph have fleeting atom or bond stereo?

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param strict: Only include strict fleeting stereocenters?
    :type strict: bool, optional
    :return: A boolean flag for each type (atom and bond)
    :rtype: (bool, bool)
    """
    keys = ts_fleeting_stereocenter_keys(tsg, strict=strict)
    akeys = {k for k in keys if isinstance(k, numbers.Number)}
    bkeys = keys - akeys
    return bool(akeys), bool(bkeys)


# # stereo correction
def stereo_corrected_geometry(
    gra, geo, geo_idx_dct=None, local_stereo: bool = False, lin_ts_bonds: bool = False
):
    """Obtain a geometry corrected for stereo parities based on a graph.

    :param gra: molecular graph with stereo parities
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param local_stereo: is this graph using local instead of canonical stereo?
    :param lin_ts_bonds: For TS graphs, apply stereo correction bonds that are linear
        for the reactants or products?
    :returns: a molecular geometry with corrected stereo
    """
    if geo is None:
        return None

    # Align the graph and the geometry keys/indices
    gra, geo, *_, idx_dct = align_with_geometry(gra, geo, (), geo_idx_dct)

    # Note: The alignment preserves graph atom ordering, so the local stereo does not
    # change, i.e. this could be done before or after
    gra = gra if local_stereo else to_local_stereo(gra)

    # Determine atoms expected to be linear for the TS, if requested
    excl_keys = set() if lin_ts_bonds else linear_atom_keys(gra)

    par_dct = stereo_parities(gra)
    atm_keys = atom_stereo_keys(gra)
    bnd_keys = bond_stereo_keys(gra)
    bnd_keys = {bk for bk in bnd_keys if not bk & excl_keys}

    # 1. Correct linear vinyl groups
    geo = geometry_correct_linear_vinyls(gra, geo, excl_keys=excl_keys)
    geo = geometry_correct_nonplanar_pi_bonds(gra, geo, excl_keys=excl_keys)

    # 2. If there is a single, wrong atom stereocenter, simply reflect the geometry
    if len(atm_keys) == 1:
        (atm_key,) = atm_keys
        curr_par = geometry_atom_parity(gra, geo, atm_key)
        if curr_par != par_dct[atm_key]:
            geo = geom_base.reflect_coordinates(geo)

    # 3. Loop over stereo-sites making corrections where needed
    for bnd_key in bnd_keys:
        curr_par = geometry_bond_parity(gra, geo, bnd_key)
        if curr_par != par_dct[bnd_key]:
            geo = geometry_rotate_bond(gra, geo, bnd_key, numpy.pi)

    for atm_key in atm_keys:
        curr_par = geometry_atom_parity(gra, geo, atm_key)
        if curr_par != par_dct[atm_key]:
            geo = geometry_pseudorotate_atom(gra, geo, atm_key)

            assert geo is not None, (
                f"Failed to correct the following geometry:"
                f"\ngeo:\n{geo}\ngra:\n{gra}"
            )

    # Restore the original atom ordering of the geometry
    return geom_base.reorder(geo, idx_dct)


def set_stereo_from_geometry(gra, geo, local_stereo=False, geo_idx_dct=None):
    """Determine stereo parities from a geometry

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :returns: molecular graph with stereo parities set from geometry;
        parities already present will be wiped out
    :rtype: automol graph data structure
    """
    # Align the graph and the geometry keys/indices
    gra, geo, _, key_dct, _ = align_with_geometry(gra, geo, (), geo_idx_dct)

    par_eval_ = parity_evaluator_measure_from_geometry_(geo, local_stereo=local_stereo)
    can_par_eval_ = None

    if local_stereo:
        # If requesting local stereo, we need an auxiliary canonical parity evaluator
        can_par_eval_ = parity_evaluator_measure_from_geometry_(geo)

    gra, *_ = calculate_stereo(gra, par_eval_=par_eval_, can_par_eval_=can_par_eval_)

    # Restore the original atom keys of the graph
    return relabel(gra, key_dct)


def reflect(gra):
    """Calculate new parities that would result from geometric reflection

    To replicate the effect of reflecting the geometry, we convert to local
    stereo, invert parities, and then convert back.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    """
    if has_atom_stereo(gra):
        loc_gra = to_local_stereo(gra)
        loc_gra = invert_atom_stereo_parities(loc_gra)
        gra = from_local_stereo(loc_gra)
    return gra


def unassigned_stereocenter_keys(
    gra,
    atom: bool = True,
    bond: bool = True,
    pri_dct: Optional[Dict[int, int]] = None,
):
    """Find keys to unassigned stereocenters in this graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param atom: Include atom stereocenters? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenters? defaults to True
    :type bond: bool, optional
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: Optional[Dict[int, int]]
    :returns: Keys to stereogenic atoms and bonds which are unassigned
    :rtype: frozenset
    """
    gra = without_dummy_atoms(gra)
    pri_dct = canonical_priorities(gra) if pri_dct is None else pri_dct
    cand_dct = stereocenter_candidates(gra, atom=atom, bond=bond)
    ste_keys = unassigned_stereocenter_keys_from_candidates(gra, cand_dct, pri_dct)
    return ste_keys


def geometry_pseudorotate_atom(
    gra, geo, key, ang=numpy.pi, degree=False, geo_idx_dct=None
):
    r"""Pseudorotate an atom in a molecular geometry by a certain amount.

    'Pseudorotate' here means to rotate all but two of the atom's neighbors, which can
    be used to invert/correct stereochemistry at an atom:

        1   2                                     1   2
         \ /                                       \ /
          C--3   = 1,4 pseudorotation by pi =>   3--C
          |                                         |
          4                                         4

    The two fixed atoms will be chosen to prevent the structural 'damage' from the
    rotation as much as possible. For example, atoms in rings will be favored to be
    fixed.

    If such a choice is not possible -- for example, if three or more neighbors are
    locked into connected rings -- then no geometry will be returned.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param key: The graph key of the atom to be rotated
    :type key: frozenset[int]
    :param ang: The angle of rotation (in radians, unless `degree = True`)
    :type ang: float
    :param degree: Is the angle of rotation in degrees?, default False
    :type degree: bool
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    ang = ang * phycon.DEG2RAD if degree else ang
    # Align the graph and the geometry keys/indices
    gra, geo, (key,), _, idx_dct = align_with_geometry(gra, geo, (key,), geo_idx_dct)

    gra_reac = ts_reactants_graph(gra)
    gra_prod = ts_products_graph(gra)

    rxn_keys = ts_reacting_atom_keys(gra)
    rsy_keys_lst = []
    for rgra in (gra_reac, gra_prod):
        rsy_keys_lst.extend(ring_systems_atom_keys(rgra, lump_spiro=False))
    rsy_keys_lst = list(set(rsy_keys_lst))
    nkeys = atom_neighbor_atom_keys(gra, key)
    # Gather neighbors connected in a ring system
    ring_nkey_sets = [nkeys & ks for ks in rsy_keys_lst if nkeys & ks]
    ring_nkey_sets = sorted(ring_nkey_sets, key=len, reverse=True)
    # Gather the remaining neighbors
    rem_nkeys = [k for k in nkeys if not any(k in ks for ks in ring_nkey_sets)]
    # Sort the remaining neighbors by branch size and atomic number
    anum_dct = atomic_numbers(gra)
    size_dct = dict_.transform_values(branch_dict(gra, key), len)
    sort_dct = {k: (k in rxn_keys, -size_dct[k], -anum_dct[k]) for k in rem_nkeys}
    rem_nkeys = sorted(rem_nkeys, key=sort_dct.get)

    # Now, put the two lists together
    nkey_sets = ring_nkey_sets + [{k} for k in rem_nkeys]

    # Now, find a pair of atoms to keep fixed
    found_pair = False
    for nkeys1, nkeys2 in mit.pairwise([*nkey_sets, set()]):
        if len(nkeys1) == 2 or len(nkeys1 | nkeys2) == 2:
            found_pair = True
            nkey1, nkey2, *_ = list(nkeys1) + list(nkeys2)
            break

    if not found_pair:
        return None

    # Determine the rotational axis as the unit bisector between the fixed pair
    xyz, nxyz1, nxyz2 = geom_base.coordinates(geo, idxs=(key, nkey1, nkey2))
    rot_axis = util.vector.unit_bisector(nxyz1, nxyz2, orig_xyz=xyz)

    # Identify the remaining keys to be rotated
    rot_nkeys = nkeys - {nkey1, nkey2}
    rot_keys = set(itertools.chain(*(branch_atom_keys(gra, key, k) for k in rot_nkeys)))

    # For bicyclic TSs, we need to avoid rotating the whole molecule when there is a
    # forming ring
    if nkey1 in rot_keys or nkey2 in rot_keys:
        rot_keys = set(
            itertools.chain(
                *(
                    branch_atom_keys(gra_reac, key, k, return_missing_neighbor=True)
                    for k in rot_nkeys
                )
            )
        )

    geo = geom_base.rotate(geo, rot_axis, ang, orig_xyz=xyz, idxs=rot_keys)

    # Restore the original atom ordering of the geometry
    return geom_base.reorder(geo, idx_dct)
