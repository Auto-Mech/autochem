""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numbers
from typing import Dict, Optional

import numpy
from automol.graph.base._00core import (
    atom_keys,
    atom_stereo_keys,
    backbone_keys,
    bond_stereo_keys,
    frozen,
    has_atom_stereo,
    has_stereo,
    invert_atom_stereo_parities,
    is_ts_graph,
    relabel,
    set_stereo_parities,
    stereo_keys,
    stereo_parities,
    ts_reagents_graph_without_stereo,
    ts_reverse,
    without_dummy_atoms,
    without_stereo,
)
from automol.graph.base._05stereo import (
    geometry_atom_parity,
    geometry_bond_parity,
    parity_evaluator_flip_from_graph,
    parity_evaluator_measure_from_geometry_,
    parity_evaluator_reactants_from_local_ts_graph_,
    parity_evaluator_read_from_graph,
    stereocenter_candidates,
    unassigned_stereocenter_keys_from_candidates,
)
from automol.graph.base._07geom import (
    geometry_correct_linear_vinyls,
    geometry_correct_nonplanar_pi_bonds,
    geometry_pseudorotate_atom,
    geometry_rotate_bond,
)
from automol.graph.base._08canon import (
    calculate_stereo,
    canonical_priorities,
    is_canonical_direction,
    is_canonical_enantiomer,
    refine_priorities,
    stereo_assignment_representation,
)
from automol.util import dict_


# # core functions
def expand_stereo(gra, symeq=False, enant=True):
    """Obtain all possible stereoisomers of a graph, ignoring its assignments

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param symeq: Include symmetrically equivalent stereoisomers?
    :type symeq: bool
    :param enant: Include all enantiomers, or only canonical ones?
    :type enant: bool
    :returns: a series of molecular graphs for the stereoisomers
    """
    # 1. Run the core stereo expansion algorithm
    gps = _expand_stereo_core(gra)

    # 3. If requested, filter out non-canonical enantiomers
    if not enant:
        gps = _remove_noncanonical_enantiomers_from_expansion(gps)

    # 4. If requested, filter out symmetry equivalents
    if not symeq:
        gps = _remove_symmetry_equivalents_from_expansion(gps)

    sgras = [sgra for sgra, _ in gps]
    sgras = tuple(sorted(sgras, key=frozen))
    return sgras


def _expand_stereo_core(gra):
    """Obtain all possible stereoisomers of a graph, ignoring its assignments

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param dir: If this is a TS graph, get priorities for the canonical direction?
    :type dir: bool, optional
    :returns: a series of molecular graphs for the stereoisomers
    """
    ts_ = is_ts_graph(gra)

    cand_dct = stereocenter_candidates(gra)

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
        gps = _select_ts_canonical_direction_priorities(gprs)
    else:
        gps = [(g, p) for g, p, _ in gprs]

    return gps


def _select_ts_canonical_direction_priorities(gprs):
    """Select select priorities for the canonical directions of each TS"""
    gps = []
    for ts_gra, pri_dct, rpri_dct in gprs:
        ts_rgra = ts_reverse(ts_gra)
        is_can_dir = is_canonical_direction(ts_gra, pri_dct, ts_rgra, rpri_dct)
        gps.append((ts_gra, pri_dct if is_can_dir else rpri_dct))
    return gps


def to_local_stereo(gra, pri_dct=None):
    """Convert canonical stereo parities to local ones

    Local parities are based directly on the key values of neighboring
    atoms, whereas canonical parities are based on their canonical
    priorities.  Consequently, local parities are specific to the
    particular way the graph is labeled, so the graph cannot be relabeled
    without corrupting stereo information, but they are useful for
    temporarily decoupling stereo parities from each other as the graph is
    manipulated in other ways.

    Note that, for consistency with InChI and other systems, hydrogen keys
    are treated as having lowest priority. This is done by setting their
    sort value to negative infinity.

    For TS graphs, canonical priorities are given with respect to the canonical
    TS direction. To avoid dependence of local parities on the canonical
    direction, we must *undo* the direction reversal that occurs during the
    canonical priority calculation when generating local parities. As far as I
    am aware, the *only* case affected by this is Sn2 reactions, so that is all
    that I have implemented here.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :returns: molecular graph with local stereo parities
    :rtype: automol graph data structure
    """
    can_gra = gra
    if has_stereo(can_gra):
        pri_dct_ = (
            None if pri_dct is None else dict_.by_key(pri_dct, backbone_keys(can_gra))
        )
        loc_gra, *_ = calculate_stereo(
            can_gra,
            par_eval_=parity_evaluator_flip_from_graph,
            can_par_eval_=parity_evaluator_read_from_graph,
            pri_dct=pri_dct_,
        )
    else:
        loc_gra = can_gra

    return loc_gra


def _remove_noncanonical_enantiomers_from_expansion(gps):
    """Remove non-canonical enantiomers from an expansion"""
    gps = list(gps)

    # a. Augment the list of graphs and priorities with local stereo graphs
    gpls = [(g, p, to_local_stereo(g, p)) for g, p in gps]

    # b. Find pairs of enantiomers and remove the non-canonical ones
    for ugpl, rgpl in itertools.combinations(gpls, r=2):
        ugra, upri_dct, uloc_gra = ugpl
        rgra, rpri_dct, rloc_gra = rgpl
        if rloc_gra == invert_atom_stereo_parities(uloc_gra):
            is_can = is_canonical_enantiomer(ugra, upri_dct, rgra, rpri_dct)

            if is_can is True:
                gps.remove((rgra, rpri_dct))
            elif is_can is False:
                gps.remove((ugra, upri_dct))

    return gps


def _remove_symmetry_equivalents_from_expansion(gps):
    """Remove symmetry-equivalent stereoisomers from an expansion"""
    gps0 = gps
    gps = []
    seen_reps = []
    for gra, pri_dct in gps0:
        rep = stereo_assignment_representation(gra, pri_dct)
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


def has_fleeting_atom_or_bond_stereo(tsg, strict: bool = True) -> (bool, bool):
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
def stereo_corrected_geometry(gra, geo, geo_idx_dct=None, local_stereo=False):
    """Obtain a geometry corrected for stereo parities based on a graph

    :param gra: molecular graph with stereo parities
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param local_stereo: is this graph using local instead of canonical
        stereo?
    :type local_stereo: bool
    :returns: a molecular geometry with corrected stereo
    """
    sgr = gra if local_stereo else to_local_stereo(gra)
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(atm_keys)} if geo_idx_dct is None else geo_idx_dct
    )
    gra = relabel(gra, geo_idx_dct)

    par_dct = stereo_parities(sgr)
    bnd_keys = bond_stereo_keys(sgr)
    atm_keys = atom_stereo_keys(sgr)

    # 1. Correct linear vinyl groups
    geo = geometry_correct_linear_vinyls(gra, geo)
    geo = geometry_correct_nonplanar_pi_bonds(gra, geo)

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

    return geo


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
    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    par_eval_ = parity_evaluator_measure_from_geometry_(
        geo, local_stereo=local_stereo, geo_idx_dct=geo_idx_dct
    )
    can_par_eval_ = None

    if local_stereo:
        # If requesting local stereo, we need an auxiliary canonical parity evaluator
        can_par_eval_ = parity_evaluator_measure_from_geometry_(
            geo, geo_idx_dct=geo_idx_dct
        )

    gra, *_ = calculate_stereo(gra, par_eval_=par_eval_, can_par_eval_=can_par_eval_)

    return gra


def from_local_stereo(gra, pri_dct=None):
    """Convert local stereo parities to canonical ones

    :param gra: molecular graph with local stereo parities
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :returns: molecular graph with canonical stereo parities
    :rtype: automol graph data structure
    """
    loc_gra = gra

    if has_stereo(loc_gra):
        pri_dct_ = (
            None if pri_dct is None else dict_.by_key(pri_dct, backbone_keys(loc_gra))
        )
        can_gra, *_ = calculate_stereo(
            loc_gra,
            par_eval_=parity_evaluator_flip_from_graph,
            pri_dct=pri_dct_,
        )
    else:
        can_gra = loc_gra

    return can_gra


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
