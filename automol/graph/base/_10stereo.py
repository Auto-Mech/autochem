""" graph functions that depend on stereo assignments

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools

import numpy
from automol.graph.base._00core import (
    atom_keys,
    atom_stereo_keys,
    bond_stereo_keys,
    frozen,
    has_stereo,
    is_ts_graph,
    relabel,
    set_stereo_parities,
    stereo_parities,
    ts_reagents_graph_without_stereo,
    ts_reverse,
    without_stereo,
)
from automol.graph.base._03kekule import stereocenter_candidate_keys
from automol.graph.base._06geom import (
    geometry_atom_parity,
    geometry_bond_parity,
    geometry_correct_linear_vinyls,
    geometry_correct_nonplanar_pi_bonds,
    geometry_pseudorotate_atom,
    geometry_rotate_bond,
)
from automol.graph.base._07canon import (
    calculate_stereo,
    is_canonical_direction,
    is_canonical_enantiomer,
    parity_evaluator_reagents_from_ts_,
    refine_priorities,
    reflect_local_stereo,
    stereo_assignment_representation,
    stereocenter_keys_from_candidates,
    to_local_stereo,
)


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

    keys_pool = stereocenter_candidate_keys(gra)

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
            keys = stereocenter_keys_from_candidates(gra0, keys_pool, pri_dct, new=True)

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


def _remove_noncanonical_enantiomers_from_expansion(gps):
    """Remove non-canonical enantiomers from an expansion"""
    gps = list(gps)

    # a. Augment the list of graphs and priorities with local stereo graphs
    gpls = [(g, p, to_local_stereo(g, p)) for g, p in gps]

    # b. Find pairs of enantiomers and remove the non-canonical ones
    for ugpl, rgpl in itertools.combinations(gpls, r=2):
        ugra, upri_dct, uloc_gra = ugpl
        rgra, rpri_dct, rloc_gra = rgpl
        if rloc_gra == reflect_local_stereo(uloc_gra):
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
    if stereo and has_stereo(tsg):
        par_eval_ = parity_evaluator_reagents_from_ts_(tsg, prod=prod)
        gra, *_ = calculate_stereo(gra, par_eval_=par_eval_)
    return gra


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
