""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import itertools
import numbers
from collections import abc
from typing import Any, Callable, Dict, List, Optional

import numpy
from automol import util
from automol.graph.base._00core import (
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    atom_stereo_sorted_neighbor_keys,
    atom_symbols,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    backbone_keys,
    bond_orders,
    bond_stereo_keys,
    bond_stereo_parities,
    bond_stereo_sorted_neighbor_keys,
    explicit,
    has_atom_stereo,
    has_stereo,
    implicit,
    is_ts_graph,
    mass_numbers,
    nonbackbone_hydrogen_keys,
    relabel as relabel_,
    set_atom_stereo_parities,
    set_stereo_parities,
    stereo_keys,
    stereo_parities,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reactants_graph_without_stereo,
    ts_reverse,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_pi_bonds,
    without_stereo,
)
from automol.graph.base._02algo import connected_components, is_connected
from automol.graph.base._03kekule import stereocenter_candidate_keys
from automol.graph.base._04ts import (
    constrained_1_2_insertion_local_parities,
    sn2_local_stereo_reversal_flips,
    vinyl_addition_local_parities,
)
from automol.graph.base._06geom import geometry_atom_parity, geometry_bond_parity
from automol.util import dict_
from phydat import ptab

# Special types
ParityEvaluator = Callable[[Any, Dict[int, int], List[int], bool], Dict[int, int]]


# # canonical key functions
def canonical_enantiomer(gra, relabel: bool = True):
    """Determine the canonical graph of the canonical enantiomer.

    Graphs with stereo will be reflected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param relabel: Relabel the graph with canonical keys? defaults to True
    :type relabel: bool, optional
    :returns: a canonicalized graph of the canonical enantiomer, along with
        a boolean flag indicating whether or not the graph has been
        reflected; `True` indicates it has been, `False` indicates it
        hasn't, and `None` indicates that it isn't an enantiomer
    :rtype: (automol graph data structure, bool)
    """
    gra, _, is_can_enant = canonical_enantiomer_with_keys(gra, relabel=relabel)
    return gra, is_can_enant


def canonical_enantiomer_with_keys(gra, relabel: bool = False):
    """Determine the canonical graph of the canonical enantiomer, along with
    the canonical key mapping.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param relabel: Relabel the graph with canonical keys? defaults to True
    :type relabel: bool, optional
    :returns: a canonicalized graph of the canonical enantiomer, along with a boolean
        flag indicating whether or not the input graph is this enantiomer;
        `True` indicates it is, `False` indicates it isn't, and `None` indicates that it
        is achiral
    :rtype: (automol graph data structure, bool)
    """
    if not has_atom_stereo(gra):
        can_key_dct = canonical_keys(gra, backbone_only=False)
        is_can_enant = None
    else:
        assert is_connected(gra), f"Requires a connected graph.\n{gra}"
        # Calculate canonical keys for the unreflected graph while converting
        # to the local stereo representation
        loc_gra, _, pri_dct, *_ = calculate_stereo(
            gra,
            par_eval_=parity_evaluator_flip_local_(),
            can_par_eval_=parity_evaluator_read_canonical_(),
            backbone_only=False,
        )

        # Reflect the graph in the local stereo representation
        rloc_gra = reflect_local_stereo(loc_gra)

        # Determine canonical keys for the reflected graph while converting
        # back to the canonical stereo representation
        rgra, _, rpri_dct, *_ = calculate_stereo(
            rloc_gra, backbone_only=False, par_eval_=parity_evaluator_flip_local_()
        )

        is_can_enant = is_canonical_enantiomer(gra, pri_dct, rgra, rpri_dct)

        if is_can_enant is False:
            gra = rgra
            pri_dct = rpri_dct

        can_key_dct = break_priority_ties(gra, pri_dct)

    if relabel:
        gra = relabel_(gra, can_key_dct)

    return gra, can_key_dct, is_can_enant


def canonical_ts_direction(gra):
    """If this is a TS graph, determine the canonical direction for it

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: a graph of the canonical TS direction, along with a boolean flag
        indicating whether or not the input graph has this direction;
        `True` indicates it is, `False` indicates it isn't, and `None` indicates that it
        isn't a TS graph
    :rtype: (automol graph data structure, bool)
    """
    if is_ts_graph(gra):
        *_, is_can_dir = calculate_stereo(gra)
        if is_can_dir is False:
            gra = ts_reverse(gra)
    else:
        is_can_dir = None

    return gra, is_can_dir


def canonical(gra):
    """A graph relabeled with canonical keys

    Stereo parities in the graph are assumed to be canonical.

    Requires a connected graph

    :param gra: a connected molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :returns: a new molecular graph with canonical keys; if explicit
        hydrogens are included, they will be relabeled as well
    """
    can_key_dct = canonical_keys(gra, backbone_only=False)
    return relabel_(gra, can_key_dct)


def canonical_keys(gra, backbone_only=True):
    """Determine canonical keys for this graph.

    Stereo parities in the graph are assumed to be canonical.

    Requires a connected graph

    :param gra: a connected molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :returns: a dictionary of canonical keys by atom key
    :rtype: dict[int: int]
    """
    assert is_connected(gra), f"Requires a connected graph.\n{gra}"
    atm_par_dct0 = atom_stereo_parities(gra)
    bnd_par_dct0 = bond_stereo_parities(gra)

    gra, _, pri_dct, *_ = calculate_stereo(gra, backbone_only=backbone_only)
    can_key_dct = break_priority_ties(gra, pri_dct)

    atm_par_dct = atom_stereo_parities(gra)
    bnd_par_dct = bond_stereo_parities(gra)

    assert atm_par_dct == atm_par_dct0, (
        f"Atom stereo parities don't match input. Something is wrong:\n"
        f"input: {atm_par_dct0}\n"
        f"return: {atm_par_dct}\n"
    )

    assert bnd_par_dct == bnd_par_dct0, (
        f"Bond stereo parities don't match input. Something is wrong:\n"
        f"input: {bnd_par_dct0}\n"
        f"return: {bnd_par_dct}\n"
    )

    return can_key_dct


def stereo_assignment_representation(gra, pri_dct: dict):
    """Generate a representation of a stereo assignment, for checking for
    symmetric equivalence or for determining a canonical enantiomer

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :returns: A canonical representation of the assignment
    """
    atm_keys = sorted(atom_stereo_keys(gra), key=pri_dct.get)
    bnd_keys = sorted(bond_stereo_keys(gra), key=lambda x: sorted(map(pri_dct.get, x)))

    rep = tuple(
        dict_.values_by_key(atom_stereo_parities(gra), atm_keys)
        + dict_.values_by_key(bond_stereo_parities(gra), bnd_keys)
    )
    return rep


def ts_direction_representation(tsg, pri_dct: dict):
    """Generate a representation of the reaction direction to determine the
    canonical direction of a TS graph

    The reaction representation consists of two pieces:

    1. Canonical representations of the forming and breaking bonds,
    respectively.
    2. Canonical representations of the stereo assignments.

    :param tsg: A TS graph
    :type tsg: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :returns: A canonical representation of the TS reaction
    """
    frm_keys = ts_forming_bond_keys(tsg)
    brk_keys = ts_breaking_bond_keys(tsg)
    # Rep value 1: Number of bonds broken and formed
    rxn_rep1 = (len(frm_keys), len(brk_keys))
    # Rep value 2: Canonical keys of bonds broken and formed
    rxn_rep2 = (
        sorted(sorted(map(pri_dct.get, k)) for k in frm_keys),
        sorted(sorted(map(pri_dct.get, k)) for k in brk_keys),
    )
    rep = (rxn_rep1, rxn_rep2)
    return rep


def is_canonical_enantiomer(ugra, upri_dct, rgra, rpri_dct):
    """Is this enantiomer the canonical one?

    :param ugra: An unreflected molecular graph
    :type ugra: automol graph data structure
    :param upri_dct: A dictionary mapping atom keys to priorities for `ugra`
    :type upri_dct: dict
    :param rgra: A reflected molecular graph
    :type rgra: automol graph data structure
    :param rpri_dct: A dictionary mapping atom keys to priorities for `rgra`
    :type rpri_dct: dict
    :returns: `True` if it is, `False` if it isn't, and `None` if it isn't an
        enantiomer
    """
    urep = stereo_assignment_representation(ugra, upri_dct)
    rrep = stereo_assignment_representation(rgra, rpri_dct)
    return True if (urep < rrep) else False if (urep > rrep) else None


def is_canonical_direction(ftsg, fpri_dct, rtsg, rpri_dct):
    """Is this TS direction the canonical one?

    :param ftsg: A TS graph in the forward direction
    :type ftsg: automol graph data structure
    :param fpri_dct: A dictionary mapping atom keys to priorities for the
        forward direction
    :type fpri_dct: dict
    :param rtsg: A TS graph in the reverse direction
    :type rtsg: automol graph data structure
    :param rpri_dct: A dictionary mapping atom keys to priorities for the
        reverse direction
    :type rpri_dct: dict
    :returns: A canonical representation of the TS reaction
    """
    frep = ts_direction_representation(ftsg, fpri_dct)
    rrep = ts_direction_representation(rtsg, rpri_dct)
    return True if (frep < rrep) else False if (frep > rrep) else None


# # canonical stereo functions
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
    pri_dct = (
        canonical_priorities(gra, backbone_only=False) if pri_dct is None else pri_dct
    )
    ste_keys = stereocenter_keys_from_priorities(
        gra, pri_dct, atom=atom, bond=bond, new=True
    )
    return ste_keys


def reflect(gra):
    """Calculate new parities that would result from geometric reflection

    To replicate the effect of reflecting the geometry, we convert to local
    stereo, invert parities, and then convert back.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    """
    if has_atom_stereo(gra):
        loc_gra = to_local_stereo(gra)
        loc_gra = reflect_local_stereo(loc_gra)
        gra = from_local_stereo(loc_gra)
    return gra


def reflect_local_stereo(gra):
    """Reflect a graph with local stereo parities.

    Assuming local stereo parities, the parities can simply be reversed.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    """
    atm_par_dct = atom_stereo_parities(gra)
    atm_par_dct = dict_.transform_values(
        atm_par_dct, lambda x: x if x is None else not x
    )
    gra = set_atom_stereo_parities(gra, atm_par_dct)
    return gra


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
            backbone_only=False,
            par_eval_=parity_evaluator_flip_local_(),
            can_par_eval_=parity_evaluator_read_canonical_(),
            pri_dct=pri_dct_,
        )
    else:
        loc_gra = can_gra

    return loc_gra


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
            backbone_only=False,
            par_eval_=parity_evaluator_flip_local_(),
            pri_dct=pri_dct_,
        )
    else:
        can_gra = loc_gra

    return can_gra


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
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(atm_keys))}
    )

    par_eval_ = parity_evaluator_from_geometry_(
        geo, local_stereo=local_stereo, geo_idx_dct=geo_idx_dct
    )
    can_par_eval_ = None

    if local_stereo:
        # If requesting local stereo, we need an auxiliary canonical parity evaluator
        can_par_eval_ = parity_evaluator_from_geometry_(geo, geo_idx_dct=geo_idx_dct)

    gra, *_ = calculate_stereo(gra, par_eval_=par_eval_, can_par_eval_=can_par_eval_)

    return gra


# # core algorithm functions
def canonical_priorities(gra, backbone_only=True, pri_dct=None):
    """Determine canonical priorities for this graph's atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: dict[int: int]
    :returns: A dictionary of canonical priorities by atom key.
    :rtype: dict[int: int]
    """
    *_, pri_dct, _ = calculate_stereo(gra, backbone_only=backbone_only, pri_dct=pri_dct)
    return pri_dct


def calculate_stereo(
    gra,
    par_eval_: Optional[ParityEvaluator] = None,
    can_par_eval_: Optional[ParityEvaluator] = None,
    pri_dct: Optional[Dict[int, int]] = None,
    backbone_only: bool = True,
):
    """Algorithm for calculating stereo parities and priorities

    Depending on the parity evaluators passed in, this function can be used to:
        1. Assign canonical parities to a graph from a geometry
        2. Assign local parities to a graph from a geometry
        3. Convert canonical parities to local parities
        4. Convert local parities to canonical parities
        5. Determine reagent stereochemistry from a TS graph

    It can also be used simply to calculate canonical stereo priorities, which are
    always returned along with the stereo-assigned graph.

    :param gra: a molecular graph
    :type gra: automol graph data structure
    :param par_eval_: A parity evaluator, used for calculating priorities
    :type par_eval_: Optional[ParityEvaluator]
    :param can_par_eval_: When `par_eval_` is non-canonical, this parity
        evaluator must return canonical parities (used in priority calculation)
    :type can_par_eval_: Optional[ParityEvaluator]
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: Optional[Dict[int, int]]
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, a canonical priority mapping, and a flag indicating
        whether a TS graph has a canonical direction (None for non-TS graphs)
    :rtype: automol graph data structure, -//-, Dict[int, int], Optional[bool]
    """
    if is_ts_graph(gra):
        gra, can_gra, pri_dct, is_can_dir = _calculate_ts_stereo(
            gra, par_eval_, can_par_eval_, pri_dct=pri_dct
        )
    else:
        gra, can_gra, pri_dct = _calculate_stereo_core(
            gra, par_eval_, can_par_eval_, pri_dct=pri_dct
        )
        is_can_dir = None

    # 3. If requested, add in priorities for explicit hydrogens.
    if not backbone_only:
        pri_dct = reassign_hydrogen_priorities(gra, pri_dct)

    return gra, can_gra, pri_dct, is_can_dir


def _calculate_ts_stereo(
    tsg,
    par_eval_: Optional[ParityEvaluator] = None,
    can_par_eval_: Optional[ParityEvaluator] = None,
    pri_dct: Optional[Dict[int, int]] = None,
):
    """Algorithm for calculating stereo parities and priorities for a TS graph

    :param gra: a molecular graph
    :type gra: automol graph data structure
    :param par_eval_: A parity evaluator, used for calculating priorities
    :type par_eval_: Optional[ParityEvaluator]
    :param can_par_eval_: When `par_eval_` is non-canonical, this parity
        evaluator must return canonical parities (used in priority calculation)
    :type can_par_eval_: Optional[ParityEvaluator]
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: Optional[Dict[int, int]]
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, a canonical priority mapping, and a flag indicating
        whether a TS graph has a canonical direction
    :rtype: automol graph data structure, -//-, Dict[int, int], bool
    """
    tsg0 = tsg
    pri_dct0 = pri_dct

    # 1. Run the core stereo calculation algorithm in the forward TS direction
    tsg, can_tsg, pri_dct = _calculate_stereo_core(
        tsg0, par_eval_, can_par_eval_, pri_dct=pri_dct0
    )

    # 2. Rerun the core stereo calculation algorithm in the reverse TS direction
    rtsg = ts_reverse(tsg0)
    rtsg, rcan_tsg, rpri_dct = _calculate_stereo_core(
        rtsg, par_eval_, can_par_eval_, pri_dct=pri_dct0, is_rev_ts=True
    )

    # 3. Determine which direction is canonical
    is_can_dir = is_canonical_direction(can_tsg, pri_dct, rcan_tsg, rpri_dct)
    if is_can_dir is False:
        tsg = ts_reverse(rtsg)
        can_tsg = ts_reverse(rcan_tsg)
        pri_dct = rpri_dct

    return tsg, can_tsg, pri_dct, is_can_dir


def _calculate_stereo_core(
    gra,
    par_eval_: Optional[ParityEvaluator] = None,
    can_par_eval_: Optional[ParityEvaluator] = None,
    pri_dct: Optional[Dict[int, int]] = None,
    is_rev_ts: bool = False,
):
    """Core algorithm for calculating stereo parities and priorities

    :param gra: a molecular graph
    :type gra: automol graph data structure
    :param par_eval_: A parity evaluator, used for calculating priorities
    :type par_eval_: Optional[ParityEvaluator]
    :param can_par_eval_: When `par_eval_` is non-canonical, this parity
        evaluator must return canonical parities (used in priority calculation)
    :type can_par_eval_: Optional[ParityEvaluator]
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: Optional[Dict[int, int]]
    :param is_rev_ts: A flag passed to the parity evaluators, indicating that the graph
        is a reversed TS; defaults to False
    :type is_rev_ts: bool, optional
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, and a canonical priority mapping
    :rtype: automol graph data structure, -//-, Dict[int, int]
    """
    gra0 = gra
    par_eval_ = parity_evaluator_read_canonical_() if par_eval_ is None else par_eval_
    can_par_eval_ = par_eval_ if can_par_eval_ is None else can_par_eval_

    # Graph 1 will be for the priority calculation, graph 2 for the parity
    # assignments that will be returned.
    gra = can_gra = without_stereo(gra0)

    pri_dct0 = 0  # Can't use None, since pri_dct can be None
    while pri_dct != pri_dct0:
        pri_dct0 = pri_dct

        # a. Refine priorities based on the canonical graph
        pri_dct = refine_priorities(can_gra, pri_dct)

        # b. Find stereogenic atoms and bonds based on current priorities
        keys = stereocenter_keys_from_priorities(can_gra, pri_dct, new=True)

        # c. If there are none, the calculation is complete. Exit the loop.
        if not keys:
            break

        # d. Assign parities to the canonical graph using the canonical parity evaluator
        can_par_dct = can_par_eval_(gra0, pri_dct, keys, is_rev_ts=is_rev_ts)
        can_gra = set_stereo_parities(can_gra, can_par_dct)

        # e. Assign parities to the auxiliary graph using the auxiliary parity evaluator
        if par_eval_ != can_par_eval_:
            par_dct = par_eval_(gra0, pri_dct, keys, is_rev_ts=is_rev_ts)
            gra = set_stereo_parities(gra, par_dct)
        else:
            gra = can_gra

    return gra, can_gra, pri_dct


def refine_priorities(gra, pri_dct=None, _backbone_only=True):
    """Refine the canonical priorities for this graph based on some sort value

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param srt_eval_: An evaluator for sort values, based on current
        priorities. Curried such that srt_val_(pri_dct)(key) returns the
        sort value.
    """
    gra = without_dummy_atoms(gra)
    cgras = connected_components(gra)
    cpri_dcts = [
        None if pri_dct is None else dict_.by_key(pri_dct, ks, fill=False)
        for ks in map(atom_keys, cgras)
    ]
    pri_dct = {}
    for cgra, cpri_dct in zip(cgras, cpri_dcts):
        cpri_dct = _refine_priorities(
            cgra, pri_dct=cpri_dct, _backbone_only=_backbone_only
        )
        pri_dct.update(cpri_dct)

    return pri_dct


def break_priority_ties(gra, pri_dct):
    """Break ties within priority classes.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param srt_eval_: An evaluator for sort values, based on current class
        indices. Curried such that srt_val_(pri_dct)(key) returns the sort
        value.
    """
    gra = without_dummy_atoms(gra)
    cgras = connected_components(gra)
    cpri_dcts = [
        None if pri_dct is None else dict_.by_key(pri_dct, ks, fill=False)
        for ks in map(atom_keys, cgras)
    ]
    pri_dct = {}
    for cgra, cpri_dct in zip(cgras, cpri_dcts):
        cpri_dct = _break_priority_ties(cgra, pri_dct=cpri_dct)
        pri_dct.update(cpri_dct)

    return pri_dct


def reassign_hydrogen_priorities(gra, pri_dct, neg=False):
    """Reassign priorities to hydrogens, negating them on request

    :param neg: Negate hydrogen keys, to give them lowest priority?
    :param neg: bool
    """
    gra = without_dummy_atoms(gra)
    cgras = connected_components(gra)
    cpri_dcts = [
        None if pri_dct is None else dict_.by_key(pri_dct, ks, fill=False)
        for ks in map(atom_keys, cgras)
    ]
    pri_dct = {}
    for cgra, cpri_dct in zip(cgras, cpri_dcts):
        cpri_dct = _reassign_hydrogen_priorities(cgra, pri_dct=cpri_dct, neg=neg)
        pri_dct.update(cpri_dct)

    return pri_dct


def _refine_priorities(gra, pri_dct: Optional[dict] = None, _backbone_only=True):
    """Refine the canonical priorities for this graph based on some sort value

    (Only for connected graphs)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param srt_eval_: An evaluator for sort values, based on current
        priorities. Curried such that srt_val_(pri_dct)(key) returns the
        sort value.
    """
    assert gra == without_dummy_atoms(
        gra
    ), f"Remove dummy atoms:\n{gra}\n{without_dummy_atoms(gra)}"

    if _backbone_only:
        gra = implicit(gra)

    pri_dct = (
        dict_.by_key({}, atom_keys(gra), fill_val=0) if pri_dct is None else pri_dct
    )
    next_pri = max(pri_dct.values()) + 1
    pri_dct = dict_.by_key(pri_dct, atom_keys(gra), fill_val=next_pri)
    srt_eval_ = sort_evaluator_atom_invariants_(gra)

    ngb_keys_dct = atoms_neighbor_atom_keys(gra)

    cla_dct = class_dict_from_priority_dict(pri_dct)

    # Set up the new_clas list, containing priorities and priority classes that
    # are up for re-evaluation.
    new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
    new_clas = sorted(new_cla_dct.items(), reverse=True)

    while new_clas:
        # Pop the next priority class for re-evaluation
        idx, cla = new_clas.pop(0)

        # Sort and partition the priority class based on sort values. After the
        # first iteration, only the neighboring priority classes cause further
        # subdivision.
        srt_val_ = srt_eval_(pri_dct)
        cla = sorted(cla, key=srt_val_)
        parts = [tuple(ks) for _, ks in itertools.groupby(cla, srt_val_)]

        # Assign new priorities to the class partitions and update pri_dct and
        # cla_dct.
        new_idx = idx
        if len(parts) > 1:
            new_idxs = []
            for part in parts:
                cla_dct[new_idx] = part
                pri_dct.update({k: new_idx for k in part})

                # Track the set of new indices
                new_idxs.append(new_idx)

                # Increment the index for the next class by the number of
                # members in this one, so that that priorities are stable.
                new_idx += len(part)

            # Identify indices of classes with neighboring atoms, as these may
            # be affected by the re-classification.
            ngb_idxs = frozenset()
            for new_idx in new_idxs:
                new_cla = cla_dct[new_idx]

                # Get neighboring keys to the members of this new class.
                ngb_keys = frozenset.union(*map(ngb_keys_dct.__getitem__, new_cla))

                # Get priorities of these neighboring atoms.
                ngb_idxs |= frozenset(map(pri_dct.get, ngb_keys))

                # Don't include classes that were already up for re-evaluation.
                ngb_idxs -= frozenset(dict(new_clas))

            for ngb_idx in sorted(ngb_idxs):
                ngb_cla = cla_dct[ngb_idx]
                if len(ngb_cla) > 1:
                    new_clas.insert(0, (ngb_idx, cla_dct[ngb_idx]))

    return pri_dct


def _break_priority_ties(gra, pri_dct):
    """Break ties within priority classes.

    (Only for connected graphs)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param srt_eval_: An evaluator for sort values, based on current class
        indices. Curried such that srt_val_(pri_dct)(key) returns the sort
        value.
    """
    assert gra == without_dummy_atoms(
        gra
    ), f"Remove dummy atoms:\n{gra}\n{without_dummy_atoms(gra)}"

    pri_dct = pri_dct.copy()

    cla_dct = class_dict_from_priority_dict(pri_dct)

    # Set up the new_clas list, containing priorities and priority classes that
    # are up for re-evaluation.
    new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
    new_clas = sorted(new_cla_dct.items(), reverse=True)

    while new_clas:
        # Pop the next priority class for tie-breaking
        idx, cla = new_clas.pop(0)
        cla = list(cla)

        # Give the last element of this priority class a new index
        new_idx = idx + len(cla) - 1
        pri_dct[cla[-1]] = new_idx

        # Now, refine partitions based on the change just made.
        pri_dct = refine_priorities(gra, pri_dct, _backbone_only=False)

        # Update the list of classes needing further tie breaking
        cla_dct = class_dict_from_priority_dict(pri_dct)
        new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
        new_clas = sorted(new_cla_dct.items(), reverse=True)

    return pri_dct


def _reassign_hydrogen_priorities(gra, pri_dct, neg=False):
    """Reassign priorities to hydrogens, negating them on request

    (Only for connected graphs)

    :param neg: Negate hydrogen keys, to give them lowest priority?
    :param neg: bool
    """
    pri_dct = dict_.by_key(pri_dct, backbone_keys(gra), fill=False)

    # Assigne all hydrogens to the same priority class
    next_pri = max(pri_dct.values()) + 1
    pri_dct.update({k: next_pri for k in nonbackbone_hydrogen_keys(gra)})
    pri_dct = _refine_priorities(gra, pri_dct, _backbone_only=False)

    # If requested, negate priorities for *all* hydrogen keys, for consistency
    if neg:
        all_hyd_keys = atom_keys(gra, symb="H")
        pri_dct = {k: -abs(p) if k in all_hyd_keys else p for k, p in pri_dct.items()}

    return pri_dct


def sort_evaluator_atom_invariants_(gra):
    """A sort function based on atom invariants with two levels of currying.

    To get the sort value for a specific key, use
        srt_val = sort_evaluator_atom_invariants_(gra)(pri_dct)(key)

    My reasoning for doing things this way is that `gra` never changes, but
    `pri_dct` does, so we need to be able to update the function with new
    index dictionaries. Ultimately, it is convenient to return a function
    of `key` only because this can be passed to standard python sorting and
    grouping functions such as `sorted()` and `itertools.groupby()`.

    The canonical order has been modified relative to the one used by
    Scheider (2015) to be more InChI-like (following Hill ordering).

    :param gra: molecular graph
    :type gra: automol graph data structure
    """

    def _replace_none(val):
        return numpy.inf if val is None else val

    def _replace_reacting_bond_order(val):
        """Ensure that reacting bonds will sort after non-reacting bonds:
        non-reacting bonds < forming bonds < breaking bonds
        """
        assert val == round(val, 1)
        return int(val * 100) if val in (0.1, 0.9) else val

    def _normalize_symbol(symb):
        """normalize atomic symbols to make them sort as follows:
        (C first, others in alphabetical order, H last)
        (Hydrogens should only appear for TS graphs)
        """
        symb = ptab.to_symbol(symb)
        if symb == "C":
            symb = ""  # Sorts first (always)
        if symb == "H":
            symb = "zz"  # Sorts last against atomic symbols
        return symb

    # For the main properties, remove reacting bonds and/or hydrogens before
    # calculating
    imp_gra = without_pi_bonds(implicit(gra))
    imp_gra_no_rxbs = without_bonds_by_orders(imp_gra, [0.1, 0.9])
    bnds_dct_no_rxbs = atoms_bond_keys(imp_gra_no_rxbs)

    symb_dct = dict_.transform_values(atom_symbols(gra), _normalize_symbol)
    hnum_dct = atom_implicit_hydrogens(imp_gra_no_rxbs)
    mnum_dct = mass_numbers(imp_gra_no_rxbs)
    apar_dct = dict_.transform_values(
        atom_stereo_parities(imp_gra_no_rxbs), _replace_none
    )

    bpar_dct = bond_stereo_parities(imp_gra_no_rxbs)

    def _sortable_bond_stereo_values(bnd_keys):
        bpars = dict_.values_by_key(bpar_dct, bnd_keys)
        bpars = tuple(sorted(map(_replace_none, bpars)))
        return bpars

    bpars_dct = dict_.transform_values(bnds_dct_no_rxbs, _sortable_bond_stereo_values)

    # For bond orders, *KEEP* reacting bonds
    bord_dct = bond_orders(imp_gra)

    def _sortable_bond_orders(bnd_keys):
        bords = dict_.values_by_key(bord_dct, bnd_keys)
        bords = tuple(sorted(map(_replace_reacting_bond_order, bords)))
        return bords

    bords_dct = dict_.transform_values(
        atoms_bond_keys(imp_gra, ts_=True), _sortable_bond_orders
    )

    # For neighboring keys, *KEEP* reacting bonds and hydrogens
    nkeys_dct = atoms_neighbor_atom_keys(imp_gra, ts_=True)

    # Fill in missing hydrogen values
    hyd_keys = nonbackbone_hydrogen_keys(gra)
    symb_dct.update({k: "zzz" for k in hyd_keys})
    bnds_dct_no_rxbs.update({k: [] for k in hyd_keys})
    hnum_dct.update({k: -1 for k in hyd_keys})
    mnum_dct.update({k: -1 for k in hyd_keys})
    apar_dct.update({k: -1 for k in hyd_keys})
    bpars_dct.update({k: [-1] for k in hyd_keys})
    bords_dct.update({k: [-1] for k in hyd_keys})
    # Use the actual neighboring atoms of hydrogens, so that they will sort
    # based on their parent atoms
    nkeys_dct.update(dict_.by_key(atoms_neighbor_atom_keys(gra, ts_=True), hyd_keys))

    def _evaluator(pri_dct: dict):
        """Sort value evaluator based on current priorities.

        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        """

        def _value(key):
            symb = symb_dct[key]  # symbol
            deg = len(bnds_dct_no_rxbs[key])  # number of bonds
            hnum = hnum_dct[key]  # number of hydrogens
            mnum = mnum_dct[key]
            apar = apar_dct[key]
            bpars = bpars_dct[key]
            bords = bords_dct[key]
            nidxs = tuple(sorted(map(pri_dct.get, nkeys_dct[key])))
            return (symb, deg, hnum, mnum, apar, bpars, bords, nidxs)

        return _value

    return _evaluator


# # parity evaluators
def parity_evaluator_from_geometry_(
    geo, local_stereo=False, geo_idx_dct=None
) -> ParityEvaluator:
    r"""Determines stereo parity from a geometry

    (For internal use by the calculate_priorities_and_assign_stereo() function)

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :returns: A parity evaluator, which takes a graph, a priority mapping, and a set of
        keys and returns a dictionary of parities for those keys
    :rtype: ParityEvaluator
    """
    geo_idx_dct_ = geo_idx_dct

    def _evaluator(
        gra, pri_dct: Dict[int, int], keys: List[int], is_rev_ts: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param is_rev_ts: Is this a reversed TS graph?
        :type is_rev_ts: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """
        atm_keys = sorted(atom_keys(gra))
        geo_idx_dct = (
            geo_idx_dct_
            if geo_idx_dct_ is not None
            else {k: i for i, k in enumerate(sorted(atm_keys))}
        )

        assert gra == explicit(
            gra
        ), "Explicit graph should be used when evaluating from geometry."

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                nkeys = atom_stereo_sorted_neighbor_keys(gra, key, pri_dct=pri_dct)

                # Get the atom parity
                par = geometry_atom_parity(
                    gra, geo, key, nkeys, geo_idx_dct=geo_idx_dct
                )

            # Otherwise, this is a bond
            else:
                assert (
                    isinstance(key, abc.Collection) and len(key) == 2
                ), f"{key} is not a valid bond key."
                key1, key2 = key
                nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(
                    gra, key1, key2, pri_dct=pri_dct
                )

                # Get the bond parity
                par = geometry_bond_parity(
                    gra, geo, (key1, key2), (nkey1s, nkey2s), geo_idx_dct=geo_idx_dct
                )

            return par

        par_dct = {k: _parity(k) for k in keys}
        if local_stereo:
            can_gra = set_stereo_parities(gra, par_dct)

            loc_par_eval_ = parity_evaluator_flip_local_()
            par_dct = loc_par_eval_(can_gra, pri_dct, keys, is_rev_ts)

        return par_dct

    return _evaluator


def parity_evaluator_read_canonical_() -> ParityEvaluator:
    """Determines stereo parity from a graph with canonical stereo

    (For internal use by the calculate_priorities_and_assign_stereo() function)

    This is the trivial case, where parities in the graph are already
    assumed to be canonical and so nothing needs to be done to calculate
    them.

    :returns: A parity evaluator, which takes a graph, a priority mapping, and a set of
        keys and returns a dictionary of parities for those keys
    :rtype: ParityEvaluator
    """

    def _evaluator(
        gra, pri_dct: Dict[int, int], keys: List[int], is_rev_ts: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param is_rev_ts: Is this a reversed TS graph?
        :type is_rev_ts: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """
        # Do-nothing lines to prevent linting complaint
        assert pri_dct or not pri_dct
        assert is_rev_ts or not is_rev_ts

        return dict_.by_key(stereo_parities(gra), keys)

    return _evaluator


def parity_evaluator_flip_local_() -> ParityEvaluator:
    """Determines canonical from local stereo parity or vice versa (same
    operation)

    (For internal use by the calculate_priorities_and_assign_stereo() function)

    :returns: A parity evaluator, which takes a graph, a priority mapping, and a set of
        keys and returns a dictionary of parities for those keys
    :rtype: ParityEvaluator
    """

    def _evaluator(
        gra, pri_dct: Dict[int, int], keys: List[int], is_rev_ts: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param is_rev_ts: Is this a reversed TS graph?
        :type is_rev_ts: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """
        gra = explicit(gra)
        par_dct = stereo_parities(gra)

        flip_dct = sn2_local_stereo_reversal_flips(gra) if is_rev_ts else {}

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                par = par_dct[key]

                loc_srt = atom_stereo_sorted_neighbor_keys(gra, key)
                can_srt = atom_stereo_sorted_neighbor_keys(gra, key, pri_dct=pri_dct)

                is_odd = util.is_odd_permutation(loc_srt, can_srt)
                ret_par = is_odd ^ par

                if key in flip_dct:
                    ret_par ^= flip_dct[key]

            # Otherwise, this is a bond
            else:
                assert (
                    isinstance(key, abc.Collection) and len(key) == 2
                ), f"{key} is not a valid bond key."
                par = par_dct[key]

                key1, key2 = key
                loc_srt1, loc_srt2 = bond_stereo_sorted_neighbor_keys(gra, key1, key2)
                can_srt1, can_srt2 = bond_stereo_sorted_neighbor_keys(
                    gra, key1, key2, pri_dct=pri_dct
                )

                if not (loc_srt1[-1] == can_srt1[-1]) ^ (loc_srt2[-1] == can_srt2[-1]):
                    ret_par = par
                else:
                    ret_par = not par

            return ret_par

        return {k: _parity(k) for k in keys}

    return _evaluator


def parity_evaluator_reagents_from_ts_(tsg, prod=False) -> ParityEvaluator:
    r"""Determines reactant or product stereochemistry from a TS graph

    (For internal use by the calculate_priorities_and_assign_stereo() function)

    Sn2 constraint is taken care of by parity_evaluator_flip_local_

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param prod: Do this for the products, instead of the reactants?
    :type prod: bool
    :returns: A parity evaluator, which takes a graph, a priority mapping, and a set of
        keys and returns a dictionary of parities for those keys
    :rtype: ParityEvaluator
    """
    tsg = ts_reverse(tsg) if prod else tsg

    # Sn2 reactions are handled by converting to local stereochemistry
    loc_tsg = to_local_stereo(tsg)

    # Handle constrained insertion/eliminations
    cpar_dct = constrained_1_2_insertion_local_parities(loc_tsg)
    loc_tsg = set_stereo_parities(loc_tsg, cpar_dct)

    # Handle vinyl radical additions
    vpar_dct = vinyl_addition_local_parities(loc_tsg)
    loc_tsg = set_stereo_parities(loc_tsg, vpar_dct)

    # Generate the reagents graph with local parities
    loc_rgra = ts_reactants_graph_without_stereo(loc_tsg, keep_stereo=True)

    # Now that we have handled the exceptions, the local parities correspond to
    # what they will be for the reactants/products graph, so we can simply flip
    # the local stereo to find canonical assignments
    par_eval_flip_ = parity_evaluator_flip_local_()

    def _evaluator(
        gra, pri_dct: Dict[int, int], keys: List[int], is_rev_ts: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param is_rev_ts: Is this a reversed TS graph?
        :type is_rev_ts: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """

        # Do-nothing line to prevent linting complaint
        assert gra or not gra
        assert is_rev_ts or not is_rev_ts
        return par_eval_flip_(loc_rgra, pri_dct, keys)

    return _evaluator


# # core algorithm helpers
def stereocenter_keys_from_priorities(
    gra, pri_dct: Dict[int, int], atom: bool = True, bond: bool = True, new: bool = True
):
    """Find stereogenic atoms and bonds in this graph, given a set of atom
    priority values

    If `new` flag is set to `True`, only unassigned stereocenters will be detected

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: Dict[int, int]
    :param atom: Include atom stereocenters? defaults to True
    :type atom: bool, optional
    :param bond: Include bond stereocenters? defaults to True
    :type bond: bool, optional
    :param new: Detect only new, unassigned stereocenters?
    :type new: bool, optional
    :returns: the stereogenic atom and bond keys
    :rtype: frozenset
    """

    def _atom_is_stereogenic(key):
        nkeys = atom_stereo_sorted_neighbor_keys(gra, key, pri_dct=pri_dct)
        assert len(nkeys) <= 4, f"Too many neighbors! {nkeys}"
        pris = list(map(pri_dct.get, nkeys))
        return len(set(pris)) == len(pris)

    def _bond_is_stereogenic(key):
        nkeys_pair = bond_stereo_sorted_neighbor_keys(gra, *key, pri_dct=pri_dct)
        ret = True
        for nkeys in nkeys_pair:
            assert len(nkeys) <= 2, f"Too many neighbors! {nkeys}"
            pris = list(map(pri_dct.get, nkeys))
            ret &= len(nkeys) == 1 or len(set(pris)) == len(pris)
        return ret

    def _is_stereogenic(key):
        return (
            _atom_is_stereogenic(key)
            if isinstance(key, numbers.Number)
            else _bond_is_stereogenic(key)
        )

    keys = stereocenter_candidate_keys(gra, atom=atom, bond=bond)
    if new:
        keys -= stereo_keys(gra)

    ste_keys = frozenset(filter(_is_stereogenic, keys))
    return ste_keys


def class_dict_from_priority_dict(pri_dct: dict):
    """Obtain a class dictionary from a priority dictionary.

    :param pri_dct: A dictionary mapping atom keys to priorities.
    :type pri_dct: dict
    :returns: A dictionary mapping priorities onto the full set of keys
        for that priority class, as a sorted tuple.
    :rtype: dict[int: tuple]
    """
    keys = sorted(pri_dct.keys())
    clas = sorted(keys, key=pri_dct.get)
    cla_dct = {i: tuple(c) for i, c in itertools.groupby(clas, key=pri_dct.get)}
    return cla_dct
