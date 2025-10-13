"""canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""

import itertools
from collections.abc import Callable
from typing import Any, Dict, List, Optional, Tuple

import numpy
from phydat import ptab

from ...util import dict_
from ._00core import (
    atom_bond_counts,
    atom_implicit_hydrogens,
    atom_keys,
    atom_stereo_keys,
    atom_stereo_parities,
    atom_symbols,
    atoms_bond_keys,
    atoms_neighbor_atom_keys,
    backbone_keys,
    bond_orders,
    bond_stereo_keys,
    bond_stereo_parities,
    explicit,
    has_stereo,
    implicit,
    invert_atom_stereo_parities,
    is_ts_graph,
    mass_numbers,
    nonbackbone_hydrogen_keys,
    set_stereo_parities,
    stereo_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reverse,
    with_explicit_stereo_hydrogens,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_pi_bonds,
    without_stereo,
)
from ._00core import (
    relabel as relabel_,
)
from ._02algo import connected_components, is_connected
from ._03kekule import bad_stereo_bond_keys_from_kekule, kekule
from ._05stereo import (
    CenterNeighborDict,
    parity_evaluator_flip_from_graph,
    parity_evaluator_read_from_graph,
    stereocenter_candidates,
    stereocenter_candidates_grouped,
    unassigned_stereocenter_keys_from_candidates,
)

# Special types
ParityEvaluator = Callable[[Any, Dict[int, int], List[int], bool], Dict[int, int]]


# # canonical key functions
def canonical(gra):
    """A graph relabeled with canonical keys

    Stereo parities in the graph are assumed to be canonical.

    Requires a connected graph

    :param gra: a connected molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :returns: a new molecular graph with canonical keys; if explicit
        hydrogens are included, they will be relabeled as well
    """
    can_key_dct = canonical_keys(gra)
    return relabel_(gra, can_key_dct)


def canonical_priorities(gra, pri_dct=None):
    """Determine canonical priorities for this graph's atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: dict[int: int]
    :returns: A dictionary of canonical priorities by atom key.
    :rtype: dict[int: int]
    """
    *_, pri_dct, _ = calculate_stereo(gra, pri_dct=pri_dct)
    return pri_dct


def canonical_keys(gra):
    """Determine canonical keys for this graph.

    Stereo parities in the graph are assumed to be canonical.

    Requires a connected graph

    :param gra: a connected molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :returns: a dictionary of canonical keys by atom key
    :rtype: dict[int: int]
    """
    assert is_connected(gra), f"Requires a connected graph.\n{gra}"
    atm_par_dct0 = atom_stereo_parities(gra)
    bnd_par_dct0 = bond_stereo_parities(gra)

    gra, _, pri_dct, *_ = calculate_stereo(gra)
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


def canonical_amchi_graph_with_numbers(
    gra, stereo: bool = True
) -> Tuple[Any, Dict[int, int], bool, bool]:
    """Put a connected graph in AMChI-canonical form

    :param gra: A connected graph
    :type gra: automol graph data structure
    :param stereo: Include stereo, if present? defaults to True
    :type stereo: bool, optional
    """
    assert is_connected(gra), f"Graph is disconnected\n{gra}"
    gra = without_dummy_atoms(gra)
    gra = gra if stereo else without_stereo(gra)
    ste_keys = stereo_keys(gra)
    atm_ste_keys = atom_stereo_keys(gra)
    cand_dct = stereocenter_candidates(gra)

    # 1. Convert to implicit-hydrogen form
    gra = implicit(gra)

    # 2. Determine the canonical TS direction
    *_, pri_dct, is_can_dir = calculate_stereo(gra, cand_dct=cand_dct)
    gra = ts_reverse(gra) if is_can_dir is False else gra
    cand_dct.update(stereocenter_candidates(gra, bond=False))

    # 3. Determine the canonical enantiomer
    is_can_enant = None
    if len(atm_ste_keys) == 1:
        rgra = invert_atom_stereo_parities(gra)
        is_can_enant = is_canonical_enantiomer(
            gra, pri_dct, rgra, pri_dct, cand_dct=cand_dct
        )
        gra = rgra if is_can_enant is False else gra
    elif len(atm_ste_keys) > 1:
        par_eval_ = parity_evaluator_flip_from_graph
        loc_gra = set_stereo_parities(gra, par_eval_(gra, ste_keys, pri_dct, cand_dct))
        rloc_gra = invert_atom_stereo_parities(loc_gra)
        rgra, _, rpri_dct, _ = calculate_stereo(rloc_gra, par_eval_=par_eval_)

        is_can_enant = is_canonical_enantiomer(
            gra, pri_dct, rgra, rpri_dct, cand_dct=cand_dct
        )

        gra = rgra if is_can_enant is False else gra
        pri_dct = rpri_dct if is_can_enant is False else pri_dct

    # 4. Find canonical numbers and relabel
    num_dct = break_priority_ties(gra, pri_dct, neg_hkeys=False)
    gra = relabel_(gra, num_dct)

    return gra, num_dct, is_can_dir, is_can_enant


def smiles_graph(
    gra: Any, res_stereo: bool = True, exp: bool = False, dummy: bool = False
) -> Any:
    """Put a connected graph in a form appropriate for writing SMILES strings.

    :param gra: molecular graph
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :param res_stereo: allow resonant double-bond stereo?
    :param dummy: Retain dummy atoms? Defaults to False
    :return: The graph, kekulized with appropriate implicit/explicit hydrogens
    """
    # Find a dominant resonance
    kgr = kekule(gra, max_stereo_overlap=True)
    return smiles_graph_from_kekule(kgr, res_stereo=res_stereo, exp=exp, dummy=dummy)


def smiles_graph_from_kekule(
    kgr: Any, res_stereo: bool = True, exp: bool = False, dummy: bool = False
) -> Any:
    """Put a connected graph in a form appropriate for writing SMILES strings.

    :param gra: molecular graph
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :param res_stereo: allow resonant double-bond stereo?
    :param dummy: Retain dummy atoms? Defaults to False
    :return: The graph, kekulized with appropriate implicit/explicit hydrogens
    """
    if not dummy:
        kgr = without_dummy_atoms(kgr)

    # Add in missing explicit hydrogens (connected to H or needed for stereo)
    kgr = explicit(kgr, atm_keys=atom_keys(kgr, symb="H"), neg=True)
    kgr = with_explicit_stereo_hydrogens(kgr, all_=False, neg=True)

    # Make all other hydrogens implicit
    if not exp:
        bbn_keys = backbone_keys(kgr, hyd=False)
        ste_keys = set(itertools.chain(*bond_stereo_keys(kgr)))
        nbnd_dct = atom_bond_counts(kgr, bond_order=False, with_implicit=True)
        imp_keys = {k for k in bbn_keys if not (k in ste_keys and nbnd_dct[k] <= 2)}
        kgr = implicit(kgr, atm_keys=imp_keys)

    # Remove stereo parities at single bonds if requested
    if not res_stereo:
        bad_bkeys = bad_stereo_bond_keys_from_kekule(kgr)
        kgr = without_stereo(kgr, bnd_keys=bad_bkeys)

    return kgr


def stereo_assignment_representation(gra, pri_dct: dict, cand_dct: dict | None = None):
    """Generate a representation of a stereo assignment, for checking for
    symmetric equivalence or for determining a canonical enantiomer.

    :param gra: molecular graph
    :param pri_dct: A dictionary mapping atom keys to priorities
    :param cand_dct: A dictionary mapping stereocenter candidates to their neighbor keys
    :returns: A canonical representation of the assignment
    """
    cand_dct = stereocenter_candidates(gra) if cand_dct is None else cand_dct

    atm_keys = sorted(atom_stereo_keys(gra), key=pri_dct.get)
    bnd_keys = sorted(bond_stereo_keys(gra), key=lambda x: sorted(map(pri_dct.get, x)))

    atm_nkeys_dct, bnd_nkeys_dct = stereocenter_candidates_grouped(cand_dct, pri_dct)

    atm_pars = dict_.values_by_key(atom_stereo_parities(gra), atm_keys)
    bnd_pars = dict_.values_by_key(bond_stereo_parities(gra), bnd_keys)
    atm_pris = tuple(map(pri_dct.get, atm_keys))
    bnd_pris = tuple(map(pri_dct.get, bnd_keys))
    atm_npris = tuple(tuple(map(pri_dct.get, atm_nkeys_dct[k])) for k in atm_keys)
    bnd_npris = tuple(tuple(map(pri_dct.get, bnd_nkeys_dct[k])) for k in bnd_keys)

    rep = (atm_pars, bnd_pars, atm_pris, bnd_pris, atm_npris, bnd_npris)
    return rep


def ts_direction_representation(tsg, pri_dct: dict, cand_dct: dict | None = None):
    """Generate a representation of the reaction direction to determine the
    canonical direction of a TS graph.

    The reaction representation consists of two pieces:

    1. Canonical representations of the forming and breaking bonds,
    respectively.
    2. Canonical representations of the stereo assignments.

    :param tsg: A TS graph
    :param pri_dct: A dictionary mapping atom keys to priorities
    :param cand_dct: A dictionary mapping stereocenter candidates to their neighbor keys
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
    # Rep value 3: Stereo assignment
    ste_rep = stereo_assignment_representation(tsg, pri_dct, cand_dct=cand_dct)
    rep = (rxn_rep1, rxn_rep2, ste_rep)
    return rep


def is_canonical_enantiomer(
    ugra, upri_dct, rgra, rpri_dct, cand_dct: dict | None = None
):
    """Determine whether this enantiomer is the canonical one.

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
    cand_dct = stereocenter_candidates(ugra) if cand_dct is None else cand_dct

    urep = stereo_assignment_representation(ugra, upri_dct, cand_dct=cand_dct)
    rrep = stereo_assignment_representation(rgra, rpri_dct, cand_dct=cand_dct)
    return True if (urep < rrep) else False if (urep > rrep) else None


def is_canonical_direction(
    ftsg,
    fpri_dct,
    rtsg,
    rpri_dct,
    fcand_dct: dict | None = None,
    rcand_dct: dict | None = None,
):
    """Determine whether this TS direction is the canonical one.

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
    frep = ts_direction_representation(ftsg, fpri_dct, cand_dct=fcand_dct)
    rrep = ts_direction_representation(rtsg, rpri_dct, cand_dct=rcand_dct)
    return True if (frep < rrep) else False if (frep > rrep) else None


# # core algorithm functions
def calculate_stereo(
    gra,
    par_eval_: Optional[ParityEvaluator] = None,
    can_par_eval_: Optional[ParityEvaluator] = None,
    pri_dct: Optional[Dict[int, int]] = None,
    cand_dct: Optional[CenterNeighborDict] = None,
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
    :param cand_dct: Optional pre-calculated stereocenter candidates
    :type cand_dct: Optional[CenterNeighborDict]
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, a canonical priority mapping, and a flag indicating
        whether a TS graph has a canonical direction (None for non-TS graphs)
    :rtype: automol graph data structure, -//-, Dict[int, int], Optional[bool]
    """
    if is_ts_graph(gra):
        gra, can_gra, pri_dct, is_can_dir = _calculate_ts_stereo(
            gra, par_eval_, can_par_eval_, pri_dct=pri_dct, cand_dct=cand_dct
        )
    else:
        gra, can_gra, pri_dct = _calculate_stereo_core(
            gra, par_eval_, can_par_eval_, pri_dct=pri_dct, cand_dct=cand_dct
        )
        is_can_dir = None

    return gra, can_gra, pri_dct, is_can_dir


def _calculate_ts_stereo(
    tsg,
    par_eval_: Optional[ParityEvaluator] = None,
    can_par_eval_: Optional[ParityEvaluator] = None,
    pri_dct: Optional[Dict[int, int]] = None,
    cand_dct: Optional[CenterNeighborDict] = None,
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
    :param cand_dct: Optional pre-calculated stereocenter candidates
    :type cand_dct: Optional[CenterNeighborDict]
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, a canonical priority mapping, and a flag indicating
        whether a TS graph has a canonical direction
    :rtype: automol graph data structure, -//-, Dict[int, int], bool
    """
    tsg0 = tsg
    pri_dct0 = pri_dct

    # 1. Run the core stereo calculation algorithm in the forward TS direction
    tsg, can_tsg, pri_dct = _calculate_stereo_core(
        tsg0, par_eval_, can_par_eval_, pri_dct=pri_dct0, cand_dct=cand_dct
    )

    # 2. Rerun the core stereo calculation algorithm in the reverse TS direction
    rtsg = ts_reverse(tsg0)
    rtsg, rcan_tsg, rpri_dct = _calculate_stereo_core(
        rtsg, par_eval_, can_par_eval_, pri_dct=pri_dct0, cand_dct=cand_dct, is_rev=True
    )

    # 3. Determine which direction is canonical
    is_can_dir = is_canonical_direction(
        can_tsg, pri_dct, rcan_tsg, rpri_dct, fcand_dct=cand_dct
    )
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
    cand_dct: Optional[CenterNeighborDict] = None,
    is_rev: bool = False,
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
    :param cand_dct: Optional pre-calculated stereocenter candidates
    :type cand_dct: Optional[CenterNeighborDict]
    :param is_rev_: A flag passed to the parity evaluators, indicating that the graph
        is a reversed TS; defaults to False
    :type is_rev: bool, optional
    :returns: A gaph with stereo assigned from `par_eval_`, a graph with stereo assigned
        from `can_par_eval_`, and a canonical priority mapping
    :rtype: automol graph data structure, -//-, Dict[int, int]
    """
    gra0 = gra
    par_eval_ = parity_evaluator_read_from_graph if par_eval_ is None else par_eval_
    can_par_eval_ = par_eval_ if can_par_eval_ is None else can_par_eval_
    cand_dct = stereocenter_candidates(gra) if cand_dct is None else cand_dct

    # Graph 1 will be for the priority calculation, graph 2 for the parity
    # assignments that will be returned.
    gra = can_gra = without_stereo(gra0)

    pri_dct0 = 0  # Can't use None, since pri_dct can be None
    while pri_dct != pri_dct0:
        pri_dct0 = pri_dct

        # a. Refine priorities based on the canonical graph
        pri_dct = refine_priorities(can_gra, pri_dct)

        # b. Find stereogenic atoms and bonds based on current priorities
        keys = unassigned_stereocenter_keys_from_candidates(can_gra, cand_dct, pri_dct)

        # c. If there are none, the calculation is complete. Exit the loop.
        if not keys:
            break

        # d. Assign parities to the canonical graph using the canonical parity evaluator
        can_par_dct = can_par_eval_(gra0, keys, pri_dct, cand_dct, is_rev_ts=is_rev)
        can_gra = set_stereo_parities(can_gra, can_par_dct)

        # e. Assign parities to the auxiliary graph using the auxiliary parity evaluator
        if par_eval_ != can_par_eval_:
            par_dct = par_eval_(gra0, keys, pri_dct, cand_dct, is_rev_ts=is_rev)
            gra = set_stereo_parities(gra, par_dct)
        else:
            gra = can_gra

    return gra, can_gra, pri_dct


def refine_priorities(gra, pri_dct: Optional[Dict[int, int]] = None):
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
    gras = connected_components(gra)
    pri_dcts = [
        None if pri_dct is None else dict_.by_key(pri_dct, ks, fill=False)
        for ks in map(atom_keys, gras)
    ]
    pri_dct = {}
    for gra_, pri_dct_ in zip(gras, pri_dcts):
        pri_dct_ = _refine_priorities(gra_, pri_dct=pri_dct_)
        pri_dct.update(pri_dct_)

    return pri_dct


def break_priority_ties(gra, pri_dct: Dict[int, int], neg_hkeys: bool = False):
    """Break ties within priority classes.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param neg_hkeys: Negate hydrogen keys? defaults to False
    :type neg_hkeys: bool, optional
    :param srt_eval_: An evaluator for sort values, based on current class
        indices. Curried such that srt_val_(pri_dct)(key) returns the sort
        value.
    """
    gra = without_dummy_atoms(gra)
    gras = connected_components(gra)
    pri_dcts = [
        None if pri_dct is None else dict_.by_key(pri_dct, ks, fill=False)
        for ks in map(atom_keys, gras)
    ]
    pri_dct = {}
    for gra_, pri_dct_ in zip(gras, pri_dcts):
        pri_dct_ = _break_priority_ties(gra_, pri_dct=pri_dct_, neg_hkeys=neg_hkeys)
        pri_dct.update(pri_dct_)

    return pri_dct


def _refine_priorities(
    gra, pri_dct: Optional[Dict[int, int]] = None, neg_hkeys: bool = True
):
    """Refine the canonical priorities for this graph based on some sort value

    (Only for connected graphs)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param neg_hkeys: Negate hydrogen keys? defaults to True
    :type neg_hkeys: bool, optional
    :param srt_eval_: An evaluator for sort values, based on current
        priorities. Curried such that srt_val_(pri_dct)(key) returns the
        sort value.
    """
    keys = sorted(atom_keys(gra))
    pri_dct = {k: 0 for k in keys} if pri_dct is None else pri_dct

    # For the refinement algorithm, remove minus signs from hydrogen priorities
    pri_dct = dict_.transform_values(pri_dct, abs)

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

    # Now that the algorithm is complete, restore minus signs to hydrogen priorities
    if neg_hkeys:
        hkeys = atom_keys(gra, symb="H")
        pri_dct = {k: -abs(p) if k in hkeys else p for k, p in pri_dct.items()}

    return pri_dct


def _break_priority_ties(gra, pri_dct: Dict[int, int], neg_hkeys: bool = False):
    """Break ties within priority classes.

    (Only for connected graphs)

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param neg_hkeys: Negate hydrogen keys? defaults to False
    :type neg_hkeys: bool, optional
    :param srt_eval_: An evaluator for sort values, based on current class
        indices. Curried such that srt_val_(pri_dct)(key) returns the sort
        value.
    """
    pri_dct = dict_.transform_values(pri_dct, abs)

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
        pri_dct = _refine_priorities(gra, pri_dct, neg_hkeys=False)

        # Update the list of classes needing further tie breaking
        cla_dct = class_dict_from_priority_dict(pri_dct)
        new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
        new_clas = sorted(new_cla_dct.items(), reverse=True)

    # Now that the algorithm is complete, restore minus signs to hydrogen priorities
    if neg_hkeys:
        hkeys = atom_keys(gra, symb="H")
        pri_dct = {k: -abs(p) if k in hkeys else p for k, p in pri_dct.items()}

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

    def _bond_type(bnd_keys):
        bords = dict_.values_by_key(bord_dct, bnd_keys)
        bords = tuple(sorted(map(_replace_reacting_bond_order, bords)))
        return bords

    btyps_dct = dict_.transform_values(atoms_bond_keys(imp_gra, ts_=True), _bond_type)

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
    btyps_dct.update({k: [-1] for k in hyd_keys})
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
            btyps = btyps_dct[key]
            nidxs = tuple(sorted(map(pri_dct.get, nkeys_dct[key])))
            return (symb, deg, hnum, mnum, apar, bpars, btyps, nidxs)

        return _value

    return _evaluator


# # core algorithm helpers
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
