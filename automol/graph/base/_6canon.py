""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import functools
import itertools
import numbers
from collections import abc
from typing import Any, Callable, Dict, List, Optional

import numpy
from phydat import ptab

from automol import util
from automol.graph.base._0core import (
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
    local_stereo_priorities,
    mass_numbers,
    nonbackbone_hydrogen_keys,
    set_atom_stereo_parities,
    set_stereo_parities,
    stereo_parities,
    tetrahedral_atom_keys,
    ts_breaking_bond_keys,
    ts_forming_bond_keys,
    ts_reverse,
    without_bonds_by_orders,
    without_dummy_atoms,
    without_pi_bonds,
    without_stereo,
)
from automol.graph.base._0core import (
    relabel as relabel_,
)
from automol.graph.base._2algo import (
    connected_components,
    is_connected,
    rings_bond_keys,
)
from automol.graph.base._3kekule import rigid_planar_bond_keys
from automol.graph.base._5geom import geometry_atom_parity, geometry_bond_parity
from automol.util import dict_

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
    ce_gra, _, is_refl = canonical_enantiomer_with_keys(gra, relabel=relabel)
    return ce_gra, is_refl


def canonical_enantiomer_with_keys(gra, relabel: bool = False):
    """Determine the canonical graph of the canonical enantiomer, along with
    the canonical key mapping.

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
    if not has_atom_stereo(gra):
        egra = gra
        ecan_key_dct = canonical_keys(gra, backbone_only=False)
        is_refl = None
    else:
        assert is_connected(gra), f"Requires a connected graph.\n{gra}"
        # Calculate canonical keys for the unreflected graph while converting
        # to the local stereo representation
        ugra = gra
        upri_dct, uloc_gra, _ = canonical_priorities_and_stereo_parities(
            ugra,
            backbone_only=False,
            par_eval_=parity_evaluator_read_canonical_(),
            par_eval2_=parity_evaluator_flip_local_(),
        )

        # Reflect the graph in the local stereo representation
        rloc_gra = reflect_local_stereo(uloc_gra)

        # Determine canonical keys for the reflected graph while converting
        # back to the canonical stereo representation
        rpri_dct, rgra, _ = canonical_priorities_and_stereo_parities(
            rloc_gra,
            backbone_only=False,
            break_ties=True,
            par_eval_=parity_evaluator_flip_local_(),
            par_eval2_=parity_evaluator_flip_local_(),
        )

        is_can = is_canonical_enantiomer(ugra, upri_dct, rgra, rpri_dct)

        if is_can in (True, None):
            egra = ugra
            epri_dct = upri_dct
        else:
            egra = rgra
            epri_dct = rpri_dct

        ecan_key_dct = break_priority_ties(egra, epri_dct)
        is_refl = None if is_can is None else (not is_can)

    if relabel:
        egra = relabel_(egra, ecan_key_dct)

    return egra, ecan_key_dct, is_refl


def canonical_ts_direction(gra):
    """If this is a TS graph, determine the canonical direction for it

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: a canonicalized graph of the canonical TS direction, along with
        a boolean flag indicating whether or not the TS graph has been
        reversed; `True` indicates it has been, `False` indicates it
        hasn't, and `None` indicates that it isn't a TS graph
    :rtype: (automol graph data structure, bool)
    """
    if not is_ts_graph(gra):
        cd_gra = gra
        is_rev = None
    else:
        assert is_connected(gra), f"Requires a connected graph.\n{gra}"
        # Calculate canonical keys for the unreflected graph while converting
        # to the local stereo representation
        ftsg = gra
        rtsg = ts_reverse(gra)

        _, _, is_can = canonical_priorities_and_stereo_parities(gra)

        is_rev = not is_can
        cd_gra = rtsg if is_rev else ftsg

    return cd_gra, is_rev


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
    :param break_ties: Break ties after keys have been refined?
    :type break_ties: bool
    :returns: a dictionary of canonical keys by atom key
    :rtype: dict[int: int]
    """
    assert is_connected(gra), f"Requires a connected graph.\n{gra}"
    atm_par_dct0 = atom_stereo_parities(gra)
    bnd_par_dct0 = bond_stereo_parities(gra)

    can_key_dct, gra, _ = canonical_priorities_and_stereo_parities(
        gra, backbone_only=backbone_only, break_ties=True
    )

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


def stereo_assignment_representation(gra, pri_dct):
    """Generate a representation of a stereo assignment, for checking for
    symmetric equivalence or for determining a canonical enantiomer

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :returns: A canonical representation of the assignment
    """
    atm_keys = sorted(atom_stereo_keys(gra), key=pri_dct.__getitem__)
    bnd_keys = sorted(
        bond_stereo_keys(gra), key=lambda x: sorted(map(pri_dct.__getitem__, x))
    )

    atm_pris = tuple([pri_dct[k]] for k in atm_keys)
    bnd_pris = tuple(sorted(map(pri_dct.__getitem__, k)) for k in bnd_keys)

    atm_pars = dict_.values_by_key(atom_stereo_parities(gra), atm_keys)
    bnd_pars = dict_.values_by_key(bond_stereo_parities(gra), bnd_keys)

    pris = atm_pris + bnd_pris
    pars = atm_pars + bnd_pars
    rep = tuple(zip(pris, pars))
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


def ts_direction_representation(tsg, pri_dct):
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
        sorted(sorted(map(pri_dct.__getitem__, k)) for k in frm_keys),
        sorted(sorted(map(pri_dct.__getitem__, k)) for k in brk_keys),
    )
    # Rep value 3: Stereochemistry
    ste_rep = stereo_assignment_representation(tsg, pri_dct)
    rep = (rxn_rep1, rxn_rep2, ste_rep)
    return rep


def ts_is_canonical_direction(ftsg, fpri_dct, rtsg, rpri_dct):
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
    ftsg = implicit(ftsg)
    fpri_dct = dict_.by_key(fpri_dct, atom_keys(ftsg))
    frep = ts_direction_representation(ftsg, fpri_dct)

    rtsg = implicit(rtsg)
    rpri_dct = dict_.by_key(rpri_dct, atom_keys(rtsg))
    rrep = ts_direction_representation(rtsg, rpri_dct)
    return frep < rrep


# # canonical stereo functions
def stereogenic_atom_keys(gra, pri_dct=None, assigned=False):
    """Find stereogenic atoms in this graph.

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    atoms will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include atoms that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic atom keys
    :rtype: frozenset
    """
    gra = without_dummy_atoms(gra)
    pri_dct = (
        canonical_priorities(gra, backbone_only=False) if pri_dct is None else pri_dct
    )
    ste_atm_keys = stereogenic_atom_keys_from_priorities(
        gra, pri_dct=pri_dct, assigned=assigned
    )
    return ste_atm_keys


def stereogenic_bond_keys(gra, pri_dct=None, assigned=False):
    """Find stereogenic bonds in this graph.

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    bonds will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include bonds that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic bond keys
    :rtype: frozenset
    """
    gra = without_dummy_atoms(gra)
    pri_dct = (
        canonical_priorities(gra, backbone_only=False) if pri_dct is None else pri_dct
    )
    ste_bnd_keys = stereogenic_bond_keys_from_priorities(
        gra, pri_dct=pri_dct, assigned=assigned
    )
    return ste_bnd_keys


def stereogenic_keys(gra, pri_dct=None, assigned=False):
    """Find stereogenic atoms and bonds in this graph.

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    atoms will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include atoms/bonds that already have assignments?
    :param assigned: bool
    :returns: keys to stereogenic atoms and bonds
    :rtype: frozenset
    """
    gra = without_dummy_atoms(gra)
    pri_dct = (
        canonical_priorities(gra, backbone_only=False) if pri_dct is None else pri_dct
    )
    ste_keys = stereogenic_keys_from_priorities(gra, pri_dct, assigned=assigned)
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
        _, loc_gra, _ = canonical_priorities_and_stereo_parities(
            can_gra,
            backbone_only=False,
            break_ties=False,
            par_eval_=parity_evaluator_read_canonical_(),
            par_eval2_=parity_evaluator_flip_local_(),
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
        _, can_gra, _ = canonical_priorities_and_stereo_parities(
            loc_gra,
            backbone_only=False,
            break_ties=False,
            par_eval_=parity_evaluator_flip_local_(),
            par_eval2_=parity_evaluator_flip_local_(),
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
        geo, local_stereo=False, geo_idx_dct=geo_idx_dct
    )

    par_eval2_ = (
        parity_evaluator_from_geometry_(geo, local_stereo=True, geo_idx_dct=geo_idx_dct)
        if local_stereo
        else None
    )

    _, gra, _ = canonical_priorities_and_stereo_parities(
        gra, par_eval_=par_eval_, par_eval2_=par_eval2_
    )

    return gra


# # core algorithm functions
def canonical_priorities(gra, backbone_only=True, break_ties=False, pri_dct=None):
    """Determine canonical priorities for this graph's atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :param break_ties: Break ties after priorities have been refined?
    :type break_ties: bool
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: dict[int: int]
    :returns: A dictionary of canonical priorities by atom key.
    :rtype: dict[int: int]
    """
    pri_dct, _, _ = canonical_priorities_and_stereo_parities(
        gra, backbone_only=backbone_only, break_ties=break_ties, pri_dct=pri_dct
    )
    return pri_dct


def canonical_priorities_and_stereo_parities(
    gra,
    par_eval_: Optional[ParityEvaluator] = None,
    par_eval2_: Optional[ParityEvaluator] = None,
    break_ties=False,
    backbone_only=True,
    pri_dct=None,
):
    """Determine canonical priorities and assign stereo parities to this graph

    This is how the parity evaluators are to be called:
    >>> par = par_eval_(gra, pri_dct)(key)  # this returns the parity

    :param gra: a molecular graph
    :type gra: automol graph data structure
    :param par_eval_: A parity evaluator for assigning parities during the
        priority calculation.
    :param par_eval2_: An optional second parity evaluator for assigning
        the parities that will be returned, if different from those to be
        used for the priority calculation.
    :param break_ties: Break ties to determine canonical keys from
        canonical priorities?
    :type break_ties: bool
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: dict[int: int]
    :returns: A dictionary of canonical priorities by atom key, a graph
        with stereo assignments, and a flag indicating for TS graphs
        whether or not they have a canonical direction.
    :rtype: dict[int: int], molecular graph data structure
    """
    orig_gra = gra
    gra = without_dummy_atoms(gra)

    par_eval_ = parity_evaluator_read_canonical_() if par_eval_ is None else par_eval_

    # 1. Iteratively assign parities and refine priorities.
    def _algo(gra_, pri_dct_, ts_rev=False):
        # Perform a TS graph reversal, if requested
        # (Note that we *keep* stereo in `gra0`, if present, because it may
        # need to be read out by the parity evaluators.)
        gra0 = ts_reverse(gra_) if ts_rev else gra_

        # Graph 1 will be for the priority calculation, graph 2 for the parity
        # assignments that will be returned.
        gra1 = gra2 = without_stereo(gra0)

        pri_dct0 = 0
        pri_dct1 = pri_dct_
        while pri_dct0 != pri_dct1:
            # a. Store the current priorities for comparison.
            pri_dct0 = pri_dct1

            # b. Refine priorities based on the assignments in graph 1.
            pri_dct1 = refine_priorities(gra1, pri_dct1)

            # c. Find stereogenic atoms and bonds based on current priorities.
            keys = stereogenic_keys_from_priorities(gra1, pri_dct1)

            # d. If there are none, the calculation is complete. Exit the loop.
            if not keys:
                break

            # e. Assign parities to graph 1 using the first parity evaluator.
            par_dct1 = par_eval_(gra0, pri_dct1, keys, ts_rev=ts_rev)
            gra1 = set_stereo_parities(gra1, par_dct1)

            # f. Assign parities to graph 2 using the second parity evaluator.
            if par_eval2_ is None:
                gra2 = gra1
            else:
                par_dct2 = par_eval2_(gra0, pri_dct1, keys, ts_rev=ts_rev)
                gra2 = set_stereo_parities(gra2, par_dct2)

        # If the graph was reversed, restore `gra2` to the original form
        gra2 = ts_reverse(gra2) if ts_rev else gra2
        return pri_dct1, gra1, gra2

    pri_dct, gra1, gra2 = _algo(gra, pri_dct)

    # 2. If this is a TS graph, rerun the algorithm on the reverse direction to
    # figure out which one is canonical
    is_can = None
    if is_ts_graph(gra):
        rpri_dct, rgra1, rgra2 = _algo(gra, pri_dct, ts_rev=True)
        if ts_is_canonical_direction(gra1, pri_dct, rgra1, rpri_dct):
            is_can = True
        else:
            pri_dct = rpri_dct
            gra1 = rgra1
            gra2 = rgra2
            is_can = False

    # 3. If requested, break priority ties to determine canonical keys.
    if break_ties:
        pri_dct = break_priority_ties(gra1, pri_dct)

    # 4. If requested, add in priorities for explicit hydrogens.
    if not backbone_only:
        pri_dct = reassign_hydrogen_priorities(gra, pri_dct, break_ties=break_ties)

    gra = set_stereo_parities(orig_gra, stereo_parities(gra2))
    return pri_dct, gra, is_can


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


def reassign_hydrogen_priorities(gra, pri_dct, break_ties=False, neg=False):
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
        cpri_dct = _reassign_hydrogen_priorities(
            cgra, pri_dct=cpri_dct, break_ties=break_ties, neg=neg
        )
        pri_dct.update(cpri_dct)

    return pri_dct


def _refine_priorities(gra, pri_dct=None, _backbone_only=True):
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
                ngb_idxs |= frozenset(map(pri_dct.__getitem__, ngb_keys))

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


def _reassign_hydrogen_priorities(gra, pri_dct, break_ties=False, neg=False):
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

    # If requested, break any ties
    if break_ties:
        pri_dct = _break_priority_ties(gra, pri_dct)

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

    def _evaluator(pri_dct):
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
            nidxs = tuple(sorted(map(pri_dct.__getitem__, nkeys_dct[key])))
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
        gra, pri_dct: Dict[int, int], keys: List[int], ts_rev: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param ts_rev: Is this a reversed TS graph?
        :type ts_rev: bool
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
        gra = without_dummy_atoms(gra)

        pri_dct = reassign_hydrogen_priorities(gra, pri_dct, break_ties=False, neg=True)

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
            par_dct = loc_par_eval_(can_gra, pri_dct, keys, ts_rev=ts_rev)

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
        gra, pri_dct: Dict[int, int], keys: List[int], ts_rev: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param ts_rev: Is this a reversed TS graph?
        :type ts_rev: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """
        # Do-nothing lines to prevent linting complaint
        assert pri_dct or not pri_dct
        assert ts_rev or not ts_rev

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
        gra, pri_dct: Dict[int, int], keys: List[int], ts_rev: bool = False
    ) -> Dict[int, int]:
        """Parity evaluator based on current priorities.

        :param gra: A molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: Dict[int, int]
        :param keys: The keys to evaluate parities for
        :type keys: List[int]
        :param ts_rev: Is this a reversed TS graph?
        :type ts_rev: bool
        :returns: A dictionary of parities, by key
        :rtype: Dict[int, int]
        """
        gra = explicit(gra)
        par_dct = stereo_parities(gra)

        loc_pri_dct = local_stereo_priorities(gra)

        pri_dct = reassign_hydrogen_priorities(gra, pri_dct, break_ties=False, neg=True)

        orig_gra = ts_reverse(gra) if ts_rev else gra

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                par = par_dct[key]

                nkeys = atom_stereo_sorted_neighbor_keys(gra, key)
                can_srt = sorted(nkeys, key=pri_dct.__getitem__)
                loc_srt = sorted(nkeys, key=loc_pri_dct.__getitem__)

                is_odd = util.is_odd_permutation(loc_srt, can_srt)
                ret_par = is_odd ^ par

                # If the original nkeys were different, assume we have an
                # Sn2-like situation. In this case, we need to account for the
                # difference in local parity between forward and reverse
                # directions.
                orig_nkeys = atom_stereo_sorted_neighbor_keys(orig_gra, key)
                if sorted(nkeys) != sorted(orig_nkeys):
                    # Identify forming and breaking keys involving this atom
                    fbnd_keys = [k for k in ts_forming_bond_keys(orig_gra) if key in k]
                    bbnd_keys = [k for k in ts_breaking_bond_keys(orig_gra) if key in k]

                    # Check that this is really an Sn2 reaction, in which case
                    # there will be exactly one of each
                    assert len(fbnd_keys) == len(bbnd_keys) == 1, (
                        f"Unforeseen case! Reversal alters stereo atom"
                        f"neighbors but TS is not Sn2-like. Investigate!"
                        f"\n{nkeys} != {orig_nkeys}\n{gra}\norig:\n{orig_gra}"
                    )

                    # Identify the forming and breaking neighbor keys
                    (fatm_key,) = fbnd_keys[0] - {key}
                    (batm_key,) = bbnd_keys[0] - {key}
                    assert fatm_key != batm_key, (
                        f"Breaking and forming bonds shouldn't match:"
                        f"{fbnd_keys[0]} == {bbnd_keys[0]}\n{gra}"
                    )

                    # Get locally sorted neighbors based on the original graph
                    orig_loc_srt = sorted(orig_nkeys, key=loc_pri_dct.__getitem__)

                    # The local sort vector of the original graph contains the
                    # breaking atom, while that of the current graph contains
                    # the forming atom.
                    assert batm_key in orig_loc_srt
                    assert fatm_key in loc_srt

                    # Replace the breaking atom with the forming atom for sort
                    # comparison
                    orig_loc_srt[orig_loc_srt.index(batm_key)] = fatm_key

                    # If local ordering for the original graph is an even
                    # permuation of the current one, we flip the parity.
                    is_even = util.is_even_permutation(loc_srt, orig_loc_srt)
                    ret_par ^= is_even

            # Otherwise, this is a bond
            else:
                assert (
                    isinstance(key, abc.Collection) and len(key) == 2
                ), f"{key} is not a valid bond key."
                par = par_dct[key]

                key1, key2 = key
                nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(
                    gra, key1, key2, pri_dct=pri_dct
                )

                can_nmax1 = nkey1s[-1]
                can_nmax2 = nkey2s[-1]
                loc_nmax1 = max(nkey1s, key=loc_pri_dct.__getitem__)
                loc_nmax2 = max(nkey2s, key=loc_pri_dct.__getitem__)

                if not (loc_nmax1 == can_nmax1) ^ (loc_nmax2 == can_nmax2):
                    ret_par = par
                else:
                    ret_par = not par

                # For bonds, the neighbors used shouldn't change upon TS
                # reversal for any of the cases we have dealt with. If this
                # isn't the case, we need to extend this code to handle the
                # unforeseen cases.
                orig_nkey1s, orig_nkey2s = bond_stereo_sorted_neighbor_keys(
                    orig_gra, key1, key2
                )
                assert sorted(nkey1s) == sorted(orig_nkey1s), (
                    f"Unforeseen case! Reversal alters stereo bond neighbors:"
                    f"\n{nkey1s} != {orig_nkey1s}\n{gra}\norig:\n{orig_gra}"
                )
                assert sorted(nkey2s) == sorted(orig_nkey2s), (
                    f"Unforeseen case! Reversal alters stereo bond neighbors:"
                    f"\n{nkey2s} != {orig_nkey2s}\n{gra}\norig:\n{orig_gra}"
                )

            return ret_par

        return {k: _parity(k) for k in keys}

    return _evaluator


# # core algorithm helpers
def stereogenic_keys_from_priorities(gra, pri_dct, assigned=False):
    """Find stereogenic atoms and bonds in this graph, given a set of atom
    priority values

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    bonds will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include bonds that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic atom and bond keys
    :rtype: frozenset
    """
    ste_atm_keys = stereogenic_atom_keys_from_priorities(
        gra, pri_dct, assigned=assigned
    )
    ste_bnd_keys = stereogenic_bond_keys_from_priorities(
        gra, pri_dct, assigned=assigned
    )
    return ste_atm_keys | ste_bnd_keys


def stereogenic_atom_keys_from_priorities(gra, pri_dct, assigned=False):
    """Find stereogenic atoms in this graph, given a set of atom priority
    values

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    atoms will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include atoms that already have stereo assignments?
    :type assigned: bool
    :returns: the stereogenic atom keys
    :rtype: frozenset[int]
    """
    gra = without_pi_bonds(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = reassign_hydrogen_priorities(gra, pri_dct, break_ties=False)

    atm_keys = tetrahedral_atom_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        atm_keys -= atom_stereo_keys(gra)

    def _is_stereogenic(key):
        nkeys = atom_stereo_sorted_neighbor_keys(gra, key, pri_dct=pri_dct)
        assert len(nkeys) <= 4, f"Too many neighbors! {nkeys}"
        pris = list(map(pri_dct.__getitem__, nkeys))
        return len(set(pris)) == len(pris)

    ste_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_atm_keys


def stereogenic_bond_keys_from_priorities(gra, pri_dct, assigned=False):
    """Find stereogenic bonds in this graph, given a set of atom priority
    values

    If the `assigned` flag is set to `False`, only  unassigned stereogenic
    bonds will be detected.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :param assigned: Include bonds that already have stereo assignments?
    :param assigned: bool
    :returns: the stereogenic bond keys
    :rtype: frozenset[frozenset[int]]
    """
    gra = without_pi_bonds(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = reassign_hydrogen_priorities(gra, pri_dct, break_ties=False)

    bnd_keys = rigid_planar_bond_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        bnd_keys -= bond_stereo_keys(gra)

    # Don't treat as a TS graph when checking for small rings
    rng_bnd_keys_lst = rings_bond_keys(gra, ts_=False)
    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union, filter(lambda x: len(x) < 8, rng_bnd_keys_lst), frozenset()
    )

    def _is_stereogenic(key):
        key1, key2 = key
        nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(
            gra, key1, key2, pri_dct=pri_dct
        )

        ret = True
        for nkeys in (nkey1s, nkey2s):
            assert len(nkeys) <= 2, f"Too many neighbors! {nkeys}"
            pris = list(map(pri_dct.__getitem__, nkeys))
            ret &= (
                False
                if not nkeys  # C=:O:
                else True
                if len(nkeys) == 1  # C=N:-X
                else len(set(pris)) == len(pris)  # C=C(-X)-Y
            )

        return ret

    ste_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_bnd_keys


def class_dict_from_priority_dict(pri_dct):
    """Obtain a class dictionary from a priority dictionary.

    :param pri_dct: A dictionary mapping atom keys to priorities.
    :type pri_dct: dict
    :returns: A dictionary mapping priorities onto the full set of keys
        for that priority class, as a sorted tuple.
    :rtype: dict[int: tuple]
    """
    keys = sorted(pri_dct.keys())
    clas = sorted(keys, key=pri_dct.__getitem__)
    cla_dct = {i: tuple(c) for i, c in itertools.groupby(clas, key=pri_dct.__getitem__)}
    return cla_dct
