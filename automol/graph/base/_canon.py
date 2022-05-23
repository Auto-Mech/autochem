""" canonicalization functions

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!

Reference:
Schneider, Sayle, Landrum. J. Chem. Inf. Model. 2015, 55, 10, 2111–2120
"""
import itertools
import functools
import numbers
from collections import abc
import numpy
from phydat import ptab
from automol import util
from automol.util import dict_
import automol.geom.base
from automol.graph.base._core import atom_keys
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import has_stereo
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import tetrahedral_atom_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import atom_implicit_hydrogen_valences
from automol.graph.base._core import atom_explicit_hydrogen_keys
from automol.graph.base._core import explicit_hydrogen_keys
from automol.graph.base._core import relabel
from automol.graph.base._core import without_bond_orders
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import union_from_sequence
from automol.graph.base._core import string as graph_string
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import connected_components
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._resonance import sp2_bond_keys


# # canonical key functions
def canonical_enantiomer(gra):
    """ Determine the canonical graph of the canonical enantiomer.

        Graphs with stereo will be reflected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: a canonicalized graph of the canonical enantiomer, along with
            a boolean flag indicating whether or not the graph has been
            reflected; `True` indicates it has been, `False` indicates it
            hasn't, and `None` indicates that it isn't an enantiomer
        :rtype: (automol graph data structure, bool)
    """
    can_enant_gra, is_reflected, _ = canonical_enantiomer_with_keys(gra)
    return can_enant_gra, is_reflected


def canonical_enantiomer_with_keys(gra):
    """ Determine the canonical graph of the canonical enantiomer, along with
        the canonical key mapping.

        Graphs with stereo will be reflected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :returns: a canonicalized graph of the canonical enantiomer, along with
            a boolean flag indicating whether or not the graph has been
            reflected; `True` indicates it has been, `False` indicates it
            hasn't, and `None` indicates that it isn't an enantiomer
        :rtype: (automol graph data structure, bool)
    """
    ste_atm_keys = atom_stereo_keys(gra)
    if not ste_atm_keys:
        is_reflected = None
        can_enant_key_dct = canonical_keys(gra, backbone_only=False)
        can_enant_gra = relabel(gra, can_enant_key_dct)
    else:
        # Calculate canonical keys for the unreflected graph while converting
        # to the local stereo representation
        ugra = gra
        uloc_gra, ucan_key_dct = _to_local_stereo_with_canonical_priorities(
            ugra, break_ties=True)

        # Reflect the graph in the local stereo representation
        rloc_gra = reflect_local_stereo(uloc_gra)

        # Determine canonical keys for the reflected graph while converting
        # back to the canonical stereo representation
        rgra, rcan_key_dct = _from_local_stereo_with_canonical_priorities(
            rloc_gra, break_ties=True)

        # Convert both to canonical graphs
        ucan_gra = relabel(ugra, ucan_key_dct)
        rcan_gra = relabel(rgra, rcan_key_dct)

        # Read and compare their parities
        ste_atm_keys = sorted(atom_stereo_keys(ucan_gra))
        assert ste_atm_keys == sorted(atom_stereo_keys(rcan_gra)), (
            "Sanity check. This should always be true.")
        uatm_par_dct = atom_stereo_parities(ucan_gra)
        ratm_par_dct = atom_stereo_parities(rcan_gra)

        uatm_pars = dict_.values_by_key(uatm_par_dct, ste_atm_keys)
        ratm_pars = dict_.values_by_key(ratm_par_dct, ste_atm_keys)

        # If the parities are the same, this is not an enantiomer.
        # If the unreflected parities have lower sort order
        if uatm_pars == ratm_pars:
            can_enant_gra = ucan_gra
            is_reflected = None
            can_enant_key_dct = ucan_key_dct
        # If the unreflected parities have lower sort order, don't reflect
        elif uatm_pars < ratm_pars:
            can_enant_gra = ucan_gra
            is_reflected = False
            can_enant_key_dct = ucan_key_dct
        # If the reflected parities have lower sort order, reflect
        else:
            can_enant_gra = rcan_gra
            is_reflected = True
            can_enant_key_dct = rcan_key_dct

    return can_enant_gra, is_reflected, can_enant_key_dct


def canonical(gra):
    """ A graph relabeled with canonical keys

        Stereo parities in the graph are assumed to be canonical.

        Requires a connected graph

        :param gra: a connected molecular graph with canonical stereo parities
        :type gra: automol graph data structure
        :returns: a new molecular graph with canonical keys; if explicit
            hydrogens are included, they will be relabeled as well
    """
    can_key_dct = canonical_keys(gra, backbone_only=False)
    return relabel(gra, can_key_dct)


def canonical_keys(gra, backbone_only=True):
    """ Determine canonical keys for this graph.

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
    atm_par_dct0 = atom_stereo_parities(gra)
    bnd_par_dct0 = bond_stereo_parities(gra)

    can_key_dct, gra = calculate_priorities_and_assign_parities(
        gra, backbone_only=backbone_only, break_ties=True)

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


# # canonical stereo functions
def stereogenic_atom_keys(gra, pri_dct=None, assigned=False):
    """ Find stereogenic atoms in this graph.

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
    gras = connected_components(gra)
    pri_dcts = [None if pri_dct is None else dict_.by_key(pri_dct, ks)
                for ks in map(atom_keys, gras)]
    ste_atm_keys = frozenset(itertools.chain(*(
        _stereogenic_atom_keys(g, pri_dct=d, assigned=assigned)
        for g, d in zip(gras, pri_dcts))))
    return ste_atm_keys


def stereogenic_bond_keys(gra, pri_dct=None, assigned=False):
    """ Find stereogenic bonds in this graph.

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
    gras = connected_components(gra)
    pri_dcts = [None if pri_dct is None else dict_.by_key(pri_dct, ks)
                for ks in map(atom_keys, gras)]
    ste_bnd_keys = frozenset(itertools.chain(*(
        _stereogenic_bond_keys(g, pri_dct=d, assigned=assigned)
        for g, d in zip(gras, pri_dcts))))
    return ste_bnd_keys


def reflect(gra):
    """ Reflect the graph, locally inverting all stereo centers.

        To replicate the effect of reflecting the geometry, we convert to local
        stereo before reflecting and then convert back.

        :param gra: molecular graph with canonical stereo parities
        :type gra: automol graph data structure
    """
    loc_gra = to_local_stereo(gra)
    loc_gra = reflect_local_stereo(loc_gra)
    gra = from_local_stereo(loc_gra)
    return gra


def reflect_local_stereo(gra):
    """ Reflect a graph with local stereo parities.

        Assuming local stereo parities, the parities can simply be reversed.

        :param gra: molecular graph with canonical stereo parities
        :type gra: automol graph data structure
    """
    atm_par_dct = atom_stereo_parities(gra)
    atm_par_dct = dict_.transform_values(
        atm_par_dct, lambda x: x if x is None else not x)
    gra = set_atom_stereo_parities(gra, atm_par_dct)
    return gra


def to_local_stereo(gra):
    """ Convert canonical stereo parities to local ones

        :param gra: molecular graph with canonical stereo parities
        :type gra: automol graph data structure
        :returns: molecular graph with local stereo parities
        :rtype: automol graph data structure
    """
    def _to_local_stereo_for_connected_component(gra):
        if has_stereo(gra):
            loc_gra, _ = _to_local_stereo_with_canonical_priorities(gra)
        else:
            loc_gra = gra
        return loc_gra

    gras = connected_components(gra)
    loc_gras = map(_to_local_stereo_for_connected_component, gras)
    loc_gra = union_from_sequence(loc_gras, shift_keys=False)
    return loc_gra


def _to_local_stereo_with_canonical_priorities(gra, break_ties=False):
    atm_par_dct0 = dict_.filter_by_value(
            atom_stereo_parities(gra), lambda x: x is not None)
    bnd_par_dct0 = dict_.filter_by_value(
            bond_stereo_parities(gra), lambda x: x is not None)

    pri_dct, gra = calculate_priorities_and_assign_parities(
            gra, backbone_only=False, break_ties=break_ties,
            par_eval_=parity_evaluator_from_canonical_stereo_(gra),
            ret_par_eval_=parity_evaluator_to_or_from_local_stereo_(gra))

    atm_par_dct = dict_.filter_by_value(
            atom_stereo_parities(gra), lambda x: x is not None)
    bnd_par_dct = dict_.filter_by_value(
            bond_stereo_parities(gra), lambda x: x is not None)

    assert set(atm_par_dct.keys()) == set(atm_par_dct0.keys()), (
        f"Something went wrong. Atom stereo keys don't match up:\n"
        f"input stereo atoms: {set(atm_par_dct0.keys())}\n"
        f"return stereo atoms: {set(atm_par_dct.keys())}\n"
    )

    assert set(bnd_par_dct.keys()) == set(bnd_par_dct0.keys()), (
        f"Something went wrong. Bond stereo keys don't match up:\n"
        f"input stereo bonds: {set(bnd_par_dct0.keys())}\n"
        f"return stereo bonds: {set(bnd_par_dct.keys())}\n"
    )

    return gra, pri_dct


def from_local_stereo(gra):
    """ Convert local stereo parities to canonical ones

        :param gra: molecular graph with local stereo parities
        :type gra: automol graph data structure
        :returns: molecular graph with canonical stereo parities
        :rtype: automol graph data structure
    """
    def _from_local_stereo_for_connected_component(gra):
        if has_stereo(gra):
            can_gra, _ = _from_local_stereo_with_canonical_priorities(gra)
        else:
            can_gra = gra
        return can_gra

    loc_gras = connected_components(gra)
    can_gras = map(_from_local_stereo_for_connected_component, loc_gras)
    can_gra = union_from_sequence(can_gras, shift_keys=False)
    return can_gra


def _from_local_stereo_with_canonical_priorities(gra, break_ties=False):
    atm_par_dct0 = dict_.filter_by_value(
            atom_stereo_parities(gra), lambda x: x is not None)
    bnd_par_dct0 = dict_.filter_by_value(
            bond_stereo_parities(gra), lambda x: x is not None)

    pri_dct, gra = calculate_priorities_and_assign_parities(
            gra, backbone_only=False, break_ties=break_ties,
            par_eval_=parity_evaluator_to_or_from_local_stereo_(gra),
            ret_par_eval_=parity_evaluator_to_or_from_local_stereo_(gra))

    atm_par_dct = dict_.filter_by_value(
            atom_stereo_parities(gra), lambda x: x is not None)
    bnd_par_dct = dict_.filter_by_value(
            bond_stereo_parities(gra), lambda x: x is not None)

    assert set(atm_par_dct.keys()) == set(atm_par_dct0.keys()), (
        f"Something went wrong. Atom stereo keys don't match up:\n"
        f"input stereo atoms: {set(atm_par_dct0.keys())}\n"
        f"return stereo atoms: {set(atm_par_dct.keys())}\n"
    )

    assert set(bnd_par_dct.keys()) == set(bnd_par_dct0.keys()), (
        f"Something went wrong. Bond stereo keys don't match up:\n"
        f"input stereo bonds: {set(bnd_par_dct0.keys())}\n"
        f"return stereo bonds: {set(bnd_par_dct.keys())}\n"
    )

    return gra, pri_dct


def set_stereo_from_geometry(gra, geo, geo_idx_dct=None):
    """ Determine stereo parities from a geometry

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :returns: molecular graph with stereo parities set from geometry;
            parities already present will be wiped out
        :rtype: automol graph data structure
    """
    ret_gra = without_stereo_parities(gra)
    gra = without_dummy_atoms(gra)

    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(atm_keys)})

    for comp in connected_components(gra):
        par_eval_ = parity_evaluator_from_geometry_(
            comp, geo, geo_idx_dct=geo_idx_dct)
        _, comp = calculate_priorities_and_assign_parities(
            comp, backbone_only=False, par_eval_=par_eval_)
        ret_gra = set_atom_stereo_parities(ret_gra, atom_stereo_parities(comp))
        ret_gra = set_bond_stereo_parities(ret_gra, bond_stereo_parities(comp))

    return ret_gra


# # core algorithm functions
def canonical_priorities(gra, backbone_only=True, break_ties=False,
                         pri_dct=None):
    """ Determine canonical priorities for this graph's atoms

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
    gras = connected_components(gra)
    pri_dcts = [None if pri_dct is None else dict_.by_key(pri_dct, ks)
                for ks in map(atom_keys, gras)]
    pri_dct = {}
    for gra_, pri_dct_ in zip(gras, pri_dcts):
        pri_dct_, _ = calculate_priorities_and_assign_parities(
            gra_, backbone_only=backbone_only, break_ties=break_ties,
            pri_dct=pri_dct_)
        pri_dct.update(pri_dct_)
    return pri_dct


def calculate_priorities_and_assign_parities(
        gra, par_eval_=None, ret_par_eval_=None, break_ties=False,
        backbone_only=True, pri_dct=None):
    """ Determine canonical priorities and assign stereo parities to this graph.

        One function to rule them all. This is where canonical ordering and
        canonical stereo parities are ensured.

        Requires a connected graph.

        The parity evaluators determine the parity of an atom or bond based on
        the current set of priorities as `par_eval_(pri_dct)(key)`, where
        `pri_dct` is the dictionary of priorities and `key` is the atom or
        bond key to be evaluated.

        :param gra: a connected molecular graph
        :type gra: automol graph data structure
        :param par_eval_: A parity evaluator that will be used to determine the
            canonical priorities. If `ret_par_eval_` is `None`, then it
            will also be used for to assign parities.
        :param ret_par_eval_: An optional alternative evaluator for the parity
            assignment, if different from the one used to determine canonical
            priorities.
        :param break_ties: Break ties after keys have been refined?
        :type break_ties: bool
        :param backbone_only: Consider backbone atoms only?
        :type backbone_only: bool
        :param pri_dct: Optional initial priorities, to be refined.
        :type pri_dct: dict[int: int]
        :returns: A dictionary of canonical priorities by atom key and a graph
            with stereo assignments.
        :rtype: dict[int: int], dict[int: bool], dict[frozenset: bool]
    """
    assert is_connected(gra), "Not for disconnected graphs."
    assert gra == without_dummy_atoms(gra), (
        ("Remove dummy atoms:\n"
         f"{graph_string(gra)}\n"
         f"{graph_string(without_dummy_atoms(gra))}\n"
        )
    )

    # Work with an implicit graph to determine canonical priorities for
    # backbone atoms
    gra0 = gra
    gra = implicit(gra)
    # Remove stereo parities for consistent canonical priority determination.
    # Parities will be iteratively introduced as the priorities are refined.
    gra = without_stereo_parities(gra)

    # By default, assume `gra0` has canonical stereo.
    do_stereo = has_stereo(gra0) or par_eval_ is not None
    par_eval_ = (parity_evaluator_from_canonical_stereo_(gra0)
                 if par_eval_ is None else par_eval_)
    ret_par_eval_ = par_eval_ if ret_par_eval_ is None else ret_par_eval_

    # If initial priorities aren't given, initialize them all to 0
    if pri_dct is None:
        pri_dct = {k: 0 for k in atom_keys(gra)}

    # 1. Refine the initial priorities based on atom invariants, without stereo
    pri_dct = refine_priorities(gra, pri_dct)

    # 2. Iteratively refine priorities while introducing stereo parities
    ret_atm_par_dct = {}
    ret_bnd_par_dct = {}
    if do_stereo:
        last_pri_dct = None
        while pri_dct != last_pri_dct:
            # a. Find stereogenic keys based on the current priorities. If
            # there aren't any, exit the loop.
            atm_keys = _stereogenic_atom_keys(gra, pri_dct)
            bnd_keys = _stereogenic_bond_keys(gra, pri_dct)
            if not atm_keys and not bnd_keys:
                break

            # b. Update parities based on the current priorities
            par_ = par_eval_(pri_dct)
            atm_par_dct = {k: par_(k) for k in atm_keys}
            bnd_par_dct = {k: par_(k) for k in bnd_keys}
            gra = set_atom_stereo_parities(gra, atm_par_dct)
            gra = set_bond_stereo_parities(gra, bnd_par_dct)

            # c. Assign parities to the graph. If ret_par_eval_ is different
            # from par_eval_, the assigned parities will be different from
            # those used for the priority determination.
            ret_par_ = ret_par_eval_(pri_dct)
            ret_atm_par_dct.update({k: ret_par_(k) for k in atm_keys})
            ret_bnd_par_dct.update({k: ret_par_(k) for k in bnd_keys})

            # e. Further refine priorities based on the new assignments
            last_pri_dct = pri_dct
            pri_dct = refine_priorities(gra, pri_dct)

    # 3. If requested, break ties based on keys.
    if break_ties:
        pri_dct = break_priority_ties(gra, pri_dct)

    # 4. If requested, add in priorities for explicit hydrogens.
    if not backbone_only:
        pri_dct = augment_priority_dict_with_hydrogen_keys(
            gra0, pri_dct, break_ties=break_ties)

    ret_gra = set_atom_stereo_parities(gra0, ret_atm_par_dct)
    ret_gra = set_bond_stereo_parities(ret_gra, ret_bnd_par_dct)
    return pri_dct, ret_gra


def refine_priorities(gra, pri_dct, srt_eval_=None):
    """ Refine the canonical priorities for this graph based on some sort value.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        :param srt_eval_: An evaluator for sort values, based on current
            priorities. Curried such that srt_val_(pri_dct)(key) returns the
            sort value.
    """
    if srt_eval_ is None:
        srt_eval_ = sort_evaluator_atom_invariants_(gra)

    pri_dct = pri_dct.copy()

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
                ngb_keys = frozenset.union(
                    *map(ngb_keys_dct.__getitem__, new_cla))

                # Get priorities of these neighboring atoms.
                ngb_idxs |= frozenset(map(pri_dct.__getitem__, ngb_keys))

                # Don't include classes that were already up for re-evaluation.
                ngb_idxs -= frozenset(dict(new_clas))

            for ngb_idx in sorted(ngb_idxs):
                ngb_cla = cla_dct[ngb_idx]
                if len(ngb_cla) > 1:
                    new_clas.insert(0, (ngb_idx, cla_dct[ngb_idx]))

    return pri_dct


def sort_evaluator_atom_invariants_(gra):
    """ A sort function based on atom invariants with two levels of currying.

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

    def _hill_normalize_symbol(symb):
        """ normalize atomic symbols to make them sort according to the hill
            system
            (C first, H second, others in alphabetical order)
        """
        symb = ptab.to_symbol(symb)
        if symb == 'C':
            symb = ''
        if symb == 'H':
            symb = '1'
        return symb

    ngb_keys_dct = atoms_neighbor_atom_keys(gra)
    bnds_dct = atoms_bond_keys(gra)

    symb_dct = dict_.transform_values(
        atom_symbols(gra), _hill_normalize_symbol)
    hnum_dct = atom_implicit_hydrogen_valences(gra)
    mnum_dct = mass_numbers(gra)
    apar_dct = dict_.transform_values(atom_stereo_parities(gra), _replace_none)

    bnd_par_dct = bond_stereo_parities(gra)

    def _sortable_bond_stereo_values(bnd_keys):
        bpars = dict_.values_by_key(bnd_par_dct, bnd_keys)
        bpars = tuple(sorted(map(_replace_none, bpars)))
        return bpars

    bpars_dct = dict_.transform_values(bnds_dct, _sortable_bond_stereo_values)

    def _evaluator(pri_dct):
        """ Sort value evaluator based on current priorities.

            :param pri_dct: A dictionary mapping atom keys to priorities
            :type pri_dct: dict
        """

        def _value(key):
            symb = symb_dct[key]        # symbol
            deg = len(bnds_dct[key])    # number of bonds
            hnum = hnum_dct[key]        # number of hydrogens
            mnum = mnum_dct[key]
            apar = apar_dct[key]
            bpars = bpars_dct[key]
            ngb_idxs = tuple(
                sorted(map(pri_dct.__getitem__, ngb_keys_dct[key])))
            return (symb, deg, hnum, mnum, apar, bpars, ngb_idxs)

        return _value

    return _evaluator


def break_priority_ties(gra, pri_dct):
    """ Break ties within priority classes.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        :param srt_eval_: An evaluator for sort values, based on current class
            indices. Curried such that srt_val_(pri_dct)(key) returns the sort
            value.
    """
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
        cla_dct = class_dict_from_priority_dict(pri_dct)
        pri_dct = refine_priorities(gra, pri_dct)

        # Update the list of classes needing further tie breaking
        cla_dct = class_dict_from_priority_dict(pri_dct)
        new_cla_dct = dict_.filter_by_value(cla_dct, lambda v: len(v) > 1)
        new_clas = sorted(new_cla_dct.items(), reverse=True)

    return pri_dct


# # parity evaluators
def parity_evaluator_from_canonical_stereo_(gra):
    """ Determines stereo parity from a graph with canonical stereo

        This is the trivial case, where parities in the graph are already
        assumed to be canonical and so nothing needs to be done to calculate
        them.

        :param gra: molecular graph with canonical stereo parities
        :type gra: automol graph data structure
        :returns: A parity evaluator, curried such that par_eval_(pri_dct)(key)
            returns the parity for a given atom, given a set of priorities.
    """
    atm_par_dct = atom_stereo_parities(gra)
    bnd_par_dct = bond_stereo_parities(gra)

    def _evaluator(pri_dct):
        """ Parity evaluator based on current priorities

            Class indices are ignored, since the parities are assumed to be
            canonical.

            :param pri_dct: A dictionary mapping atom keys to priorities
            :type pri_dct: dict
        """
        # Do-nothing line to prevent linting complaint
        assert pri_dct or not pri_dct

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                par = atm_par_dct[key]
            # Otherwise, this is a bond
            else:
                assert isinstance(key, abc.Collection) and len(key) == 2, (
                    f"{key} is not a valid bond key.")
                par = bnd_par_dct[key]
            return par

        return _parity

    return _evaluator


def parity_evaluator_from_geometry_(gra, geo=None, geo_idx_dct=None):
    r""" Determines stereo parity from a geometry

        Atom parity is defined as follows:

        The four keys passed in are apices of a tetrahedron. Looking at 2, 3,
        and 4 from 1, they will either ascend in clockwise or counterclockwise
        order.

        If ascending in counterclockwise order, the parity is False ('-').
        If ascending in clockwise order, the parity is True ('+').

              2                   2
             /1\                 /1\
            3---4               4---3

            counterclockwise    clockwise
            False               True
            '-'                 '+'

        (Viewed looking down from 1)

        If only three keys are passed in, they will be treated as keys 2, 3,
        and 4 above and it will be assumed that there is a lone pair at 1.

        Bond parity is defined as follows:

        For each atom in the double bond, find the heavy-atom neighbor with the
        higher canonical number. Althrough hydrogen atoms have higher canonical
        numbers, they are always given lowest priority.

        If the neighbors are cis to each other, the parity is False ('-').
        If the neighbors are trans to each other, the parity is True ('+').

            max    max      max    min
              \   /           \   /
               A=B             A=B
              /   \           /   \
            min    min      min    max

            cis             trans
            False           True
            '-'             '+'

        If one side only has a single neighbor, then it is compared with the
        maximum neighbor on the other side.

            max    nei      max
              \   /           \
               A=B             A=B
              /               /   \
            min             min    nei

            cis             trans
            False           True
            '-'             '+'

        If both sides have only single neighbors, then they are compared to
        each other.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :returns: A parity evaluator, curried such that par_eval_(pri_dct)(key)
            returns the parity for a given atom, given a set of priorities.
    """
    assert gra == explicit(gra), (
        "Explicit graph should be used when getting parities from geometry.")

    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(atm_keys)})

    orig_geo = geo

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _evaluator(pri_dct, geo=None):
        """ Parity evaluator based on current priorities.

            :param pri_dct: A dictionary mapping atom keys to priorities.
            :type pri_dct: dict
            :param geo: optionally, update the geometry
            :type geo: automol molecular geometry data structure
        """
        pri_dct = augment_priority_dict_with_hydrogen_keys(
            gra, pri_dct, break_ties=False, neg=True)

        geo = geo if geo is not None else orig_geo

        xyzs = automol.geom.base.coordinates(geo)
        xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in atm_keys}

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                # Get the neighboring keys
                nkeys = nkeys_dct[key]

                # Sort them by priority
                nkeys = sorted(nkeys, key=pri_dct.__getitem__)

                # If there are only three groups, use the stereo atom itself as
                # the top apex of the tetrahedron.
                if len(nkeys) == 4:
                    keys = nkeys
                else:
                    assert len(nkeys) == 3
                    keys = [key] + list(nkeys)

                xyzs = list(map(list, map(xyz_dct.__getitem__, keys)))
                det_mat = numpy.ones((4, 4))
                det_mat[:, 1:] = xyzs
                det_val = numpy.linalg.det(det_mat)
                assert det_val != 0.  # for now, assume no four-atom planes
                par = det_val > 0.
            # Otherwise, this is a bond
            else:
                assert isinstance(key, abc.Collection) and len(key) == 2, (
                    f"{key} is not a valid bond key.")
                key1, key2 = key
                nkey1s = nkeys_dct[key1] - {key2}
                nkey2s = nkeys_dct[key2] - {key1}

                nmax1 = max(nkey1s, key=pri_dct.__getitem__)
                nmax2 = max(nkey2s, key=pri_dct.__getitem__)

                xyz1 = xyz_dct[key1]
                xyz2 = xyz_dct[key2]
                nxyz1 = xyz_dct[nmax1]
                nxyz2 = xyz_dct[nmax2]

                bnd1_vec = numpy.subtract(nxyz1, xyz1)
                bnd2_vec = numpy.subtract(nxyz2, xyz2)

                dot_val = numpy.vdot(bnd1_vec, bnd2_vec)
                assert dot_val != 0.    # for now, assume not collinear
                par = dot_val < 0.
            return par

        return _parity

    return _evaluator


def parity_evaluator_to_or_from_local_stereo_(gra):
    """ Determines canonical from local stereo parity or vice versa (same
    operation)

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

        :param gra: molecular graph with local stereo parities
        :type gra: automol graph data structure
        :returns: A parity evaluator, curried such that par_eval_(pri_dct)(key)
            returns the parity for a given atom, given a set of priorities.
    """
    gra = explicit(gra)
    atm_par_dct = atom_stereo_parities(gra)
    bnd_par_dct = bond_stereo_parities(gra)
    nkeys_dct = atoms_neighbor_atom_keys(gra)

    loc_pri_dct = {}
    loc_pri_dct.update({k: k for k in backbone_keys(gra)})
    loc_pri_dct.update({k: -numpy.inf for k in explicit_hydrogen_keys(gra)})

    def _evaluator(pri_dct):
        """ Parity evaluator based on current priorities

            Class indices are ignored, since the parities are assumed to be
            local.

            :param pri_dct: A dictionary mapping atom keys to priorities.
            :type pri_dct: dict
        """
        pri_dct = augment_priority_dict_with_hydrogen_keys(
            gra, pri_dct, break_ties=False, neg=True)

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                par = atm_par_dct[key]

                loc_srt = sorted(nkeys_dct[key], key=loc_pri_dct.__getitem__)
                can_srt = sorted(nkeys_dct[key], key=pri_dct.__getitem__)

                if util.is_even_permutation(loc_srt, can_srt):
                    ret_par = par
                else:
                    ret_par = not par
            # Otherwise, this is a bond
            else:
                assert isinstance(key, abc.Collection) and len(key) == 2, (
                    f"{key} is not a valid bond key.")
                par = bnd_par_dct[key]

                key1, key2 = key
                nkey1s = nkeys_dct[key1] - {key2}
                nkey2s = nkeys_dct[key2] - {key1}

                loc_nmax1 = max(nkey1s, key=loc_pri_dct.__getitem__)
                loc_nmax2 = max(nkey2s, key=loc_pri_dct.__getitem__)
                can_nmax1 = max(nkey1s, key=pri_dct.__getitem__)
                can_nmax2 = max(nkey2s, key=pri_dct.__getitem__)

                if not (loc_nmax1 == can_nmax1) ^ (loc_nmax2 == can_nmax2):
                    ret_par = par
                else:
                    ret_par = not par

            return ret_par

        return _parity

    return _evaluator


# # core algorithm helpers
def _stereogenic_atom_keys(gra, pri_dct=None, assigned=False):
    """ Find stereogenic atoms in this graph, given a set of canonical priorities.

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
    pri_dct = (canonical_priorities(gra, backbone_only=False)
               if pri_dct is None else pri_dct)
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = augment_priority_dict_with_hydrogen_keys(gra, pri_dct,
                                                       break_ties=False)

    atm_keys = tetrahedral_atom_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        atm_keys -= atom_stereo_keys(gra)

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(key):
        nkeys = list(nkeys_dct[key])
        idxs = list(map(pri_dct.__getitem__, nkeys))
        return len(set(idxs)) == len(idxs)

    ste_atm_keys = frozenset(filter(_is_stereogenic, atm_keys))
    return ste_atm_keys


def _stereogenic_bond_keys(gra, pri_dct=None, assigned=False):
    """ Find stereogenic bonds in this graph, given a set of canonical priorities.

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
    pri_dct = (canonical_priorities(gra, backbone_only=False)
               if pri_dct is None else pri_dct)
    gra = without_bond_orders(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = augment_priority_dict_with_hydrogen_keys(gra, pri_dct,
                                                       break_ties=False)

    bnd_keys = sp2_bond_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        bnd_keys -= bond_stereo_keys(gra)

    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, rings_bond_keys(gra)), frozenset())

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _is_stereogenic(key):

        def _is_asymmetric_on_bond(atm1_key, atm2_key):
            nkeys = list(nkeys_dct[atm1_key] - {atm2_key})

            if not nkeys:                # C=:O:
                # Atoms without neighbors are automatically symmetric
                ret = False
            elif len(nkeys) == 1:        # C=N:-X
                # Atoms without 1 neighbor are automatically asymmetric
                ret = True
            else:
                # For atoms with 2 neighbors, we need to determine whether or
                # not they are symmetric from the priorities.
                assert len(nkeys) == 2   # C=C(-X)-Y
                ret = pri_dct[nkeys[0]] != pri_dct[nkeys[1]]

            return ret

        atm1_key, atm2_key = key
        return (_is_asymmetric_on_bond(atm1_key, atm2_key) and
                _is_asymmetric_on_bond(atm2_key, atm1_key))

    ste_bnd_keys = frozenset(filter(_is_stereogenic, bnd_keys))
    return ste_bnd_keys


def augment_priority_dict_with_hydrogen_keys(gra, pri_dct, break_ties=False,
                                             neg=False):
    """ Add explicit hydrogen keys to the priority dictionary

        :param neg: Negate the keys, to give them minimum (rather than maximum)
            priority?
        :param neg: bool
    """
    pri_dct = pri_dct.copy()

    hyd_keys_pool = explicit_hydrogen_keys(gra)

    sgn = -1 if neg else +1

    if hyd_keys_pool:
        bbn_keys = sorted(backbone_keys(gra), key=pri_dct.__getitem__)
        gra = explicit(gra)

        hyd_keys_dct = atom_explicit_hydrogen_keys(gra)

        next_idx = len(bbn_keys)

        # If not breaking ties, assign equal priority to hydrogens bonded to
        # atoms from the same priority class.
        hyd_pri_dct = {}
        if not break_ties:
            # Partition hydrogens into classes based on their parent atoms.
            bbn_parts = sorted_classes_from_priority_dict(pri_dct)
            get_hyd_key_ = hyd_keys_dct.__getitem__
            hyd_parts = [
                list(itertools.chain(*map(get_hyd_key_, bbn_part)))
                for bbn_part in bbn_parts]
            for hyd_part in hyd_parts:
                hyd_pri_dct.update({k: sgn*next_idx for k in hyd_part})
                next_idx += len(hyd_part)
        # Otherwise, give each hydrogen a unique label.
        else:
            srt_hyd_keys = itertools.chain(
                *map(hyd_keys_dct.__getitem__, bbn_keys))
            hyd_pri_dct.update({k: sgn*(i+next_idx) for i, k in
                                enumerate(srt_hyd_keys)})

        hyd_pri_dct = dict_.by_key(hyd_pri_dct, hyd_keys_pool)
        pri_dct.update(hyd_pri_dct)

    return pri_dct


def class_dict_from_priority_dict(pri_dct):
    """ Obtain a class dictionary from a priority dictionary.

        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: dict
        :returns: A dictionary mapping priorities onto the full set of keys
            for that priority class, as a sorted tuple.
        :rtype: dict[int: tuple]
    """
    keys = sorted(pri_dct.keys())
    clas = sorted(keys, key=pri_dct.__getitem__)
    cla_dct = {i: tuple(c)
               for i, c in itertools.groupby(clas, key=pri_dct.__getitem__)}
    return cla_dct


def sorted_classes_from_priority_dict(pri_dct):
    """ Obtain classes from index dict, sorted by priority.

        :param pri_dct: A dictionary mapping atom keys to priorities.
        :type pri_dct: dict
        :returns: A tuple of tuples of keys for each priority class, sorted by
            priority value.
        :rtype: tuple[tuple[int]]
    """
    keys = sorted(pri_dct.keys())
    clas = sorted(keys, key=pri_dct.__getitem__)
    cla_dct = tuple(
        tuple(c) for _, c in itertools.groupby(clas, key=pri_dct.__getitem__))
    return cla_dct
