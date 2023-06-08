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
from automol.graph.base._core import atom_keys
from automol.graph.base._core import backbone_keys
from automol.graph.base._core import bond_orders
from automol.graph.base._core import atom_stereo_keys
from automol.graph.base._core import bond_stereo_keys
from automol.graph.base._core import stereo_parities
from automol.graph.base._core import atom_stereo_parities
from automol.graph.base._core import bond_stereo_parities
from automol.graph.base._core import set_atom_stereo_parities
from automol.graph.base._core import set_bond_stereo_parities
from automol.graph.base._core import set_stereo_parities
from automol.graph.base._core import has_stereo
from automol.graph.base._core import has_atom_stereo
from automol.graph.base._core import atom_symbols
from automol.graph.base._core import mass_numbers
from automol.graph.base._core import stereo_candidate_atom_keys
from automol.graph.base._core import atoms_neighbor_atom_keys
from automol.graph.base._core import atoms_bond_keys
from automol.graph.base._core import implicit
from automol.graph.base._core import explicit
from automol.graph.base._core import atom_implicit_hydrogens
from automol.graph.base._core import backbone_hydrogen_keys
from automol.graph.base._core import nonbackbone_hydrogen_keys
from automol.graph.base._core import atom_nonbackbone_hydrogen_keys
from automol.graph.base._core import relabel
from automol.graph.base._core import without_pi_bonds
from automol.graph.base._core import without_stereo
from automol.graph.base._core import without_dummy_atoms
from automol.graph.base._core import ts_reactants_graph
from automol.graph.base._core import string as graph_string
from automol.graph.base._core import is_ts_graph
from automol.graph.base._algo import is_connected
from automol.graph.base._algo import connected_components
from automol.graph.base._algo import rings_bond_keys
from automol.graph.base._kekule import rigid_planar_bond_keys
from automol.graph.base._geom import geometry_atom_parity
from automol.graph.base._geom import geometry_bond_parity


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
    ce_gra, _, is_refl = canonical_enantiomer_with_keys(gra)
    return ce_gra, is_refl


def canonical_enantiomer_with_keys(gra):
    """ Determine the canonical graph of the canonical enantiomer, along with
    the canonical key mapping.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :returns: a canonicalized graph of the canonical enantiomer, along with
        a boolean flag indicating whether or not the graph has been
        reflected; `True` indicates it has been, `False` indicates it
        hasn't, and `None` indicates that it isn't an enantiomer
    :rtype: (automol graph data structure, bool)
    """
    if not has_atom_stereo(gra):
        ce_gra = gra
        is_refl = None
        ce_can_key_dct = canonical_keys(gra, backbone_only=False)
    else:
        # Calculate canonical keys for the unreflected graph while converting
        # to the local stereo representation
        ugra = gra
        ucan_key_dct, uloc_gra = calculate_priorities_and_assign_parities(
                ugra, backbone_only=False, break_ties=True,
                par_eval_=parity_evaluator_read_canonical_(ugra),
                par_eval2_=parity_evaluator_flip_local_(ugra))

        # Reflect the graph in the local stereo representation
        rloc_gra = reflect_local_stereo(uloc_gra)

        # Determine canonical keys for the reflected graph while converting
        # back to the canonical stereo representation
        rcan_key_dct, rgra = calculate_priorities_and_assign_parities(
                rloc_gra, backbone_only=False, break_ties=True,
                par_eval_=parity_evaluator_flip_local_(rloc_gra),
                par_eval2_=parity_evaluator_flip_local_(rloc_gra))

        urep = canonical_assignment_representation(ugra, ucan_key_dct)
        rrep = canonical_assignment_representation(rgra, rcan_key_dct)

        # If the reflected parities have a lower sort order, the molecule is
        # chiral and this is the mirror image of the canonical enantiomer.
        if urep > rrep:
            ce_gra = rgra
            ce_can_key_dct = rcan_key_dct
            is_refl = True
        else:
            ce_gra = ugra
            ce_can_key_dct = ucan_key_dct
            # If the sorted parities are the same, the molecule is achiral.
            if urep == rrep:
                is_refl = None
            # Otherwise, the unreflected parities have a lower sort order, so
            # the molecule is chiral and this is the canonical enantiomer.
            else:
                is_refl = False

    return ce_gra, ce_can_key_dct, is_refl


def canonical_assignment_representation(gra, pri_dct):
    """ Generate a canonical representation of a stereo assignment, for
    checking for symmetric equivalence or for determining a canonical
    enantiomer

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param pri_dct: A dictionary mapping atom keys to priorities
    :type pri_dct: dict
    :param par_dct: A dictionary mapping atom and bond keys to parities
    :type par_dct: dict
    :returns: A canonical representation of the assignment
    """
    atm_keys = sorted(atom_stereo_keys(gra), key=pri_dct.__getitem__)
    bnd_keys = sorted(bond_stereo_keys(gra),
                      key=lambda x: sorted(map(pri_dct.__getitem__, x)))

    atm_pris = tuple([pri_dct[k]] for k in atm_keys)
    bnd_pris = tuple(sorted(map(pri_dct.__getitem__, k)) for k in bnd_keys)

    atm_pars = dict_.values_by_key(atom_stereo_parities(gra), atm_keys)
    bnd_pars = dict_.values_by_key(bond_stereo_parities(gra), bnd_keys)

    pris = atm_pris + bnd_pris
    pars = atm_pars + bnd_pars
    rep = tuple(zip(pris, pars))
    return rep


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
    gra = without_dummy_atoms(gra)
    pri_dct = (canonical_priorities(gra, backbone_only=False)
               if pri_dct is None else pri_dct)
    ste_atm_keys = stereogenic_atom_keys_from_priorities(
        gra, pri_dct=pri_dct, assigned=assigned)
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
    gra = without_dummy_atoms(gra)
    pri_dct = (canonical_priorities(gra, backbone_only=False)
               if pri_dct is None else pri_dct)
    ste_bnd_keys = stereogenic_bond_keys_from_priorities(
        gra, pri_dct=pri_dct, assigned=assigned)
    return ste_bnd_keys


def stereogenic_keys(gra, pri_dct=None, assigned=False):
    """ Find stereogenic atoms and bonds in this graph.

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
    pri_dct = (canonical_priorities(gra, backbone_only=False)
               if pri_dct is None else pri_dct)
    ste_atm_keys = stereogenic_atom_keys(gra, pri_dct=pri_dct,
                                         assigned=assigned)
    ste_bnd_keys = stereogenic_bond_keys(gra, pri_dct=pri_dct,
                                         assigned=assigned)
    ste_keys = ste_atm_keys | ste_bnd_keys
    return ste_keys


def reflect(gra):
    """ Reflect the graph, locally inverting all stereo centers.

    To replicate the effect of reflecting the geometry, we convert to local
    stereo before reflecting and then convert back.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    """
    if has_atom_stereo(gra):
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


def to_local_stereo(gra, pri_dct=None):
    """ Convert canonical stereo parities to local ones

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :param pri_dct: priorities, to avoid recalculating
    :type pri_dct: dict[int: int]
    :returns: molecular graph with local stereo parities
    :rtype: automol graph data structure
    """
    loc_gra = without_stereo(gra)
    comps = connected_components(ts_reactants_graph(gra))
    for comp in comps:
        if has_stereo(comp):
            pri_dct_ = (None if pri_dct is None else
                        dict_.by_key(pri_dct, backbone_keys(comp)))
            _, loc_comp = calculate_priorities_and_assign_parities(
                    comp, backbone_only=False, break_ties=False,
                    par_eval_=parity_evaluator_read_canonical_(comp),
                    par_eval2_=parity_evaluator_flip_local_(comp),
                    pri_dct=pri_dct_)
        else:
            loc_comp = comp

        loc_gra = set_stereo_parities(loc_gra, stereo_parities(loc_comp))

    return loc_gra


def from_local_stereo(gra, pri_dct=None):
    """ Convert local stereo parities to canonical ones

        :param gra: molecular graph with local stereo parities
        :type gra: automol graph data structure
        :param pri_dct: priorities, to avoid recalculating
        :type pri_dct: dict[int: int]
        :returns: molecular graph with canonical stereo parities
        :rtype: automol graph data structure
    """
    can_gra = without_stereo(gra)
    loc_comps = connected_components(ts_reactants_graph(gra))
    for loc_comp in loc_comps:
        if has_stereo(loc_comp):
            pri_dct_ = (None if pri_dct is None else
                        dict_.by_key(pri_dct, backbone_keys(loc_comp)))
            _, can_comp = calculate_priorities_and_assign_parities(
                    loc_comp, backbone_only=False, break_ties=False,
                    par_eval_=parity_evaluator_flip_local_(loc_comp),
                    par_eval2_=parity_evaluator_flip_local_(loc_comp),
                    pri_dct=pri_dct_)
        else:
            can_comp = loc_comp

        can_gra = set_stereo_parities(can_gra, stereo_parities(can_comp))

    return can_gra


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
    ret_gra = without_stereo(gra)
    gra = without_dummy_atoms(gra)

    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(atm_keys))})

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
                         pri_dct=None, ts_=False):
    """ Determine canonical priorities for this graph's atoms

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param backbone_only: Consider backbone atoms only?
    :type backbone_only: bool
    :param break_ties: Break ties after priorities have been refined?
    :type break_ties: bool
    :param pri_dct: Optional initial priorities, to be refined.
    :type pri_dct: dict[int: int]
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
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
            pri_dct=pri_dct_, ts_=ts_)
        pri_dct.update(pri_dct_)
    return pri_dct


def calculate_priorities_and_assign_parities(
        gra, par_eval_=None, par_eval2_=None, break_ties=False,
        backbone_only=True, pri_dct=None, ts_=False):
    """ Determine canonical priorities and assign stereo parities to this graph

    This is how the parity evaluators are to be called:
    >>> par = par_eval_(pri_dct)(key)               # this returns the parity
    or, for TS graphs, we can pass in ts_select = 'R', 'P', or 'T' to select
    which type of stereochemistry we want:
    >>> par = par_eval_(pri_dct, ts_select)(key)    # this returns the parity

    :param gra: a connected molecular graph
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
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
    :returns: A dictionary of canonical priorities by atom key and a graph
        with stereo assignments.
    :rtype: dict[int: int], molecular graph data structure
    """
    assert is_connected(gra), "Not for disconnected graphs."
    assert gra == without_dummy_atoms(gra), (
        ("Remove dummy atoms:\n"
         f"{graph_string(gra)}\n"
         f"{graph_string(without_dummy_atoms(gra))}\n")
    )

    gra0 = gra

    par_eval_ = (parity_evaluator_read_canonical_(gra)
                 if par_eval_ is None else par_eval_)
    par_eval2_ = par_eval_ if par_eval2_ is None else par_eval2_

    # Graph 1 will be for the priority calculation, graph 2 for the parity
    # assignments that will be returned.
    gra1 = without_stereo(gra)
    gra2 = without_stereo(gra)

    # Work with an implicit graph for the priority calculation
    gra1 = implicit(gra1)

    # 1. Initially, put all atoms in the same priority class with priority 0.
    if pri_dct is None:
        pri_dct = {k: 0 for k in atom_keys(gra1)}

    # 2. Refine the initial priorities based on atom invariants, without stereo
    pri_dct = refine_priorities(gra1, pri_dct)

    # 3. Iteratively assign parities and refine priorities.
    pri_dct0 = None
    while pri_dct0 != pri_dct:
        # a. Find stereogenic atoms and bonds based on current priorities.
        atm_keys = stereogenic_atom_keys_from_priorities(
            gra1, pri_dct, ts_=ts_)
        bnd_keys = stereogenic_bond_keys_from_priorities(
            gra1, pri_dct, ts_=ts_)

        # b. If there are none, the calculation is complete. Exit the loop.
        if not atm_keys and not bnd_keys:
            break

        # These functions should take 'R', 'P', 'T' arguments, and we
        # will need to call each separately for ts graphs
        # c. Assign parities to graph 1 using the first parity evaluator.
        p1_ = par_eval_(pri_dct)
        gra1 = set_atom_stereo_parities(gra1, {k: p1_(k) for k in atm_keys})
        gra1 = set_bond_stereo_parities(gra1, {k: p1_(k) for k in bnd_keys})

        # These functions should take 'R', 'P', 'T' arguments, and we
        # will need to call each separately for ts graphs
        # d. Assign parities to graph 2 using the second parity evaluator.
        p2_ = par_eval2_(pri_dct)
        gra2 = set_atom_stereo_parities(gra2, {k: p2_(k) for k in atm_keys})
        gra2 = set_bond_stereo_parities(gra2, {k: p2_(k) for k in bnd_keys})

        # e. Store the current priorities for comparison.
        pri_dct0 = pri_dct

        # f. Refine priorities based on the assignments in graph 1.
        pri_dct = refine_priorities(gra1, pri_dct)

    # 4. If requested, break priority ties to determine canonical keys.
    if break_ties:
        pri_dct = break_priority_ties(gra1, pri_dct)

    # If requested, add in priorities for explicit hydrogens.
    if not backbone_only:
        pri_dct = assign_hydrogen_priorities(
            gra0, pri_dct, break_ties=break_ties)

    # Return the priorities calculated for graph 1, and return graph 2 with its
    # stereo assignments.
    return pri_dct, gra2


def refine_priorities(gra, pri_dct=None, srt_eval_=None):
    """ Refine the canonical priorities for this graph based on some sort value

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        :param srt_eval_: An evaluator for sort values, based on current
            priorities. Curried such that srt_val_(pri_dct)(key) returns the
            sort value.
    """
    if pri_dct is None:
        pri_dct = {k: 0 for k in atom_keys(gra)}

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

    ts_ = is_ts_graph(gra)

    nkeys_dct = atoms_neighbor_atom_keys(gra)
    bnds_dct = atoms_bond_keys(gra)

    symb_dct = dict_.transform_values(
        atom_symbols(gra), _hill_normalize_symbol)
    hnum_dct = atom_implicit_hydrogens(gra)
    mnum_dct = mass_numbers(gra)
    apar_dct = dict_.transform_values(atom_stereo_parities(gra), _replace_none)
    if ts_:
        apar_prd_dct = dict_.transform_values(
            atom_stereo_parities(gra, ts_select='P'), _replace_none)
        apar_ts_dct = dict_.transform_values(
            atom_stereo_parities(gra, ts_select='T'), _replace_none)

    def _sortable_bond_parities(bnd_par_dct):
        def _rep(bnd_keys):
            bpars = dict_.values_by_key(bnd_par_dct, bnd_keys)
            bpars = tuple(sorted(map(_replace_none, bpars)))
            return bpars

        bpars_dct = dict_.transform_values(bnds_dct, _rep)
        return bpars_dct

    bpars_dct = _sortable_bond_parities(bond_stereo_parities(gra))
    if ts_:
        bpars_prd_dct = _sortable_bond_parities(
            bond_stereo_parities(gra, ts_select='P'))
        bpars_ts_dct = _sortable_bond_parities(
            bond_stereo_parities(gra, ts_select='T'))

    bnd_ord_dct = bond_orders(gra)

    def _sortable_bond_orders(bnd_keys):
        bords = dict_.values_by_key(bnd_ord_dct, bnd_keys)
        bords = tuple(sorted(map(_replace_none, bords)))
        return bords

    bords_dct = dict_.transform_values(bnds_dct, _sortable_bond_orders)

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
            bords = bords_dct[key]
            nidxs = tuple(sorted(map(pri_dct.__getitem__, nkeys_dct[key])))
            val = (symb, deg, hnum, mnum, apar, bpars, bords, nidxs)
            if ts_:
                apar_prd = apar_prd_dct[key]
                apar_ts = apar_ts_dct[key]
                bpars_prd = bpars_prd_dct[key]
                bpars_ts = bpars_ts_dct[key]
                val += (apar_prd, apar_ts, bpars_prd, bpars_ts)
            return val

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
def parity_evaluator_from_geometry_(gra, geo=None, geo_idx_dct=None):
    r""" Determines stereo parity from a geometry

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
    gra = without_dummy_atoms(ts_reactants_graph(gra))

    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(atm_keys))})

    orig_geo = geo

    nkeys_dct = atoms_neighbor_atom_keys(gra)

    def _evaluator(pri_dct, geo=None):
        """ Parity evaluator based on current priorities.

            :param pri_dct: A dictionary mapping atom keys to priorities.
            :type pri_dct: dict
            :param geo: optionally, update the geometry
            :type geo: automol molecular geometry data structure
        """
        pri_dct = assign_hydrogen_priorities(
            gra, pri_dct, break_ties=False, neg=True)

        geo = geo if geo is not None else orig_geo

        def _parity(key):
            # If the key is a number, this is an atom
            if isinstance(key, numbers.Number):
                # Get the neighboring keys
                nkeys = nkeys_dct[key]

                # Sort them by priority
                nkeys = sorted(nkeys, key=pri_dct.__getitem__)

                # Get the atom parity
                par = geometry_atom_parity(
                    gra, geo, key, nkeys, geo_idx_dct=geo_idx_dct)
            # Otherwise, this is a bond
            else:
                assert isinstance(key, abc.Collection) and len(key) == 2, (
                    f"{key} is not a valid bond key.")
                key1, key2 = key
                nkey1s = nkeys_dct[key1] - {key2}
                nkey2s = nkeys_dct[key2] - {key1}

                nkey1s = sorted(nkey1s, key=pri_dct.__getitem__)
                nkey2s = sorted(nkey2s, key=pri_dct.__getitem__)

                # Get the bond parity
                par = geometry_bond_parity(
                    gra, geo, (key1, key2), (nkey1s, nkey2s),
                    geo_idx_dct=geo_idx_dct)

            return par

        return _parity

    return _evaluator


def parity_evaluator_read_canonical_(gra):
    """ Determines stereo parity from a graph with canonical stereo

    This is the trivial case, where parities in the graph are already
    assumed to be canonical and so nothing needs to be done to calculate
    them.

    :param gra: molecular graph with canonical stereo parities
    :type gra: automol graph data structure
    :returns: A parity evaluator, curried such that par_eval_(pri_dct)(key)
        returns the parity for a given atom, given a set of priorities.
    """

    def _evaluator(pri_dct, ts_select=None):
        """ Parity evaluator based on current priorities

        Class indices are ignored, since the parities are assumed to be
        canonical.

        :param pri_dct: A dictionary mapping atom keys to priorities
        :type pri_dct: dict
        """
        # Do-nothing line to prevent linting complaint
        assert pri_dct or not pri_dct

        atm_par_dct = atom_stereo_parities(gra, ts_select=ts_select)
        bnd_par_dct = bond_stereo_parities(gra, ts_select=ts_select)

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


def parity_evaluator_flip_local_(gra):
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

    loc_pri_dct = local_priority_dict(gra)

    def _evaluator(pri_dct):
        """ Parity evaluator based on current priorities

            Class indices are ignored, since the parities are assumed to be
            local.

            :param pri_dct: A dictionary mapping atom keys to priorities.
            :type pri_dct: dict
        """
        pri_dct = assign_hydrogen_priorities(
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
def stereogenic_atom_keys_from_priorities(gra, pri_dct, assigned=False,
                                          ts_=False):
    """ Find stereogenic atoms in this graph, given a set of canonical
        priorities

        If the `assigned` flag is set to `False`, only  unassigned stereogenic
        atoms will be detected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pri_dct: priorities, to avoid recalculating
        :type pri_dct: dict[int: int]
        :param assigned: Include atoms that already have stereo assignments?
        :type assigned: bool
        :param ts_: If this is a TS graph, treat it as such
        :type ts_: bool
        :returns: the stereogenic atom keys
        :rtype: frozenset
    """
    if not ts_:
        gra = ts_reactants_graph(gra)

    gra = without_pi_bonds(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = assign_hydrogen_priorities(
        gra, pri_dct, break_ties=False)

    atm_keys = stereo_candidate_atom_keys(gra)
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


def stereogenic_bond_keys_from_priorities(gra, pri_dct, assigned=False,
                                          ts_=False):
    """ Find stereogenic bonds in this graph, given a set of canonical
        priorities

        If the `assigned` flag is set to `False`, only  unassigned stereogenic
        bonds will be detected.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param pri_dct: priorities, to avoid recalculating
        :type pri_dct: dict[int: int]
        :param assigned: Include bonds that already have stereo assignments?
        :param assigned: bool
        :param ts_: If this is a TS graph, treat it as such
        :type ts_: bool
        :returns: the stereogenic bond keys
        :rtype: frozenset
    """
    if not ts_:
        gra = ts_reactants_graph(gra)

    gra = without_pi_bonds(gra)
    gra = explicit(gra)  # for simplicity, add the explicit hydrogens back in
    pri_dct = assign_hydrogen_priorities(
        gra, pri_dct, break_ties=False)

    bnd_keys = rigid_planar_bond_keys(gra)
    if not assigned:
        # Remove assigned stereo keys
        bnd_keys -= bond_stereo_keys(gra)

    # Don't treat as a TS graph when checking for small rings
    rng_bnd_keys_lst = rings_bond_keys(gra, ts_=False)
    bnd_keys -= functools.reduce(  # remove double bonds in small rings
        frozenset.union,
        filter(lambda x: len(x) < 8, rng_bnd_keys_lst), frozenset())

    nkeys_dct = atoms_neighbor_atom_keys(gra, ts_=False)

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


def assign_hydrogen_priorities(gra, pri_dct, break_ties=False, neg=False):
    """ Add explicit hydrogen keys to the priority dictionary

    :param neg: Negate hydrogen keys, to give them lowest priority?
    :param neg: bool
    :param ts_: If this is a TS graph, treat it as such
    :type ts_: bool
    """
    pri_dct = pri_dct.copy()

    # Backbone hydrogens are always included, so only non-backbone hydrogens
    # need to be assigned
    hyd_keys_pool = nonbackbone_hydrogen_keys(gra)

    if hyd_keys_pool:
        bbn_keys = sorted(backbone_keys(gra), key=pri_dct.__getitem__)
        bbn_pri_dct = dict_.by_key(pri_dct, bbn_keys)
        gra = explicit(gra)

        hyd_keys_dct = atom_nonbackbone_hydrogen_keys(gra)

        next_idx = max(bbn_pri_dct.values()) + 1

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
                hyd_pri_dct.update({k: next_idx for k in hyd_part})
                next_idx += len(hyd_part)
        # Otherwise, give each hydrogen a unique label.
        else:
            srt_hyd_keys = itertools.chain(
                *map(hyd_keys_dct.__getitem__, bbn_keys))
            hyd_pri_dct.update({k: (i+next_idx) for i, k in
                                enumerate(srt_hyd_keys)})

        hyd_pri_dct = dict_.by_key(hyd_pri_dct, hyd_keys_pool)
        pri_dct.update(hyd_pri_dct)

    # If requested, negate priorities for *all* hydrogen keys, for consistency
    if neg:
        all_hyd_keys = atom_keys(gra, symb='H')
        pri_dct = {k: -abs(p) if k in all_hyd_keys else p
                   for k, p in pri_dct.items()}

    return pri_dct


def local_priority_dict(gra):
    """ Generate a local ``priority'' dictionary
    """
    loc_pri_dct = {}
    loc_pri_dct.update(
        {k: k for k in backbone_keys(gra, hyd=False)})
    loc_pri_dct.update(
        {k: -abs(k) for k in backbone_hydrogen_keys(gra)})
    loc_pri_dct.update(
        {k: -numpy.inf for k in nonbackbone_hydrogen_keys(gra)})
    return loc_pri_dct


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
