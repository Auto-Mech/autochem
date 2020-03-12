""" molecular graph reaction identifiers

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.

Function return values:
    Each function returns a list of "transformations" (see graph.trans)
    describing bonds broken and formed to transform reactants into products,
    along with sort indices for putting the reactants and products in a
    standard sort order for each reaction.
"""
import itertools
from automol import formula
import automol.graph.trans as trans
import automol.convert.graph as graph_convert
from automol.graph._graph import atom_count as _atom_count
from automol.graph._graph import heavy_atom_count as _heavy_atom_count
from automol.graph._graph import electron_count as _electron_count
from automol.graph._graph import atom_keys as _atom_keys
from automol.graph._graph import bond_keys as _bond_keys
from automol.graph._graph import explicit as _explicit
from automol.graph._graph import union as _union
from automol.graph._graph import union_from_sequence as _union_from_sequence
from automol.graph._graph import connected_components as _connected_components
from automol.graph._graph import full_isomorphism as _full_isomorphism
from automol.graph._graph import add_bonds as _add_bonds
from automol.graph._graph import remove_bonds as _remove_bonds
from automol.graph._graph import (add_atom_explicit_hydrogen_keys as
                                  _add_atom_explicit_hydrogen_keys)
from automol.graph._graph import (unsaturated_atom_keys as
                                  _unsaturated_atom_keys)


def hydrogen_migration(rct_xgrs, prd_xgrs):
    """ find a hydrogen migration transformation

    Hydrogen migrations are identified by adding a hydrogen to an unsaturated
    site of the reactant and adding a hydrogen to an unsaturated site of the
    product and seeing if they match up. If so, we have a hydrogen migration
    between these two sites.
    """
    assert _is_valid_reagent_graph_list(rct_xgrs)
    assert _is_valid_reagent_graph_list(prd_xgrs)

    tras = []
    rct_idxs = None
    prd_idxs = None

    if len(rct_xgrs) == 1 and len(prd_xgrs) == 1:
        xgr1, = rct_xgrs
        xgr2, = prd_xgrs

        h_atm_key1 = max(_atom_keys(xgr1)) + 1
        h_atm_key2 = max(_atom_keys(xgr2)) + 1

        atm_keys1 = _unsaturated_atom_keys(xgr1)
        atm_keys2 = _unsaturated_atom_keys(xgr2)
        for atm_key1, atm_key2 in itertools.product(atm_keys1, atm_keys2):
            xgr1_h = _add_atom_explicit_hydrogen_keys(
                xgr1, {atm_key1: [h_atm_key1]})
            xgr2_h = _add_atom_explicit_hydrogen_keys(
                xgr2, {atm_key2: [h_atm_key2]})

            inv_atm_key_dct = _full_isomorphism(xgr2_h, xgr1_h)
            if inv_atm_key_dct:
                tra = trans.from_data(
                    frm_bnd_keys=[{atm_key1,
                                   inv_atm_key_dct[h_atm_key2]}],
                    brk_bnd_keys=[{inv_atm_key_dct[atm_key2],
                                   inv_atm_key_dct[h_atm_key2]}])
                tras.append(tra)

                rct_idxs = (0,)
                prd_idxs = (0,)

    tras = tuple(tras)
    return tras, rct_idxs, prd_idxs


def hydrogen_abstraction(rct_xgrs, prd_xgrs):
    """ find a hydrogen abstraction transformation

    Hydrogen abstractions are identified first by checking whether the
    molecular formulas are consistent with a reaction of the form R1H + R2 =>
    R2H + R1. If they do, we identify the abstraction sites by adding hydrogens
    to unsaturated sites of the R1 product to see if we get the R1H reactant.
    We then do the same for the R2 reactant and the R2H product.
    """
    assert _is_valid_reagent_graph_list(rct_xgrs)
    assert _is_valid_reagent_graph_list(prd_xgrs)

    tras = []
    rct_idxs = None
    prd_idxs = None

    if len(rct_xgrs) == 2 and len(prd_xgrs) == 2:
        rct_fmls = list(map(graph_convert.formula, rct_xgrs))
        prd_fmls = list(map(graph_convert.formula, prd_xgrs))

        ret = formula.reac.argsort_hydrogen_abstraction(rct_fmls, prd_fmls)
        if ret:
            rct_idxs_, prd_idxs_ = ret

            q1h_xgr, q2_xgr = list(map(rct_xgrs.__getitem__, rct_idxs_))
            q2h_xgr, q1_xgr = list(map(prd_xgrs.__getitem__, prd_idxs_))

            q1_tras = _partial_hydrogen_abstraction(q1h_xgr, q1_xgr)
            q2_rev_tras = _partial_hydrogen_abstraction(q2h_xgr, q2_xgr)
            if q1_tras and q2_rev_tras:
                for q1_tra, q2_rev_tra in itertools.product(q1_tras,
                                                            q2_rev_tras):
                    xgr1_ = _union(trans.apply(q1_tra, q1h_xgr), q2_xgr)
                    xgr2_ = _union(q1_xgr, q2h_xgr)

                    q2_tra = trans.reverse(q2_rev_tra, xgr2_, xgr1_)
                    tra = trans.from_data(
                        frm_bnd_keys=trans.formed_bond_keys(q2_tra),
                        brk_bnd_keys=trans.broken_bond_keys(q1_tra))

                    tras.append(tra)

                    rct_idxs = rct_idxs_
                    prd_idxs = prd_idxs_

    tras = tuple(tras)
    return tras, rct_idxs, prd_idxs


def _partial_hydrogen_abstraction(qh_xgr, q_xgr):
    tras = []
    h_atm_key = max(_atom_keys(q_xgr)) + 1
    uns_atm_keys = _unsaturated_atom_keys(q_xgr)
    for atm_key in uns_atm_keys:
        q_xgr_h = _add_atom_explicit_hydrogen_keys(
            q_xgr, {atm_key: [h_atm_key]})
        inv_atm_key_dct = _full_isomorphism(q_xgr_h, qh_xgr)
        if inv_atm_key_dct:
            brk_bnd_keys = [frozenset(
                {inv_atm_key_dct[atm_key], inv_atm_key_dct[h_atm_key]})]
            tra = trans.from_data(frm_bnd_keys=[], brk_bnd_keys=brk_bnd_keys)
            tras.append(tra)

    tras = tuple(tras)
    return tras


def addition(rct_xgrs, prd_xgrs):
    """ find an addition transformation

    Additions are identified by joining an unsaturated site on one reactant to
    an unsaturated site on the other. If the result matches the products, this
    is an addition reaction.
    """
    assert _is_valid_reagent_graph_list(rct_xgrs)
    assert _is_valid_reagent_graph_list(prd_xgrs)

    tras = []
    rct_idxs = None
    prd_idxs = None

    if len(rct_xgrs) == 2 and len(prd_xgrs) == 1:
        x_xgr, y_xgr = rct_xgrs
        prd_xgr, = prd_xgrs
        x_atm_keys = _unsaturated_atom_keys(x_xgr)
        y_atm_keys = _unsaturated_atom_keys(y_xgr)
        for x_atm_key, y_atm_key in itertools.product(x_atm_keys, y_atm_keys):
            xy_xgr = _add_bonds(
                _union(x_xgr, y_xgr), [{x_atm_key, y_atm_key}])

            atm_key_dct = _full_isomorphism(xy_xgr, prd_xgr)
            if atm_key_dct:
                tra = trans.from_data(frm_bnd_keys=[{x_atm_key, y_atm_key}],
                                      brk_bnd_keys=[])
                tras.append(tra)

                # sort the reactants so that the largest species is first
                rct_idxs = _argsort_reactants(rct_xgrs)
                prd_idxs = (0,)

    tras = tuple(tras)
    return tras, rct_idxs, prd_idxs


def beta_scission(rct_xgrs, prd_xgrs):
    """ find a beta scission transformation

    Implemented as the reverse of an addition reaction.
    """
    tras = []

    rev_tras, prd_idxs, rct_idxs = addition(prd_xgrs, rct_xgrs)
    if rev_tras:
        rct_xgr = _union_from_sequence(rct_xgrs)
        prd_xgr = _union_from_sequence(prd_xgrs)
        tras = [trans.reverse(tra, prd_xgr, rct_xgr) for tra in rev_tras]

    tras = tuple(set(tras))
    return tras, rct_idxs, prd_idxs


def elimination(rct_xgrs, prd_xgrs):
    """ find an elimination transformation

    Eliminations are identified by breaking two bonds from the reactant,
    forming three fragments. This will form one "central fragment" with two
    break sites and two "end fragments" with one break site each. If the
    central fragment plus the two end fragments, joined at their break sites,
    matches the products, this is an elimination reaction.
    """
    assert _is_valid_reagent_graph_list(rct_xgrs)
    assert _is_valid_reagent_graph_list(prd_xgrs)

    tras = []
    rct_idxs = None
    prd_idxs = None

    if len(rct_xgrs) == 1 and len(prd_xgrs) == 2:
        rct_xgr, = rct_xgrs
        rct_bnd_keys = _bond_keys(rct_xgr)
        # Loop over pairs of bonds and break them. Then, if this forms three
        # fragments, join the two end fragments and compare the result to the
        # products.
        for brk_bnd_key1, brk_bnd_key2 in itertools.combinations(
                rct_bnd_keys, r=2):
            rct_xgr_ = _remove_bonds(rct_xgr, [brk_bnd_key1, brk_bnd_key2])

            # Find the central fragment, which is the one connected to both
            # break sites. If there's a loop there may not be a central
            # fragment, in which case this function will return None.
            cent_frag_atm_keys = _central_fragment_atom_keys(
                rct_xgr_, brk_bnd_key1, brk_bnd_key2)
            if cent_frag_atm_keys is not None:
                atm1_key, = brk_bnd_key1 - cent_frag_atm_keys
                atm2_key, = brk_bnd_key2 - cent_frag_atm_keys
                frm_bnd_key = frozenset({atm1_key, atm2_key})

                rct_xgr_ = _add_bonds(rct_xgr_, [frm_bnd_key])

                prd_xgr = _union_from_sequence(prd_xgrs)
                atm_key_dct = _full_isomorphism(rct_xgr_, prd_xgr)
                if atm_key_dct:
                    tra = trans.from_data(
                        frm_bnd_keys=[frm_bnd_key],
                        brk_bnd_keys=[brk_bnd_key1, brk_bnd_key2])
                    tras.append(tra)

                    rct_idxs = (0,)

                    cent_prd_atm_keys = frozenset(
                        map(atm_key_dct.__getitem__, cent_frag_atm_keys))

                    if cent_prd_atm_keys <= _atom_keys(prd_xgrs[0]):
                        prd_idxs = (0, 1)
                    else:
                        assert cent_prd_atm_keys <= _atom_keys(prd_xgrs[1])
                        prd_idxs = (1, 0)

    tras = tuple(tras)
    return tras, rct_idxs, prd_idxs


def _central_fragment_atom_keys(xgr, brk_bnd_key1, brk_bnd_key2):
    """ Determine atom keys for the central fragment after breaking two bonds.

    The central fragment is the one connected to both break sites.  If there's
    a loop there may not be a central fragment, in which case this function
    will return None.
    """
    xgrs = _connected_components(xgr)
    atm_keys = None
    if len(xgrs) == 3:
        for atm_keys_ in map(_atom_keys, xgrs):
            if (len(brk_bnd_key1 - atm_keys_) == 1 and
                    len(brk_bnd_key2 - atm_keys_) == 1):
                atm_keys = atm_keys_
    return atm_keys


def insertion(rct_xgrs, prd_xgrs):
    """ find a insertion transformation

    Implemented as the reverse of an elimination transformation.
    """
    tras = []

    rev_tras, prd_idxs, rct_idxs = elimination(prd_xgrs, rct_xgrs)
    if rev_tras:
        rct_xgr = _union_from_sequence(rct_xgrs)
        prd_xgr = _union_from_sequence(prd_xgrs)
        tras = [trans.reverse(tra, prd_xgr, rct_xgr) for tra in rev_tras]

    tras = tuple(set(tras))
    return tras, rct_idxs, prd_idxs


def substitution(rct_xgrs, prd_xgrs):
    """ find an substitution transformation

    Substitutions are identified by breaking one bond in the reactants and one
    bond from the products and checking for isomorphism.
    """
    assert _is_valid_reagent_graph_list(rct_xgrs)
    assert _is_valid_reagent_graph_list(prd_xgrs)

    tras = []
    rct_idxs = None
    prd_idxs = None

    if len(rct_xgrs) == 2 and len(prd_xgrs) == 2:
        rct_xgr = _union_from_sequence(rct_xgrs)
        prd_xgr = _union_from_sequence(prd_xgrs)

        rct_bnd_keys = _bond_keys(rct_xgr)
        prd_bnd_keys = _bond_keys(prd_xgr)
        for rct_bnd_key, prd_bnd_key in itertools.product(
                rct_bnd_keys, prd_bnd_keys):
            rct_xgr_ = _remove_bonds(rct_xgr, [rct_bnd_key])
            prd_xgr_ = _remove_bonds(prd_xgr, [prd_bnd_key])

            inv_atm_key_dct = _full_isomorphism(prd_xgr_, rct_xgr_)
            if inv_atm_key_dct:
                brk_bnd_key = rct_bnd_key
                frm_bnd_key = frozenset(
                    map(inv_atm_key_dct.__getitem__, prd_bnd_key))

                tra = trans.from_data(
                    frm_bnd_keys=[frm_bnd_key],
                    brk_bnd_keys=[brk_bnd_key])
                tras.append(tra)

                rct_idxs = _argsort_reactants(rct_xgrs)
                prd_idxs = _argsort_reactants(prd_xgrs)

    tras = tuple(set(tras))
    return tras, rct_idxs, prd_idxs


def _is_valid_reagent_graph_list(xgrs):
    return _are_all_explicit(xgrs) and _have_no_common_atom_keys(xgrs)


def _are_all_explicit(xgrs):
    return all(xgr == _explicit(xgr) for xgr in xgrs)


def _have_no_common_atom_keys(xgrs):
    atm_keys = list(itertools.chain(*map(_atom_keys, xgrs)))
    return len(atm_keys) == len(set(atm_keys))


def _argsort_reactants(xgrs):

    def __sort_value(args):
        _, xgr = args
        val = (-_heavy_atom_count(xgr),
               -_atom_count(xgr),
               -_electron_count(xgr))
        return val

    idxs = tuple(idx for idx, xgr in sorted(enumerate(xgrs), key=__sort_value))
    return idxs
