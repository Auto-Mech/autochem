""" reaction finders

Function arguments:
    Each function takes a list of reactant graphs and a list of product graphs.
    Note that the reactant graphs *cannot* have overlapping atom keys, and
    likewise for the product graphs. Otherwise, there would be no way to
    express the bonds broken and formed between reactants.

Reaction ID Classes:

1 Reac, 1 Prod:
    Hydrogen Migration: 1 frm, 1 brk
2 Reac, 1 Prod:
    Addition: 1 frm
    Addition: 2 frm (NEED) (CH2OO + SO2 = CH2SO3)
    Insertion: 2 frm, 1 brk (BROKEN)
1 Reac, 2 Prod:
    Beta-Scission: 1 brk
    Ring-Forming-Scisson: 1 frm, 1 brk
    Eliminations: 1 frm, 2 brk
2 Reac, 2 Prod
    Hydrogen Abstraction: 1 frm, 1 brk
    Substitutions: 1 frm, 1 brk
    (?) Catalyzed Isom: 2 frm, 2 brk (NEED) (CH2OO + HCOOH = CHOOH + HCOOH)
    ^ double Hydrogen abstraction
"""

import itertools

from .. import form, graph
from ..const import ReactionClass
from ..graph import ts
from ._0core import (
    from_forward_reverse,
    reverse_without_recalculating,
    unique,
)
from ._1util import assert_is_valid_reagent_graph_list
from ._2stereo import expand_stereo_to_match_reagents


def trivial(rct_gras, prd_gras):
    """find a trivial reaction, with the same reactants and products"""
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == len(prd_gras):
        rct_idxs = []
        prd_idxs = []

        # One at a time, find matches for each reactant; track the positions to
        # get the right sort order
        prd_gras_pool = list(prd_gras)
        for rct_idx, rct_gra in enumerate(rct_gras):
            prd_idx = next(
                (
                    idx
                    for idx, prd_gra in enumerate(prd_gras_pool)
                    if graph.isomorphism(rct_gra, prd_gra)
                ),
                None,
            )

            if prd_idx is not None:
                rct_idxs.append(rct_idx)
                prd_idxs.append(prd_idx)
                prd_gras_pool.pop(prd_idx)
            else:
                break

        if rct_idxs and prd_idxs:
            # reorder the reactants and products
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs))

            rcts_gra = graph.union_from_sequence(rct_gras)
            prds_gra = graph.union_from_sequence(prd_gras)

            rxn = from_forward_reverse(
                cla=ReactionClass.TRIVIAL,
                ftsg=ts.graph(rcts_gra, [], []),
                rtsg=ts.graph(prds_gra, [], []),
                rcts_keys=list(map(graph.atom_keys, rct_gras)),
                prds_keys=list(map(graph.atom_keys, prd_gras)),
            )
            rxns.append(rxn)

    return tuple(rxns)


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migrations(rct_gras, prd_gras):
    """find hydrogen migrations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen migrations are identified by adding a hydrogen to an unsaturated
    site of the reactant and adding a hydrogen to an unsaturated site of the
    product and seeing if they match up. If so, we have a hydrogen migration
    between these two sites.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 1:
        (rct_gra,) = rct_gras
        (prd_gra,) = prd_gras

        # Find keys for reactant graph
        rct_h_key = max(graph.atom_keys(rct_gra)) + 1
        rct_rad_keys = graph.unsaturated_atom_keys(rct_gra)

        # Find keys for product graph
        prd_h_key = max(graph.atom_keys(prd_gra)) + 1
        prd_rad_keys = graph.unsaturated_atom_keys(prd_gra)

        for rct_rad_key, prd_rad_key in itertools.product(rct_rad_keys, prd_rad_keys):
            # Add hydrogens to each radical site and see if the result matches
            rct_h_gra = graph.add_bonded_atom(
                rct_gra, "H", rct_rad_key, bnd_atm_key=rct_h_key
            )
            prd_h_gra = graph.add_bonded_atom(
                prd_gra, "H", prd_rad_key, bnd_atm_key=prd_h_key
            )

            iso_dct = graph.isomorphism(rct_h_gra, prd_h_gra, stereo=False)
            if iso_dct:
                inv_dct = dict(map(reversed, iso_dct.items()))

                rct_don_key = inv_dct[prd_rad_key]
                prd_don_key = iso_dct[rct_rad_key]

                # Check equivalent donor atoms for other possible TSs
                rct_don_keys = graph.equivalent_atoms(rct_h_gra, rct_don_key)
                prd_don_keys = graph.equivalent_atoms(prd_h_gra, prd_don_key)

                for rct_don_key, prd_don_key in itertools.product(
                    rct_don_keys, prd_don_keys
                ):
                    # Check equivalent hydrogen atoms on these donors
                    rct_hyd_keys = graph.atom_neighbor_atom_keys(
                        rct_gra, rct_don_key, symb="H"
                    )
                    prd_hyd_keys = graph.atom_neighbor_atom_keys(
                        prd_gra, prd_don_key, symb="H"
                    )

                    for rct_hyd_key, prd_hyd_key in itertools.product(
                        rct_hyd_keys, prd_hyd_keys
                    ):
                        forw_tsg = ts.graph(
                            rct_gra,
                            frm_bnd_keys=[(rct_rad_key, rct_hyd_key)],
                            brk_bnd_keys=[(rct_don_key, rct_hyd_key)],
                        )

                        back_tsg = ts.graph(
                            prd_gra,
                            frm_bnd_keys=[(prd_rad_key, prd_hyd_key)],
                            brk_bnd_keys=[(prd_don_key, prd_hyd_key)],
                        )

                        if graph.isomorphism(
                            forw_tsg, ts.reverse(back_tsg), stereo=False
                        ):
                            # Here, find TSs with stereochemistry that are
                            # consistent
                            rxn = from_forward_reverse(
                                cla=ReactionClass.HYDROGEN_MIGRATION,
                                ftsg=forw_tsg,
                                rtsg=back_tsg,
                                rcts_keys=[graph.atom_keys(rct_gra)],
                                prds_keys=[graph.atom_keys(prd_gra)],
                            )
                            rxns.append(rxn)

    return rxns


# 2. Beta scissions
def beta_scissions(rct_gras, prd_gras):
    """find beta scission reactions

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of additions.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = tuple(map(reverse_without_recalculating, additions(prd_gras, rct_gras)))
    return rxns


# 3. Ring-forming scissions
def ring_forming_scissions(rct_gras, prd_gras):
    """find ring-forming scissions consistent with these reactants and
        products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Ring-forming scissions are found by breaking ring-bonds on one product and
    joining the ends to unsaturated sites on the other product
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        (rgra,) = rct_gras
        pgra = graph.union_from_sequence(prd_gras)
        for pgra1, pgra2 in itertools.permutations(prd_gras):
            bnd_keys = list(itertools.chain(*graph.rings_bond_keys(pgra1)))
            atm_keys = graph.unsaturated_atom_keys(pgra2)

            for bnd_key, atm_key in itertools.product(bnd_keys, atm_keys):
                # Break a ring bond
                gra = graph.remove_bonds(pgra, [bnd_key])

                for end_key in bnd_key:
                    # Add to one end of the broken ring
                    fgra = graph.add_bonds(gra, [(atm_key, end_key)])
                    inv_dct = graph.isomorphism(fgra, rgra, stereo=False)
                    if inv_dct:
                        (other_end_key,) = bnd_key - {end_key}
                        f_frm_bnd_key = (inv_dct[end_key], inv_dct[other_end_key])
                        f_brk_bnd_key = (inv_dct[end_key], inv_dct[atm_key])
                        b_frm_bnd_key = (end_key, atm_key)
                        b_brk_bnd_key = (end_key, other_end_key)
                        forw_tsg = ts.graph(
                            rgra,
                            frm_bnd_keys=[f_frm_bnd_key],
                            brk_bnd_keys=[f_brk_bnd_key],
                        )
                        back_tsg = ts.graph(
                            pgra,
                            frm_bnd_keys=[b_frm_bnd_key],
                            brk_bnd_keys=[b_brk_bnd_key],
                        )

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.RING_FORM_SCISSION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=[graph.atom_keys(rgra)],
                            prds_keys=[graph.atom_keys(pgra1), graph.atom_keys(pgra2)],
                        )
                        rxns.append(rxn)

    return rxns


# 4. Eliminations
def eliminations(rct_gras, prd_gras):
    """find eliminations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Eliminations are identified by forming a bond between an attacking heavy
    atom and another atom not initially bonded to it, forming a ring. The bond
    adjacent to the attacked atom is then broken, along with a second bond in
    the ring, downstream of the attacking heavy atom, away from the attacked
    atom.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        (rct_gra,) = rct_gras
        prds_gra = graph.union_from_sequence(prd_gras)

        # Generate keys all bonds and 1/2 the forming bond
        frm1_keys = graph.atom_keys(rct_gra)
        bnd_keys = graph.bond_keys(rct_gra)

        rct_symbs = graph.atom_symbols(rct_gra)
        frm2_keys = graph.unsaturated_atom_keys(rct_gra)
        frm2_keys = frozenset(key for key in frm2_keys if rct_symbs[key] == "O")

        frm_bnd_keys = [
            (frm1_key, frm2_key)
            for frm1_key, frm2_key in itertools.product(frm1_keys, frm2_keys)
            if frm1_key != frm2_key and frozenset({frm1_key, frm2_key}) not in bnd_keys
        ]

        for frm1_key, frm2_key in frm_bnd_keys:
            prds_gra_ = graph.add_bonds(rct_gra, [(frm2_key, frm1_key)])

            # Get keys of all bonds in the ring formed by this extra bond
            rng_bnd_keys = next(
                (
                    ks
                    for ks in graph.rings_bond_keys(prds_gra_)
                    if frozenset({frm1_key, frm2_key}) in ks
                ),
                None,
            )

            if rng_bnd_keys is not None:
                # Elims break two bonds of the ring formed by the forming bond
                # Loop over all ring bond-pairs, break bonds, see if prods form
                # Ensure to preclude the forming-bond from this set
                brk_bnds = tuple(
                    bond
                    for bond in itertools.combinations(rng_bnd_keys, 2)
                    if frozenset({frm1_key, frm2_key}) not in bond
                )

                for brk_bnd_1, brk_bnd_2 in brk_bnds:
                    prds_gra_2_ = prds_gra_
                    prds_gra_2_ = graph.remove_bonds(prds_gra_2_, [brk_bnd_1])
                    prds_gra_2_ = graph.remove_bonds(prds_gra_2_, [brk_bnd_2])

                    inv_dct = graph.isomorphism(prds_gra_2_, prds_gra, stereo=False)
                    if inv_dct:
                        f_frm_bnd_key = (frm1_key, frm2_key)

                        inv_ = inv_dct.__getitem__
                        b_frm_bnd_key1 = tuple(map(inv_, brk_bnd_1))
                        b_frm_bnd_key2 = tuple(map(inv_, brk_bnd_2))
                        b_brk_bnd_key = tuple(map(inv_, f_frm_bnd_key))

                        forw_tsg = ts.graph(
                            rct_gra,
                            frm_bnd_keys=[f_frm_bnd_key],
                            brk_bnd_keys=[brk_bnd_1, brk_bnd_2],
                        )
                        back_tsg = ts.graph(
                            prds_gra,
                            frm_bnd_keys=[b_frm_bnd_key1, b_frm_bnd_key2],
                            brk_bnd_keys=[b_brk_bnd_key],
                        )

                        rcts_atm_keys = list(map(graph.atom_keys, rct_gras))
                        prds_atm_keys = list(map(graph.atom_keys, prd_gras))

                        if inv_dct[frm1_key] not in prds_atm_keys[1]:
                            prds_atm_keys = list(reversed(prds_atm_keys))

                        assert inv_dct[frm1_key] in prds_atm_keys[1]
                        assert inv_dct[frm2_key] in prds_atm_keys[1]

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.ELIMINATION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        )
                        rxns.append(rxn)

    return rxns


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstractions(rct_gras, prd_gras):
    """find hydrogen abstractions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen abstractions are identified first by checking whether the
    molecular formulas are consistent with a reaction of the form R1H + R2 =>
    R2H + R1. If they do, we identify the abstraction sites by adding hydrogens
    to unsaturated sites of the R1 product to see if we get the R1H reactant.
    We then do the same for the R2 reactant and the R2H product.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_fmls = list(map(graph.formula, rct_gras))
        prd_fmls = list(map(graph.formula, prd_gras))

        ret = form.reac.argsort_hydrogen_abstraction(rct_fmls, prd_fmls)
        if ret:
            rct_idxs_, prd_idxs_ = ret
            rct_gras = list(map(rct_gras.__getitem__, rct_idxs_))
            prd_gras = list(map(prd_gras.__getitem__, prd_idxs_))

            q1h_gra, q2_gra = rct_gras
            q2h_gra, q1_gra = prd_gras

            rets1 = _partial_hydrogen_abstraction(q1h_gra, q1_gra)
            rets2 = _partial_hydrogen_abstraction(q2h_gra, q2_gra)
            for ret1, ret2 in itertools.product(rets1, rets2):
                f_q1h_q_atm_key, f_q1h_h_atm_key, b_q2_q_atm_key = ret1
                b_q1h_q_atm_key, b_q1h_h_atm_key, f_q2_q_atm_key = ret2

                # Create the forward/backward ts graphs
                rcts_gra = graph.union_from_sequence(rct_gras)
                prds_gra = graph.union_from_sequence(prd_gras)
                f_frm_bnd_key = (f_q2_q_atm_key, f_q1h_h_atm_key)
                f_brk_bnd_key = (f_q1h_q_atm_key, f_q1h_h_atm_key)
                b_frm_bnd_key = (b_q2_q_atm_key, b_q1h_h_atm_key)
                b_brk_bnd_key = (b_q1h_q_atm_key, b_q1h_h_atm_key)
                forw_tsg = ts.graph(
                    rcts_gra, frm_bnd_keys=[f_frm_bnd_key], brk_bnd_keys=[f_brk_bnd_key]
                )
                back_tsg = ts.graph(
                    prds_gra, frm_bnd_keys=[b_frm_bnd_key], brk_bnd_keys=[b_brk_bnd_key]
                )

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.HYDROGEN_ABSTRACTION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(graph.atom_keys, rct_gras)),
                    prds_keys=list(map(graph.atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


# 2. Additions
def additions(rct_gras, prd_gras):
    """find additions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Additions are identified by joining an unsaturated site on one reactant to
    an unsaturated site on the other. If the result matches the products, this
    is an addition reaction.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)
    rcts_gra = graph.union_from_sequence(rct_gras)
    prds_gra = graph.union_from_sequence(prd_gras)

    rxns = []

    frm_bnd_pairs = set()
    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct1_gra, rct2_gra = rct_gras
        atm1_keys = graph.unsaturated_atom_keys(rct1_gra)
        atm2_keys = graph.unsaturated_atom_keys(rct2_gra)
        frm_bnd_pairs = set(map(frozenset, itertools.product(atm1_keys, atm2_keys)))
    elif len(rct_gras) == 1 and len(prd_gras) == 1:
        atm_keys = graph.unsaturated_atom_keys(rcts_gra)
        frm_bnd_pairs = set(map(frozenset, itertools.combinations(atm_keys, 2)))
        frm_bnd_pairs -= graph.bond_keys(rcts_gra)

    for atm1_key, atm2_key in frm_bnd_pairs:
        rcts_gra_ = graph.add_bonds(rcts_gra, [(atm1_key, atm2_key)])

        iso_dct = graph.isomorphism(rcts_gra_, prds_gra, stereo=False)
        if iso_dct:
            f_frm_bnd_key = (atm1_key, atm2_key)
            b_brk_bnd_key0 = (iso_dct[atm1_key], iso_dct[atm2_key])
            b_brk_bnd_keys = graph.equivalent_bonds(
                prds_gra, b_brk_bnd_key0, stereo=False, dummy=False
            )

            for b_brk_bnd_key in b_brk_bnd_keys:
                forw_tsg = ts.graph(
                    rcts_gra, frm_bnd_keys=[f_frm_bnd_key], brk_bnd_keys=[]
                )
                back_tsg = ts.graph(
                    prds_gra, frm_bnd_keys=[], brk_bnd_keys=[b_brk_bnd_key]
                )

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.ADDITION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(graph.atom_keys, rct_gras)),
                    prds_keys=list(map(graph.atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


def double_insertion(rct_gras, prd_gras):
    """two atoms inserting

    DOES NOT WORK and breaks regular insertions

    Commit message says this is supposed to handle HOOOH => H2O + O2, but it ***does not
    find anything for this reaction.***
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct_gras = graph.sort_by_size(rct_gras)
        (prd_gra,) = prd_gras
        x_gra, y_gra = rct_gras
        x_atm_keys = graph.unsaturated_atom_keys(x_gra)

        if not len(x_atm_keys) > 1:
            x_gra, y_gra = y_gra, x_gra
            x_atm_keys = graph.unsaturated_atom_keys(x_gra)
        if len(x_atm_keys) > 1:
            brk_bnd_candidates = graph.bond_keys(y_gra)
            for brk_bnd_key in brk_bnd_candidates:
                y_gra_tmp = graph.remove_bonds(y_gra, (brk_bnd_key,))
                permut = itertools.permutations(x_atm_keys, 2)
                frm_bnd_pairs_lst = []
                for comb in permut:
                    zipped = zip(comb, brk_bnd_key)
                    frm_bnd_pairs_lst.append(list(zipped))
                for frm_bnd_pairs in frm_bnd_pairs_lst:
                    xy_gra = graph.add_bonds(
                        graph.union(x_gra, y_gra_tmp), frm_bnd_pairs
                    )

                    iso_dct = graph.isomorphism(xy_gra, prd_gra, stereo=False)
                    if iso_dct:
                        rcts_gra = graph.union_from_sequence(rct_gras)
                        prds_gra = prd_gra
                        brk_atm_i, brk_atm_j = brk_bnd_key
                        b_brk_bnd_key_i = (
                            iso_dct[frm_bnd_pairs[0][0]],
                            iso_dct[frm_bnd_pairs[0][1]],
                        )
                        b_brk_bnd_key_j = (
                            iso_dct[frm_bnd_pairs[1][0]],
                            iso_dct[frm_bnd_pairs[1][1]],
                        )
                        b_frm_bnd_key = (iso_dct[brk_atm_i], iso_dct[brk_atm_j])
                        forw_tsg = ts.graph(
                            rcts_gra,
                            frm_bnd_keys=frm_bnd_pairs,
                            brk_bnd_keys=(brk_bnd_key,),
                        )
                        back_tsg = ts.graph(
                            prds_gra,
                            brk_bnd_keys=[b_brk_bnd_key_i, b_brk_bnd_key_j],
                            frm_bnd_keys=(b_frm_bnd_key,),
                        )
                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.DOUBLE_INSERTION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=list(map(graph.atom_keys, rct_gras)),
                            prds_keys=list(map(graph.atom_keys, prd_gras)),
                        )
                        rxns.append(rxn)

    return rxns


def two_bond_additions(rct_gras, prd_gras):
    """two bond additions"""
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct_gras = graph.sort_by_size(rct_gras)
        x_gra, y_gra = rct_gras
        (prd_gra,) = prd_gras
        x_atm_keys = frozenset().union(
            graph.unsaturated_atom_keys(x_gra), graph.lone_pair_atom_keys(x_gra)
        )
        y_atm_keys = frozenset().union(
            graph.unsaturated_atom_keys(y_gra), graph.lone_pair_atom_keys(y_gra)
        )

        # Generate pairs of forming bonds, where each is a pair of idxs
        # describing the atoms making up the forming bond:
        # (frm1, frm2) = ((idx1, idx2), (idx1, idx2))
        frm_bnd_pairs = tuple(itertools.product(x_atm_keys, y_atm_keys))
        frm_bnds_lst = ()
        for pair in itertools.product(frm_bnd_pairs, frm_bnd_pairs):
            # Preclude pairs with same idxs (formind same bond twice)
            if pair[0] != pair[1]:
                # Preclude multiple bonds formed to same atom X---A---Y
                if pair[0][0] != pair[1][0] and pair[0][1] != pair[1][1]:
                    # Preclude the reverse
                    if pair[::-1] not in frm_bnds_lst:
                        frm_bnds_lst += (pair,)

        for frm_bnd_keys in frm_bnds_lst:
            xy_gra = graph.add_bonds(
                graph.union(x_gra, y_gra), [set(frm_bnd_keys[0]), set(frm_bnd_keys[1])]
            )

            iso_dct = graph.isomorphism(xy_gra, prd_gra, stereo=False)
            if iso_dct:
                rcts_gra = graph.union_from_sequence(rct_gras)
                prds_gra = prd_gra
                b_brk_bnd_keys = [
                    [iso_dct[frm_bnd_keys[0][0]], iso_dct[frm_bnd_keys[0][1]]],
                    [iso_dct[frm_bnd_keys[1][0]], iso_dct[frm_bnd_keys[1][1]]],
                ]
                forw_tsg = ts.graph(
                    rcts_gra, frm_bnd_keys=frm_bnd_keys, brk_bnd_keys=[]
                )
                back_tsg = ts.graph(
                    prds_gra, frm_bnd_keys=[], brk_bnd_keys=b_brk_bnd_keys
                )

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.ADDITION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(graph.atom_keys, rct_gras)),
                    prds_keys=list(map(graph.atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


# 3. Insertions
def insertions(rct_gras, prd_gras):
    """find insertions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of an elimination reaction.

    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = tuple(map(reverse_without_recalculating, eliminations(prd_gras, rct_gras)))
    return rxns


# 4. Substitutions
def substitutions(rct_gras, prd_gras):
    """find substitutions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Substitutions are identified by breaking one bond in the reactants and one
    bond from the products and checking for isomorphism.

    Currently it assumes that one of the reactants has a radical site that
    can attack the other reactants, forming a bond and breaking another.

    From the perspective of breaking and forming breaking bonds, substitutions
    are equivalent with hydrogen abstractions. Hence, we remove all cases where
    the forming bond involves a hydrogen atom off the reactant in which a bond
    is breaking.
    """
    rct_gras, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_gra = graph.union_from_sequence(rct_gras)
        prd_gra = graph.union_from_sequence(prd_gras)

        # Loop over both orders of reactants: A+B and B+A
        for rgra1, rgra2 in itertools.permutations(rct_gras):
            tet_dct = graph.tetrahedral_atoms(rgra1, min_ncount=0)
            rad_keys = graph.unsaturated_atom_keys(rgra2)

            for (key, nkeys), rad_key in itertools.product(tet_dct.items(), rad_keys):
                # Form a bond to the radical atom
                f_frm_bkey = (key, rad_key)
                gra = graph.add_bonds(rct_gra, [f_frm_bkey])
                for nkey in nkeys:
                    # Remove a bond to the neighbor atom
                    f_brk_bkey = (key, nkey)
                    gra_ = graph.remove_bonds(gra, [f_brk_bkey])

                    inv_dct = graph.isomorphism(gra_, prd_gra, stereo=False)
                    if inv_dct:
                        b_brk_bkey = tuple(map(inv_dct.get, f_frm_bkey))
                        b_frm_bkey = tuple(map(inv_dct.get, f_brk_bkey))

                        forw_tsg = ts.graph(
                            rct_gra,
                            frm_bnd_keys=[f_frm_bkey],
                            brk_bnd_keys=[f_brk_bkey],
                        )
                        back_tsg = ts.graph(
                            prd_gra,
                            frm_bnd_keys=[b_frm_bkey],
                            brk_bnd_keys=[b_brk_bkey],
                        )

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.SUBSTITUTION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=list(map(graph.atom_keys, rct_gras)),
                            prds_keys=list(map(graph.atom_keys, prd_gras)),
                        )
                        rxns.append(rxn)

    return rxns


def _is_graph(gra):
    """Determine whether the argument is a graph, instead of a sequence"""
    if hasattr(gra, "__len__") and len(gra) == 2:
        if all(isinstance(o, dict) for o in gra):
            return True

    return False


def find(rct_gras, prd_gras, stereo=False):
    """find all reactions consistent with these reactants and products

    :param rct_gras: graphs for the reactants without overlapping keys
    :param prd_gras: graphs for the products without overlapping keys
    :param stereo: Find stereo-specified reaction objects?
    :type stereo: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]
    """
    if _is_graph(rct_gras):
        assert _is_graph(prd_gras)
        rct_gras = graph.connected_components(rct_gras)
        prd_gras = graph.connected_components(prd_gras)

    rct_gras0, _ = graph.standard_keys_for_sequence(rct_gras)
    prd_gras0, _ = graph.standard_keys_for_sequence(prd_gras)

    rct_gras = list(map(graph.without_stereo, rct_gras0))
    prd_gras = list(map(graph.without_stereo, prd_gras0))

    assert_is_valid_reagent_graph_list(rct_gras)
    assert_is_valid_reagent_graph_list(prd_gras)

    # check whether this is a valid reaction
    rct_fmls = list(map(graph.formula, rct_gras))
    prd_fmls = list(map(graph.formula, prd_gras))
    if not form.reac.is_valid_reaction(rct_fmls, prd_fmls):
        return ()

    # Cycle through the different finders and gather all possible reactions
    finders_ = [
        # trivial,
        # unimolecular reactions
        hydrogen_migrations,
        beta_scissions,
        ring_forming_scissions,
        eliminations,
        # bimolecular reactions
        hydrogen_abstractions,
        additions,
        # double_insertion,
        # two_bond_additions,
        insertions,  # not fully functional if elims broken
        substitutions,
    ]

    all_rxns = []
    for finder_ in finders_:
        rxns = finder_(rct_gras, prd_gras)
        for rxn in rxns:
            if not stereo:
                all_rxns.append(rxn)
            else:
                srxns = expand_stereo_to_match_reagents(
                    rxn, rct_gras0, prd_gras0, shift_keys=True
                )
                all_rxns.extend(srxns)
    # Check for uniqueness *after* stereochemistry is assigned
    all_rxns = unique(all_rxns)
    return tuple(all_rxns)


# helpers
def _partial_hydrogen_abstraction(qh_gra, q_gra):
    rets = []

    h_atm_key = max(graph.atom_keys(q_gra)) + 1
    uns_atm_keys = graph.unsaturated_atom_keys(q_gra)
    for atm_key in uns_atm_keys:
        q_gra_h = graph.add_atom_explicit_hydrogens(q_gra, {atm_key: [h_atm_key]})
        inv_atm_key_dct = graph.isomorphism(q_gra_h, qh_gra, stereo=False)
        if inv_atm_key_dct:
            qh_q_atm_key = inv_atm_key_dct[atm_key]
            qh_h_atm_keys = graph.atom_neighbor_atom_keys(
                qh_gra, qh_q_atm_key, symb="H"
            )
            q_q_atm_key = atm_key
            for qh_h_atm_key in qh_h_atm_keys:
                rets.append((qh_q_atm_key, qh_h_atm_key, q_q_atm_key))

    return rets


def arbitrary_reactions(
    rct_gras, prd_gras, max_brk=2, max_frm=2, min_brk=0, min_frm=0, max_trans=3):
    """find all transformations consistent with these reactants and products

    :param rct_gras: graphs for the reactants without overlapping keys
    :param prd_gras: graphs for the products without overlapping keys
    :param max_brk: (int) maximum number of breakable bonds
    :param max_frm: (int) maximum number of forming bonds
    :param min_brk: (int) minimum number of breakable bonds
    :param min_frm: (int) minimum number of forming bonds
    :param max_trans: (int) maximum number of combined transformations
    :returns: the isomorphism dictionary between reactants and products
        alongwith the lists of breaking bond keys and forming keys
    :rtype: tuple[dict, tuple[frozenset], tuple[frozenset]
    """
    def _set_bimol_direction(rct_gras, prd_gras):
        # we can reduce dimensionality of bXfy matrix if it is
        # not unimolecular, by enforcing that there is a forming
        # bond between species of bimoleculars

        bimol_dir = len(rct_gras) - len(prd_gras)
        # swap reaction direction to have more reactants than products
        if bimol_dir < 0:
            prd_gra = graph.union_from_sequence(rct_gras)
            rct_gra = graph.union_from_sequence(prd_gras)
        else:
            rct_gra = graph.union_from_sequence(rct_gras)
            prd_gra = graph.union_from_sequence(prd_gras)

        return bimol_dir, rct_gra, prd_gra

    def _greatest_overlap():
        # largest_common_fragment_isomorphism doesn't always
        # return the best overlap, this function has different
        # constraints
        common_fragments = graph.common_fragments_isomorphisms(rct_gra, prd_gra)
        stable_key_pairs = {}
        stable_bnd_keys = ()
        stable_atm_keys = ()
        stable_rct_gra = ()
        stable_prd_gra = ()
        overlap_perc = -1
        if common_fragments:
            stable_key_pairs = common_fragments[0]
        for common_iso_dct in common_fragments:
            if len(common_iso_dct) < len(stable_key_pairs):
                break
            # # does it have to be the longest, or can we exclude all?
            # subgras = graph.connected_components(
            #     graph.subgraph(rct_gra, stable_key_pairs.keys()))
            # longest_subgra = max(subgras, key=lambda x: len(x[0]))
            longest_subgra = graph.subgraph(rct_gra, common_iso_dct.keys())
            common_atm_keys = graph.atom_keys(longest_subgra)
            common_bnd_keys = graph.bond_keys(longest_subgra)
            # atom keys are only truely common if they have the same bonds
            # in the reactants and products and there is no radical
            true_common_atm_keys = []
            for key in common_atm_keys:
                if key in graph.unsaturated_atom_keys(rct_gra):
                    continue
                adj_atms = graph.atom_neighbor_atom_keys(rct_gra, key)
                if all(adj_atm in common_iso_dct.keys() for adj_atm in adj_atms):
                    prd_adj_atms = [common_iso_dct[adj_atm] for adj_atm in adj_atms]
                    if set(prd_adj_atms) == set(graph.atom_neighbor_atom_keys(
                            prd_gra, common_iso_dct[key])):
                        true_common_atm_keys.append(key)

            # set the hydrogens to explicit so we can break/form them,
            # but not if they are bound to a common atm
            common_rct_gra = graph.explicit(
                rct_gra, atm_keys=set(
                    graph.atom_keys(rct_gra))-set(true_common_atm_keys))
            common_prd_gra = graph.explicit(
                prd_gra,
                atm_keys=set(graph.atom_keys(prd_gra))-set(
                    [common_iso_dct[atm] for atm in true_common_atm_keys]))

            # get the overlapping fragment (the graph from the common atm keys)
            shared_rct_gra = graph.subgraph(rct_gra, atm_keys=set(true_common_atm_keys))
            shared_bnd_keys = graph.bond_keys(shared_rct_gra)
            true_common_bnd_keys = frozenset(
                [key for key in shared_bnd_keys if len(
                    graph.unsaturated_atom_keys(common_rct_gra) & key) == 0])

            cur_overlap_perc = len(common_bnd_keys)/len(
                graph.bond_keys(graph.implicit(common_rct_gra)))
            if overlap_perc < cur_overlap_perc:
                overlap_perc = cur_overlap_perc
                stable_key_pairs = common_iso_dct
                stable_bnd_keys = true_common_bnd_keys
                stable_atm_keys = true_common_atm_keys
                stable_rct_gra = common_rct_gra
                stable_prd_gra = common_prd_gra
        return (
            overlap_perc, stable_key_pairs, stable_bnd_keys,
            stable_atm_keys, stable_rct_gra, stable_prd_gra)

    def _edge_connects_different_groups(frm_bnd):
        # Iterate over all pairs of groups
        atom_sets = [set(group) for group in connected_keys]
        for i, atom_set_i in enumerate(atom_sets):
            for j in range(i + 1, len(atom_sets)):
                if (frm_bnd & atom_set_i) and (frm_bnd & atom_sets[j]):
                    return True
                if (frm_bnd & atom_sets[j]) and (frm_bnd & atom_set_i):
                    return True
        return False

    def _multiple_transitions_per_atom(bnd_lst):
        flat_bnds = [atm for bnd in bnd_lst for atm in bnd]
        return len(flat_bnds) - len(set(flat_bnds))

    def _radical_is_involved(frm_bnd_lst, brk_bnd_lst):
        # include adjacent atoms to radical for beta-scission etc
        # rad_pos_lst = graph.radical_atom_keys(gra, sing_res=False)
        frm_bnd_set_flat = frozenset().union(*frm_bnd_lst)
        brk_bnd_set_flat = frozenset().union(*brk_bnd_lst)
        return (
            len(frm_bnd_set_flat & rad_incl_set) +
            len(brk_bnd_set_flat & rad_incl_set))

    inv_dct = {}
    frm_bnd_combo = ()
    brk_bnd_combo = ()
    bimol_dir, rct_gra, prd_gra = _set_bimol_direction(
        rct_gras, prd_gras)

    # find the overlap between the two species
    (
        overlap_perc,
        stable_key_pairs,
        stable_bnd_keys,
        stable_atm_keys,
        stable_rct_gra,
        stable_prd_gra
    ) = _greatest_overlap()
    print('overlap', overlap_perc)
    print('stable_bnd_keys', stable_bnd_keys)
    print('stable_key_pairs', stable_key_pairs)

    if overlap_perc > 0.3 or len(list(stable_key_pairs.keys())) > 6:
        # Get some info for checks later
        # (1) identify sites near the radical
        # (2) identify connected components of bimol well
        ngb_dct = graph.atoms_neighbor_atom_keys(stable_rct_gra)
        rad_pos_lst = graph.unsaturated_atom_keys(stable_rct_gra)
        rad_incl_set = rad_pos_lst.union(frozenset({
            adj_atm for rad_atm in rad_pos_lst
            for adj_atm in ngb_dct[rad_atm]
            if graph.atom_symbols(stable_rct_gra)[adj_atm] != 'H'}))
        if abs(bimol_dir) != 0:
            connected_keys = graph.connected_components_atom_keys(stable_rct_gra)

        # build list of breakable bonds
        bnd_keys = graph.bond_keys(stable_rct_gra)
        brk_bnd_keys = [
             key for key in bnd_keys
             if key not in stable_bnd_keys]

        # build list of formable bonds
        # ... all possible formations between atoms
        # ... but only including one of each equivalent atom
        # ... where they are only considered equivalent if
        # ... bound to the same atom, aka C1-C2=C3-C4 won't
        # ... have C2 and C3 equiv but will have all Hs on C1
        # ... as equiv, we can adjust it to need just one
        # ... same neighbor if we want to equivalate methyl groups etc
        frm_atm_keys_w_dups = graph.atom_keys(stable_rct_gra)
        frm_atm_keys = ()
        for atm in frm_atm_keys_w_dups:
            atm_ngbs = ngb_dct[atm]
            equiv_atms = [
                equiv_atm for equiv_atm in graph.equivalent_atoms(stable_rct_gra, atm)
                if ngb_dct[equiv_atm] == atm_ngbs]
            if not any(eq_atm in frm_atm_keys for eq_atm in equiv_atms):
                frm_atm_keys += (atm,)
        frm_bnd_keys = [
            frozenset({frm1_key, frm2_key})
            for frm1_key, frm2_key in itertools.combinations(frm_atm_keys, 2)
            if (
                frozenset({frm1_key, frm2_key}) not in bnd_keys
            )]

        # get allowed combinations of NUMBER of breaking/forming bonds
        max_brk = min(max_brk, len(brk_bnd_keys))
        max_frm = min(max_frm, len(frm_bnd_keys))
        brk_frm_combos = itertools.product(
                range(min_brk, max_brk + 1), range(min_frm, max_frm + 1))
        brk_frm_combos = sorted(brk_frm_combos, key=lambda x: sum(x))
        num_checks = 0
        for (num_brk, num_frm) in brk_frm_combos:
            if num_brk + num_frm > max_trans:
                continue
            # loop over forming bond key combinations
            for pot_frm_bnd_combo in itertools.combinations(frm_bnd_keys, num_frm):
                if inv_dct:
                    break

                frm_bnd_combo = pot_frm_bnd_combo
                if _multiple_transitions_per_atom(frm_bnd_combo):
                    continue
                if abs(bimol_dir) != 0 and not any(_edge_connects_different_groups(
                        frm_bnd) for frm_bnd in frm_bnd_combo):
                    continue

                # loop over breaking bond key combinations
                for brk_bnd_combo in itertools.combinations(brk_bnd_keys, num_brk):
                    if _multiple_transitions_per_atom(brk_bnd_combo):
                        continue
                    if not _radical_is_involved(frm_bnd_combo, brk_bnd_combo):
                        continue

                    num_checks += 1
                    pot_rct_gra = graph.remove_bonds(stable_rct_gra, brk_bnd_combo)
                    pot_rct_gra = graph.add_bonds(pot_rct_gra, frm_bnd_combo)
                    inv_dct = graph.isomorphism(
                        pot_rct_gra, stable_prd_gra, stereo=False)
                    if inv_dct:
                        if bimol_dir < 0:
                            # reverse back
                            tmp_frm_bnd_combo =  tuple(frozenset(
                                inv_dct[atm] for atm in bnd) for bnd in brk_bnd_combo)
                            brk_bnd_combo =  tuple(frozenset(
                                inv_dct[atm] for atm in bnd) for bnd in frm_bnd_combo)
                            frm_bnd_combo = tmp_frm_bnd_combo
                            inv_dct = {value: key for key, value in inv_dct.items()}
                        print(
                            'found! breaking',
                             brk_bnd_combo,
                            'forming',
                            frm_bnd_combo)
                        break
    return inv_dct, frm_bnd_combo, brk_bnd_combo
