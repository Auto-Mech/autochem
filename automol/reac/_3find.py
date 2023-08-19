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
import automol.geom
import automol.chi
from automol.par import ReactionClass
from automol.graph import ts
from automol.graph import atom_symbols
from automol.graph import atom_keys
from automol.graph import bond_keys
from automol.graph import formula
from automol.graph import union
from automol.graph import add_bonds
from automol.graph import remove_bonds
from automol.graph import isomorphism
from automol.graph import equivalent_atoms
from automol.graph import union_from_sequence
from automol.graph import unsaturated_atom_keys
from automol.graph import lone_pair_atom_keys
from automol.graph import atom_neighbor_atom_keys
from automol.graph import add_bonded_atom
from automol.graph import add_atom_explicit_hydrogens
from automol.graph import rings_bond_keys
from automol.reac._0core import from_forward_reverse
from automol.reac._0core import reverse
from automol.reac._0core import unique
from automol.reac._1util import assert_is_valid_reagent_graph_list
from automol.reac._1util import sort_reagents
from automol.reac._2stereo import expand_stereo_for_reaction


def trivial(rct_gras, prd_gras):
    """ find a trivial reaction, with the same reactants and products
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == len(prd_gras):
        rct_idxs = []
        prd_idxs = []

        # One at a time, find matches for each reactant; track the positions to
        # get the right sort order
        prd_gras_pool = list(prd_gras)
        for rct_idx, rct_gra in enumerate(rct_gras):
            prd_idx = next((idx for idx, prd_gra in enumerate(prd_gras_pool)
                            if isomorphism(rct_gra, prd_gra)), None)

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

            rcts_gra = union_from_sequence(rct_gras)
            prds_gra = union_from_sequence(prd_gras)

            rxn = from_forward_reverse(
                cla=ReactionClass.Typ.TRIVIAL,
                ftsg=ts.graph(rcts_gra, [], []),
                rtsg=ts.graph(prds_gra, [], []),
                rcts_keys=list(map(atom_keys, rct_gras)),
                prds_keys=list(map(atom_keys, prd_gras)),
            )
            rxns.append(rxn)

    return tuple(rxns)


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migrations(rct_gras, prd_gras):
    """ find hydrogen migrations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen migrations are identified by adding a hydrogen to an unsaturated
    site of the reactant and adding a hydrogen to an unsaturated site of the
    product and seeing if they match up. If so, we have a hydrogen migration
    between these two sites.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 1:
        rct_gra, = rct_gras
        prd_gra, = prd_gras

        # Find keys for reactant graph
        rct_h_key = max(atom_keys(rct_gra)) + 1
        rct_rad_keys = unsaturated_atom_keys(rct_gra)

        # Find keys for product graph
        prd_h_key = max(atom_keys(prd_gra)) + 1
        prd_rad_keys = unsaturated_atom_keys(prd_gra)

        for rct_rad_key, prd_rad_key in (
                itertools.product(rct_rad_keys, prd_rad_keys)):
            # Add hydrogens to each radical site and see if the result matches
            rct_h_gra = add_bonded_atom(
                rct_gra, 'H', rct_rad_key, bnd_atm_key=rct_h_key)
            prd_h_gra = add_bonded_atom(
                prd_gra, 'H', prd_rad_key, bnd_atm_key=prd_h_key)

            iso_dct = isomorphism(rct_h_gra, prd_h_gra, stereo=False)
            if iso_dct:
                inv_dct = dict(map(reversed, iso_dct.items()))

                rct_don_key = inv_dct[prd_rad_key]
                prd_don_key = iso_dct[rct_rad_key]

                # Check equivalent donor atoms for other possible TSs
                rct_don_keys = equivalent_atoms(rct_h_gra, rct_don_key)
                prd_don_keys = equivalent_atoms(prd_h_gra, prd_don_key)

                for rct_don_key, prd_don_key in (
                        itertools.product(rct_don_keys, prd_don_keys)):

                    # Check equivalent hydrogen atoms on these donors
                    rct_hyd_keys = atom_neighbor_atom_keys(
                        rct_gra, rct_don_key, symb='H')
                    prd_hyd_keys = atom_neighbor_atom_keys(
                        prd_gra, prd_don_key, symb='H')

                    for rct_hyd_key, prd_hyd_key in (
                            itertools.product(rct_hyd_keys, prd_hyd_keys)):
                        forw_tsg = ts.graph(
                            rct_gra,
                            frm_bnd_keys=[(rct_rad_key, rct_hyd_key)],
                            brk_bnd_keys=[(rct_don_key, rct_hyd_key)])

                        back_tsg = ts.graph(
                            prd_gra,
                            frm_bnd_keys=[(prd_rad_key, prd_hyd_key)],
                            brk_bnd_keys=[(prd_don_key, prd_hyd_key)])

                        if isomorphism(forw_tsg, ts.reverse(back_tsg),
                                       stereo=False):
                            # Here, find TSs with stereochemistry that are
                            # consistent
                            rxn = from_forward_reverse(
                                cla=ReactionClass.Typ.HYDROGEN_MIGRATION,
                                ftsg=forw_tsg,
                                rtsg=back_tsg,
                                rcts_keys=[atom_keys(rct_gra)],
                                prds_keys=[atom_keys(prd_gra)],
                            )
                            rxns.append(rxn)

    return rxns


# 2. Beta scissions
def beta_scissions(rct_gras, prd_gras):
    """ find beta scission reactions

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of additions.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = tuple(map(reverse, additions(prd_gras, rct_gras)))
    return rxns


# 3. Ring-forming scissions
def ring_forming_scissions(rct_gras, prd_gras):
    """ find ring-forming scissions consistent with these reactants and
        products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Ring-forming scissions are found by breaking ring-bonds on one product and
    joining the ends to unsaturated sites on the other product
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        rgra, = rct_gras
        pgra = union_from_sequence(prd_gras)
        for pgra1, pgra2 in itertools.permutations(prd_gras):
            bnd_keys = list(itertools.chain(*rings_bond_keys(pgra1)))
            atm_keys = unsaturated_atom_keys(pgra2)

            for bnd_key, atm_key in itertools.product(bnd_keys, atm_keys):
                # Break a ring bond
                gra = remove_bonds(pgra, [bnd_key])

                for end_key in bnd_key:
                    # Add to one end of the broken ring
                    fgra = add_bonds(gra, [(atm_key, end_key)])
                    inv_dct = isomorphism(fgra, rgra, stereo=False)
                    if inv_dct:
                        other_end_key, = bnd_key - {end_key}
                        f_frm_bnd_key = (inv_dct[end_key],
                                         inv_dct[other_end_key])
                        f_brk_bnd_key = (inv_dct[end_key], inv_dct[atm_key])
                        b_frm_bnd_key = (end_key, atm_key)
                        b_brk_bnd_key = (end_key, other_end_key)
                        forw_tsg = ts.graph(rgra,
                                            frm_bnd_keys=[f_frm_bnd_key],
                                            brk_bnd_keys=[f_brk_bnd_key])
                        back_tsg = ts.graph(pgra,
                                            frm_bnd_keys=[b_frm_bnd_key],
                                            brk_bnd_keys=[b_brk_bnd_key])

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.Typ.RING_FORM_SCISSION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=[atom_keys(rgra)],
                            prds_keys=[atom_keys(pgra1), atom_keys(pgra2)],
                        )
                        rxns.append(rxn)

    return rxns


# 4. Eliminations
def eliminations(rct_gras, prd_gras):
    """ find eliminations consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Eliminations are identified by forming a bond between an attacking heavy
    atom and another atom not initially bonded to it, forming a ring. The bond
    adjacent to the attacked atom is then broken, along with a second bond in
    the ring, downstream of the attacking heavy atom, away from the attacked
    atom.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    def _identify(frm1_keys, frm2_keys, bnd_keys):
        """ Try and identify elmination from some set of keys
        """
        _rxns = []

        frm_bnd_keys = [(frm1_key, frm2_key) for frm1_key, frm2_key
                        in itertools.product(frm1_keys, frm2_keys)
                        if frm1_key != frm2_key and
                        frozenset({frm1_key, frm2_key}) not in bnd_keys]

        for frm1_key, frm2_key in frm_bnd_keys:

            prds_gra_ = add_bonds(rct_gra, [(frm2_key, frm1_key)])

            # Get keys of all bonds in the ring formed by this extra bond
            rng_bnd_keys = next(
                (ks for ks in rings_bond_keys(prds_gra_)
                 if frozenset({frm1_key, frm2_key}) in ks), None)

            if rng_bnd_keys is not None:

                # Elims break two bonds of the ring formed by the forming bond
                # Loop over all ring bond-pairs, break bonds, see if prods form
                # Ensure to preclude the forming-bond from this set
                brk_bnds = tuple(
                    bond for bond in itertools.combinations(rng_bnd_keys, 2)
                    if frozenset({frm1_key, frm2_key}) not in bond)

                for brk_bnd_1, brk_bnd_2 in brk_bnds:
                    prds_gra_2_ = prds_gra_
                    prds_gra_2_ = remove_bonds(prds_gra_2_, [brk_bnd_1])
                    prds_gra_2_ = remove_bonds(prds_gra_2_, [brk_bnd_2])

                    inv_dct = isomorphism(prds_gra_2_, prds_gra, stereo=False)
                    if inv_dct:
                        f_frm_bnd_key = (frm1_key, frm2_key)

                        inv_ = inv_dct.__getitem__
                        b_frm_bnd_key1 = tuple(map(inv_, brk_bnd_1))
                        b_frm_bnd_key2 = tuple(map(inv_, brk_bnd_2))
                        b_brk_bnd_key = tuple(map(inv_, f_frm_bnd_key))

                        forw_tsg = ts.graph(rct_gra,
                                            frm_bnd_keys=[f_frm_bnd_key],
                                            brk_bnd_keys=[brk_bnd_1,
                                                          brk_bnd_2])
                        back_tsg = ts.graph(prds_gra,
                                            frm_bnd_keys=[b_frm_bnd_key1,
                                                          b_frm_bnd_key2],
                                            brk_bnd_keys=[b_brk_bnd_key])

                        rcts_atm_keys = list(map(atom_keys, rct_gras))
                        prds_atm_keys = list(map(atom_keys, prd_gras))

                        if inv_dct[frm1_key] not in prds_atm_keys[1]:
                            prds_atm_keys = list(reversed(prds_atm_keys))

                        assert inv_dct[frm1_key] in prds_atm_keys[1]
                        assert inv_dct[frm2_key] in prds_atm_keys[1]

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.Typ.ELIMINATION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        )
                        rxns.append(rxn)

        return _rxns

    rxns = []

    if len(rct_gras) == 1 and len(prd_gras) == 2:
        rct_gra, = rct_gras
        prds_gra = union_from_sequence(prd_gras)

        # ngb_keys_dct = atoms_neighbor_atom_keys(rct_gra)

        # Generate keys all bonds and 1/2 the forming bond
        frm1_keys = atom_keys(rct_gra)
        bnd_keys = bond_keys(rct_gra)

        frm2_keys = unsaturated_atom_keys(rct_gra)
        rct_symbs = atom_symbols(rct_gra)
        frm2_keys_o = frozenset(key for key in frm2_keys
                                if rct_symbs[key] == 'O')

        rxns.extend(_identify(frm1_keys, frm2_keys_o, bnd_keys))

        # OLD WAY. More IDs but more mistakes
        # To make the function general, try to ID reaction
        # with different types of keys for the attacking atom
        # (1) unsaturated atom sites
        # frm2_keys = unsaturated_atom_keys(rct_gra)
        # rxns.extend(_identify(frm1_keys, frm2_keys, bnd_keys))
        # if not rxns:
        #     # (2) remaining saturated atom sites
        #     frm2_keys = atom_keys(rct_gra, excl_syms=('H',)) - frm2_keys
        #     rxns.extend(_identify(frm1_keys, frm2_keys, bnd_keys))
        #     # if not rxns:  # Ignoring H2 formation for now for speed
        #     #     # (3) H atoms
        #     #     frm1_keys = atom_keys(rct_gra, sym='H')
        #     #     rxns.extend(_identify(frm1_keys, frm2_keys, bnd_keys))

    return rxns


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstractions(rct_gras, prd_gras):
    """ find hydrogen abstractions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Hydrogen abstractions are identified first by checking whether the
    molecular formulas are consistent with a reaction of the form R1H + R2 =>
    R2H + R1. If they do, we identify the abstraction sites by adding hydrogens
    to unsaturated sites of the R1 product to see if we get the R1H reactant.
    We then do the same for the R2 reactant and the R2H product.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_fmls = list(map(formula, rct_gras))
        prd_fmls = list(map(formula, prd_gras))

        ret = automol.formula.reac.argsort_hydrogen_abstraction(
            rct_fmls, prd_fmls)
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
                rcts_gra = union_from_sequence(rct_gras)
                prds_gra = union_from_sequence(prd_gras)
                f_frm_bnd_key = (f_q2_q_atm_key, f_q1h_h_atm_key)
                f_brk_bnd_key = (f_q1h_q_atm_key, f_q1h_h_atm_key)
                b_frm_bnd_key = (b_q2_q_atm_key, b_q1h_h_atm_key)
                b_brk_bnd_key = (b_q1h_q_atm_key, b_q1h_h_atm_key)
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[f_brk_bnd_key])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[b_frm_bnd_key],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.Typ.HYDROGEN_ABSTRACTION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


# 2. Additions
def additions(rct_gras, prd_gras):
    """ find additions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Additions are identified by joining an unsaturated site on one reactant to
    an unsaturated site on the other. If the result matches the products, this
    is an addition reaction.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct_gras = sort_reagents(rct_gras)
        x_gra, y_gra = rct_gras
        prd_gra, = prd_gras
        x_atm_keys = unsaturated_atom_keys(x_gra)
        y_atm_keys = unsaturated_atom_keys(y_gra)

        frm_bnd_pairs = tuple(itertools.product(x_atm_keys, y_atm_keys))
        for x_atm_key, y_atm_key in frm_bnd_pairs:
            xy_gra = add_bonds(
                union(x_gra, y_gra), [{x_atm_key, y_atm_key}])

            iso_dct = isomorphism(xy_gra, prd_gra, stereo=False)
            if iso_dct:
                rcts_gra = union_from_sequence(rct_gras)
                prds_gra = prd_gra
                f_frm_bnd_key = (x_atm_key, y_atm_key)
                b_brk_bnd_key = (iso_dct[x_atm_key], iso_dct[y_atm_key])
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=[f_frm_bnd_key],
                                    brk_bnd_keys=[])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[],
                                    brk_bnd_keys=[b_brk_bnd_key])

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.Typ.ADDITION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


def double_insertion(rct_gras, prd_gras):
    """ two atoms inserting

    DOES NOT WORK and breaks regular insertions

    Commit message says this is supposed to handle HOOOH => H2O + O2, but it ***does not
    find anything for this reaction.***
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct_gras = sort_reagents(rct_gras)
        prd_gra, = prd_gras
        x_gra, y_gra = rct_gras
        x_atm_keys = unsaturated_atom_keys(x_gra)

        if not len(x_atm_keys) > 1:
            x_gra, y_gra = y_gra, x_gra
            x_atm_keys = unsaturated_atom_keys(x_gra)
        if len(x_atm_keys) > 1:
            brk_bnd_candidates = bond_keys(y_gra)
            for brk_bnd_key in brk_bnd_candidates:
                y_gra_tmp = remove_bonds(y_gra, (brk_bnd_key,))
                permut = itertools.permutations(x_atm_keys, 2)
                frm_bnd_pairs_lst = []
                for comb in permut:
                    zipped = zip(comb, brk_bnd_key)
                    frm_bnd_pairs_lst.append(list(zipped))
                for frm_bnd_pairs in frm_bnd_pairs_lst:
                    xy_gra = add_bonds(
                        union(x_gra, y_gra_tmp), frm_bnd_pairs)

                    iso_dct = isomorphism(xy_gra, prd_gra, stereo=False)
                    if iso_dct:
                        rcts_gra = union_from_sequence(rct_gras)
                        prds_gra = prd_gra
                        brk_atm_i, brk_atm_j = brk_bnd_key
                        b_brk_bnd_key_i = (
                            iso_dct[frm_bnd_pairs[0][0]],
                            iso_dct[frm_bnd_pairs[0][1]])
                        b_brk_bnd_key_j = (
                            iso_dct[frm_bnd_pairs[1][0]],
                            iso_dct[frm_bnd_pairs[1][1]])
                        b_frm_bnd_key = (
                            iso_dct[brk_atm_i],
                            iso_dct[brk_atm_j])
                        forw_tsg = ts.graph(rcts_gra,
                                            frm_bnd_keys=frm_bnd_pairs,
                                            brk_bnd_keys=(brk_bnd_key,))
                        back_tsg = ts.graph(prds_gra,
                                            brk_bnd_keys=[
                                                b_brk_bnd_key_i,
                                                b_brk_bnd_key_j],
                                            frm_bnd_keys=(b_frm_bnd_key,))
                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.Typ.DOUBLE_INSERTION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=list(map(atom_keys, rct_gras)),
                            prds_keys=list(map(atom_keys, prd_gras)),
                        )
                        rxns.append(rxn)

    return rxns


def two_bond_additions(rct_gras, prd_gras):
    """ two bond additions
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 1:
        rct_gras = sort_reagents(rct_gras)
        x_gra, y_gra = rct_gras
        prd_gra, = prd_gras
        x_atm_keys = frozenset().union(
            unsaturated_atom_keys(x_gra),
            lone_pair_atom_keys(x_gra)
        )
        y_atm_keys = frozenset().union(
            unsaturated_atom_keys(y_gra),
            lone_pair_atom_keys(y_gra)
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
            xy_gra = add_bonds(
                union(x_gra, y_gra),
                [set(frm_bnd_keys[0]), set(frm_bnd_keys[1])])

            iso_dct = isomorphism(xy_gra, prd_gra, stereo=False)
            if iso_dct:
                rcts_gra = union_from_sequence(rct_gras)
                prds_gra = prd_gra
                b_brk_bnd_keys = [
                    [iso_dct[frm_bnd_keys[0][0]], iso_dct[frm_bnd_keys[0][1]]],
                    [iso_dct[frm_bnd_keys[1][0]], iso_dct[frm_bnd_keys[1][1]]]
                ]
                forw_tsg = ts.graph(rcts_gra,
                                    frm_bnd_keys=frm_bnd_keys,
                                    brk_bnd_keys=[])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[],
                                    brk_bnd_keys=b_brk_bnd_keys)

                # Create the reaction object
                rxn = from_forward_reverse(
                    cla=ReactionClass.Typ.ADDITION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                )
                rxns.append(rxn)

    return rxns


# 3. Insertions
def insertions(rct_gras, prd_gras):
    """ find insertions consistent with these reactants and products

    :param rct_gras: reactant graphs (must have non-overlapping keys)
    :param prd_gras: product graphs (must have non-overlapping keys)

    Implemented as the reverse of an elimination reaction.

    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = tuple(map(reverse, eliminations(prd_gras, rct_gras)))
    return rxns


# 4. Substitutions
def substitutions(rct_gras, prd_gras):
    """ find substitutions consistent with these reactants and products

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
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = []

    if len(rct_gras) == 2 and len(prd_gras) == 2:
        rct_gra = union_from_sequence(rct_gras)
        prd_gra = union_from_sequence(prd_gras)

        # Loop over both orders of reactants: A+B and B+A
        for rgra1, rgra2 in itertools.permutations(rct_gras):
            bnd_keys = bond_keys(rgra1)
            atom_symb_dct = automol.graph.atom_symbols(rgra1)
            rad_keys = unsaturated_atom_keys(rgra2)

            # Break all possible bonds in total reactant
            for bnd_key, rad_key in itertools.product(bnd_keys, rad_keys):
                gra = remove_bonds(rct_gra, [bnd_key])

                # Form all possible bonds between rad site and non-H atoms
                frm_keys = ()
                for key in bnd_key:
                    frm_symb = atom_symb_dct[key]
                    if frm_symb != 'H':
                        frm_keys += (key,)

                for frm_key in frm_keys:
                    gra_ = add_bonds(gra, [(frm_key, rad_key)])

                    inv_dct = isomorphism(gra_, prd_gra, stereo=False)
                    if inv_dct:
                        brk_key2, = bnd_key - {frm_key}
                        f_frm_bnd_key = (frm_key, rad_key)
                        f_brk_bnd_key = (frm_key, brk_key2)
                        b_frm_bnd_key = (inv_dct[frm_key], inv_dct[brk_key2])
                        b_brk_bnd_key = (inv_dct[frm_key], inv_dct[rad_key])

                        forw_tsg = ts.graph(rct_gra,
                                            frm_bnd_keys=[f_frm_bnd_key],
                                            brk_bnd_keys=[f_brk_bnd_key])
                        back_tsg = ts.graph(prd_gra,
                                            frm_bnd_keys=[b_frm_bnd_key],
                                            brk_bnd_keys=[b_brk_bnd_key])

                        rcts_atm_keys = [atom_keys(rgra1), atom_keys(rgra2)]

                        prds_atm_keys = list(map(atom_keys, prd_gras))
                        if inv_dct[rad_key] not in prds_atm_keys[0]:
                            prds_atm_keys = list(reversed(prds_atm_keys))

                        # Create the reaction object
                        rxn = from_forward_reverse(
                            cla=ReactionClass.Typ.SUBSTITUTION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        )
                        rxns.append(rxn)

    return rxns


def _is_graph(gra):
    """Determine whether the argument is a graph, instead of a sequence
    """
    if hasattr(gra, '__len__') and len(gra) == 2:
        if all(isinstance(o, dict) for o in gra):
            return True

    return False


def find(rct_gras, prd_gras, stereo=False):
    """ find all reactions consistent with these reactants and products

    :param rct_gras: graphs for the reactants without overlapping keys
    :param prd_gras: graphs for the products without overlapping keys
    :param stereo: Find stereo-specified reaction objects?
    :type stereo: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]
    """
    if _is_graph(rct_gras):
        assert _is_graph(prd_gras)
        rct_gras = automol.graph.connected_components(rct_gras)
        prd_gras = automol.graph.connected_components(prd_gras)

    rct_gras0, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras0, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rct_gras = list(map(automol.graph.without_stereo, rct_gras0))
    prd_gras = list(map(automol.graph.without_stereo, prd_gras0))

    assert_is_valid_reagent_graph_list(rct_gras)
    assert_is_valid_reagent_graph_list(prd_gras)

    # check whether this is a valid reaction
    rct_fmls = list(map(formula, rct_gras))
    prd_fmls = list(map(formula, prd_gras))
    rct_strs = list(map(automol.formula.string, rct_fmls))
    prd_strs = list(map(automol.formula.string, prd_fmls))
    assert automol.formula.reac.is_valid_reaction(rct_fmls, prd_fmls), (
        f"Invalid reaction: {str(rct_strs):s} -> {str(prd_strs):s}")

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
        rxns = unique(rxns)
        for rxn in rxns:
            if not stereo:
                all_rxns.append(rxn)
            else:
                srxns = expand_stereo_for_reaction(rxn, rct_gras0, prd_gras0)
                all_rxns.extend(srxns)

    return tuple(all_rxns)


def find_from_chi(rct_chis, prd_chis):
    """ find all reaction classes consistent with these reactants and products

    :param rct_chis: inchis for the reactants
    :param prd_chis: inchis for the products
    :returns: a list of reaction classes
    :rtype: tuple[str]
    """
    rct_geos = list(map(automol.chi.geometry, rct_chis))
    prd_geos = list(map(automol.chi.geometry, prd_chis))
    rct_gras = list(map(automol.geom.graph_without_stereo, rct_geos))
    prd_gras = list(map(automol.geom.graph_without_stereo, prd_geos))
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)
    return find(rct_gras, prd_gras)


# helpers
def _partial_hydrogen_abstraction(qh_gra, q_gra):
    rets = []

    h_atm_key = max(atom_keys(q_gra)) + 1
    uns_atm_keys = unsaturated_atom_keys(q_gra)
    for atm_key in uns_atm_keys:
        q_gra_h = add_atom_explicit_hydrogens(
            q_gra, {atm_key: [h_atm_key]})
        inv_atm_key_dct = isomorphism(q_gra_h, qh_gra, stereo=False)
        if inv_atm_key_dct:
            qh_q_atm_key = inv_atm_key_dct[atm_key]
            qh_h_atm_keys = atom_neighbor_atom_keys(
                qh_gra, qh_q_atm_key, symb='H')
            q_q_atm_key = atm_key
            for qh_h_atm_key in qh_h_atm_keys:
                rets.append((qh_q_atm_key, qh_h_atm_key, q_q_atm_key))

    return rets


# Analyze changes in the spin state to ID the spin crossing
def intersystem_crossing(rxn_muls):
    """ Assess if there is a difference between the reactant and
        product multiplicities to see if there is a change in spin

         NOT CORRECT
    """

    rct_spin_sum, prd_spin_sum = 0, 0
    for rct_mul in rxn_muls[0]:
        rct_spin_sum += (rct_mul - 1.)/2.
    for prd_mul in rxn_muls[1]:
        prd_spin_sum += (prd_mul - 1.)/2.

    return (rct_spin_sum != prd_spin_sum)
