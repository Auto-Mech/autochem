"""
    Enumerate possible reactions
"""

import itertools

from ..const import ReactionClass, ReactionInfo, ReactionSpin
from ..graph import (
    add_bonded_atom,
    add_bonds,
    addition_atom_keys,
    are_equivalent_atoms,
    atom_equivalence_class_reps,
    atom_keys,
    atom_neighbor_atom_key,
    atom_symbols,
    atoms_neighbor_atom_keys,
    backbone_keys,
    beta_scission_bond_keys,
    bond_equivalence_class_reps,
    bond_keys,
    connected_components,
    hydroperoxy_groups,
    kekules_bond_orders_averaged,
    radical_atom_keys,
    relabel,
    remove_atoms,
    remove_bonds,
    rings_atom_keys,
    sort_by_size,
    ts,
    union,
    unsaturated_atom_keys,
)
from ..util import dict_
from ._0core import filter_viable_reactions, from_forward_reverse, unique
from ._1util import assert_is_valid_reagent_graph_list
from ._instab import instability_product_graphs


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migrations(rct_gras, viable_only=True):
    """find all possible hydrogen migration reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Hydrogen migrations are enumerated looping over unsaturated sites, adding
    hydrogens to them, and looping over non-equivalent heavy atoms and removing
    hydrgens from them.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        (rct_gra,) = rct_gras

        # Identify unsaturated sites
        rct_add_key = max(atom_keys(rct_gra)) + 1
        rct_rad_keys = radical_atom_keys(rct_gra)
        rct_hyd_keys = atom_keys(rct_gra, symb="H")

        for rct_rad_key in rct_rad_keys:
            # Add a hydrogen to the radical/unsaturated site
            rct_h_gra = add_bonded_atom(
                rct_gra, "H", rct_rad_key, bnd_atm_key=rct_add_key
            )

            # Identify donor sites
            rct_don_keys = backbone_keys(rct_h_gra) - rct_rad_keys
            for rct_don_key in rct_don_keys:
                rct_hyd_key = atom_neighbor_atom_key(
                    rct_gra, rct_don_key, symbs_first=["H"], symbs_last=[]
                )
                if rct_hyd_key in rct_hyd_keys:
                    prd_gra = remove_atoms(rct_h_gra, {rct_hyd_key}, stereo=True)
                    prd_gra = relabel(prd_gra, {rct_add_key: rct_hyd_key})

                    forw_tsg = ts.graph(
                        rct_gra,
                        frm_bnd_keys=[(rct_rad_key, rct_hyd_key)],
                        brk_bnd_keys=[(rct_don_key, rct_hyd_key)],
                    )

                    back_tsg = ts.graph(
                        prd_gra,
                        frm_bnd_keys=[(rct_don_key, rct_hyd_key)],
                        brk_bnd_keys=[(rct_rad_key, rct_hyd_key)],
                    )

                    rxns.append(
                        from_forward_reverse(
                            cla=ReactionClass.HYDROGEN_MIGRATION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=[atom_keys(rct_gra)],
                            prds_keys=[atom_keys(prd_gra)],
                        )
                    )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 2. Homolytic scissions
# AVC comment: replaced commented-out if statement by changing the default here
# to False. Is the viability check broken for this case?
# def homolytic_scissions(rct_gras, viable_only=True):
def homolytic_scissions(rct_gras, viable_only=False):
    """find all possible homolytic scission reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Homolytic scissions are enumerated by identifying all pure single bonds
    (single bonds with no resonances), and looping over the results of
    breaking each of them. If this gives rise to two distinct
    fragments, the reaction is added to the list.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        (rct_gra,) = rct_gras

        # Identify all pure single bonds involving radical site neighbor
        avg_bnd_ord_dct = kekules_bond_orders_averaged(rct_gra)
        brk_bnd_keys = dict_.keys_by_value(avg_bnd_ord_dct, lambda x: x == 1)

        for brk_bnd_key in brk_bnd_keys:
            prds_gra = remove_bonds(rct_gra, [brk_bnd_key])
            prd_gras = connected_components(prds_gra)

            if len(prd_gras) == 2:
                prd_gras = sort_by_size(prd_gras)

                forw_tsg = ts.graph(
                    rct_gra, frm_bnd_keys=[], brk_bnd_keys=[brk_bnd_key]
                )
                back_tsg = ts.graph(
                    prds_gra, frm_bnd_keys=[brk_bnd_key], brk_bnd_keys=[]
                )

                # Create the reaction object
                rxns.append(
                    from_forward_reverse(
                        cla=ReactionClass.HOMOLYT_SCISSION,
                        ftsg=forw_tsg,
                        rtsg=back_tsg,
                        rcts_keys=list(map(atom_keys, rct_gras)),
                        prds_keys=list(map(atom_keys, prd_gras)),
                    )
                )

    # Dummy line to fix linting checks
    assert viable_only or not viable_only
    # filter removes all reactions
    # if viable_only:
    #    rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 3. Beta scissions
def beta_scissions(rct_gras, viable_only=True):
    """find all possible beta scission reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    FIX DESCRIPTION:
    Beta scissions are enumerated by identifying all pure single bonds (single
    bonds with no resonances), and looping over the results of breaking each of
    them. If this gives rise to two distinct fragments, the reaction is added
    to the list.

    AVC CHANGE:
    Don't require two distinct fragments -- allow for ring opening.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        (rct_gra,) = rct_gras

        beta_bnd_keys = beta_scission_bond_keys(rct_gra)
        for brk_bnd_key in beta_bnd_keys:
            prds_gra = remove_bonds(rct_gra, [brk_bnd_key])
            prd_gras = connected_components(prds_gra)

            prd_gras = sort_by_size(prd_gras)

            forw_tsg = ts.graph(rct_gra, frm_bnd_keys=[], brk_bnd_keys=[brk_bnd_key])
            back_tsg = ts.graph(prds_gra, frm_bnd_keys=[brk_bnd_key], brk_bnd_keys=[])

            # Create the reaction object
            rxns.append(
                from_forward_reverse(
                    cla=ReactionClass.BETA_SCISSION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                )
            )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 4. Ring-forming scissions (skip for now)
def ring_forming_scissions(rct_gras, viable_only=True):
    """find all possible ring-forming scission reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Right now it takes the lazy, chemically specific approach of finding
    C-O-O-H groups and forming a bond between the O of the C-O bond
    and radical sites of the species, while breaking the O-O bond.
    """

    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        (rct_gra,) = rct_gras

        # Identify the radical sites and COOH groups
        rad_keys = radical_atom_keys(rct_gra)
        cooh_grps = hydroperoxy_groups(rct_gra)

        # Get the bnd keys for filtering
        bnd_keys = bond_keys(rct_gra)

        # Set the forming and breaking bonds by looping over COOH groups
        rxn_bnd_keys = ()
        for cooh_grp in cooh_grps:
            brk_bnd_key = frozenset(cooh_grp[1:3])
            for rad_key in rad_keys:
                frm_bnd_key = frozenset({rad_key, cooh_grp[1]})
                # Only includ frm bnd if it does not exist
                # e.g., CC[C]OO already has frm bnd -> no rxn possible
                if frm_bnd_key not in bnd_keys:
                    rxn_bnd_keys += ((frm_bnd_key, brk_bnd_key),)

        # Form reactions with all combinations of frm and brk bnds
        for frm_bnd_key, brk_bnd_key in rxn_bnd_keys:
            prds_gra = rct_gra
            prds_gra = add_bonds(prds_gra, [frm_bnd_key])
            prds_gra = remove_bonds(prds_gra, [brk_bnd_key])
            prd_gras = connected_components(prds_gra)

            if len(prd_gras) == 2:
                prd_gras = sort_by_size(prd_gras)

                forw_tsg = ts.graph(
                    rct_gra, frm_bnd_keys=[frm_bnd_key], brk_bnd_keys=[brk_bnd_key]
                )
                back_tsg = ts.graph(
                    prds_gra, frm_bnd_keys=[brk_bnd_key], brk_bnd_keys=[frm_bnd_key]
                )
                # Create the reaction object
                rxns.append(
                    from_forward_reverse(
                        cla=ReactionClass.RING_FORM_SCISSION,
                        ftsg=forw_tsg,
                        rtsg=back_tsg,
                        rcts_keys=list(map(atom_keys, rct_gras)),
                        prds_keys=list(map(atom_keys, prd_gras)),
                    )
                )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 5. Eliminations
def eliminations(rct_gras, viable_only=True):
    """find all possible elimination reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Eliminations are enumerated by forming a bond between an attacking heavy
    atom and another atom not initially bonded to it, forming a ring. The bond
    adjacent to the attacked atom is then broken, along with a second bond in
    the ring, downstream from the attacking heavy atom, away from the attacked
    atom.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        (rct_gra,) = rct_gras

        ngb_keys_dct = atoms_neighbor_atom_keys(rct_gra)

        # frm1_keys = atom_keys(rct_gra, excl_syms=('H',))
        frm1_keys = unsaturated_atom_keys(rct_gra)
        frm2_keys = atom_keys(rct_gra)
        rct_symbs = atom_symbols(rct_gra)
        frm1_keys_o = frozenset(key for key in frm1_keys if rct_symbs[key] == "O")
        frm2_keys_h = frozenset(key for key in frm2_keys if rct_symbs[key] == "H")
        bnd_keys = bond_keys(rct_gra)

        frm_bnd_keys = [
            (frm1_key, frm2_key)
            for frm1_key, frm2_key in itertools.product(frm1_keys_o, frm2_keys_h)
            if frm1_key != frm2_key and not frozenset({frm1_key, frm2_key}) in bnd_keys
        ]

        for frm1_key, frm2_key in frm_bnd_keys:
            # Bond the radical atom to the hydrogen atom
            prds_gra = add_bonds(rct_gra, [(frm2_key, frm1_key)])

            # Get keys to the ring formed by this extra bond
            rng_keys = next(
                (
                    ks
                    for ks in rings_atom_keys(prds_gra)
                    if frm2_key in ks and frm1_key in ks
                ),
                None,
            )

            # HO2 Eliminations only happen through TSs with
            # 5-membered rings
            if rng_keys is not None and len(rng_keys) == 5:
                # Get breaking bondsRemove bond
                # Endpoints of ring
                brk_keys1 = (rng_keys[0], rng_keys[-1])
                # two ring
                brk_keys2 = ()
                for neigh1 in ngb_keys_dct[rng_keys[0]]:
                    for neigh2 in ngb_keys_dct[neigh1]:
                        if neigh2 in rng_keys and neigh2 not in brk_keys1:
                            brk_keys2 = (neigh1, neigh2)
                            break

                # Remove bond
                prds_gra = remove_bonds(prds_gra, [brk_keys1, brk_keys2])
                prd_gras = connected_components(prds_gra)

                if len(prd_gras) == 2:
                    forw_tsg = ts.graph(
                        rct_gra,
                        frm_bnd_keys=[(frm1_key, frm2_key)],
                        brk_bnd_keys=[brk_keys1, brk_keys2],
                    )
                    back_tsg = ts.graph(
                        prds_gra,
                        frm_bnd_keys=[brk_keys1, brk_keys2],
                        brk_bnd_keys=[(frm1_key, frm2_key)],
                    )

                    rcts_atm_keys = list(map(atom_keys, rct_gras))
                    prds_atm_keys = list(map(atom_keys, prd_gras))

                    if frm2_key not in prds_atm_keys[1]:
                        prds_atm_keys = list(reversed(prds_atm_keys))

                    # Create the reaction object
                    rxns.append(
                        from_forward_reverse(
                            cla=ReactionClass.ELIMINATION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        )
                    )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstractions(rct_gras, viable_only=True):
    """find hydrogen abstraction reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Hydrogen abstractions are enumerated by looping over unique unsaturated
    atoms on one molecule and abstracting from unique atoms on the other.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 2:
        for q1h_gra, q2_gra in itertools.permutations(rct_gras):
            hyd_keys = atom_keys(q1h_gra, symb="H")

            # Identify unique heavy atoms as potential donors
            don_keys = atom_keys(q1h_gra, excl_symbs=("H",))
            don_keys = atom_equivalence_class_reps(q1h_gra, don_keys)

            # Identify unique unsaturated atoms as potential attackers
            att_keys = unsaturated_atom_keys(q2_gra)
            att_keys = atom_equivalence_class_reps(q2_gra, att_keys)

            for don_key, att_key in itertools.product(don_keys, att_keys):
                hyd_key = atom_neighbor_atom_key(
                    q1h_gra, don_key, symbs_first=["H"], symbs_last=[]
                )
                if hyd_key in hyd_keys:
                    # Remove a hydrogen from the donor site
                    q1_gra = remove_atoms(q1h_gra, {hyd_key}, stereo=True)
                    # Add a hydrogen atom to the attacker site
                    q2h_gra = add_bonded_atom(q2_gra, "H", att_key, bnd_atm_key=hyd_key)

                    rcts_gra = union(q1h_gra, q2_gra)
                    prds_gra = union(q2h_gra, q1_gra)

                    forw_tsg = ts.graph(
                        rcts_gra,
                        frm_bnd_keys=[(att_key, hyd_key)],
                        brk_bnd_keys=[(don_key, hyd_key)],
                    )

                    back_tsg = ts.graph(
                        prds_gra,
                        frm_bnd_keys=[(don_key, hyd_key)],
                        brk_bnd_keys=[(att_key, hyd_key)],
                    )

                    rcts_atm_keys = list(map(atom_keys, [q1h_gra, q2_gra]))
                    prds_atm_keys = list(map(atom_keys, [q2h_gra, q1_gra]))

                    # Create the reaction object
                    rxns.append(
                        from_forward_reverse(
                            cla=ReactionClass.HYDROGEN_ABSTRACTION,
                            ftsg=forw_tsg,
                            rtsg=back_tsg,
                            rcts_keys=rcts_atm_keys,
                            prds_keys=prds_atm_keys,
                        )
                    )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 2. Additions
def additions(rct_gras, viable_only=True):
    """find all possible addition reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Additions are enumerated by joining an unsaturated site on one reactant to
    an unsaturated site on the other.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 2:
        rct_gras = sort_by_size(rct_gras)
        rct1_gra, rct2_gra = rct_gras

        rct1_atm_keys = addition_atom_keys(rct1_gra)
        rct2_atm_keys = addition_atom_keys(rct2_gra)

        for frm_bnd_key in itertools.product(rct1_atm_keys, rct2_atm_keys):
            rcts_gra = union(rct1_gra, rct2_gra)
            prd_gra = add_bonds(rcts_gra, [frm_bnd_key])
            prd_gras = [prd_gra]

            forw_tsg = ts.graph(rcts_gra, frm_bnd_keys=[frm_bnd_key], brk_bnd_keys=[])
            back_tsg = ts.graph(prd_gra, frm_bnd_keys=[], brk_bnd_keys=[frm_bnd_key])

            # Create the reaction object
            rxns.append(
                from_forward_reverse(
                    cla=ReactionClass.ADDITION,
                    ftsg=forw_tsg,
                    rtsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                )
            )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# 3. Insertions
def insertions(rct_gras, viable_only=True):
    """find all possible insertion reactions for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Insertions are enumerated by looping over carbenes and multiple bonds on
    one reactant, which serve as a source of potential "attacking" atoms for
    the insertion, and looping over single bonds that could be inserted into on
    the other reactant. For lack of a better term, we can call these "donating
    atoms". The insertion then looks as follows:

        A1=A2         A1
        .  .    or    .
        .  .         .  .
        D1-D2        D1-D2

    where two bonds are formed between the A and D atoms and the bond between
    the two D atoms is broken.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 2:
        for rct1_gra, rct2_gra in itertools.permutations(rct_gras):
            rcts_gra = union(rct1_gra, rct2_gra)

            # Carbenes on R1 are potential attacking atoms
            atm_keys = radical_atom_keys(rct1_gra, min_valence=2.0)
            atm_keys = tuple(atm_keys)
            # So are atoms on either side of a multiple bond
            bnd_keys = dict_.keys_by_value(
                kekules_bond_orders_averaged(rct1_gra), lambda x: x > 1.0
            )
            bnd_keys = bond_equivalence_class_reps(rct1_gra, bnd_keys)
            # Use this to form a list of attacking atom pairs for R1
            att_pairs = list(map(tuple, map(sorted, bnd_keys)))
            att_pairs += list(zip(atm_keys, atm_keys))

            # As donor pairs, consider single bonds on R2
            don_bnd_keys = dict_.keys_by_value(
                kekules_bond_orders_averaged(rct2_gra), lambda x: x == 1.0
            )
            don_bnd_keys = bond_equivalence_class_reps(rct2_gra, don_bnd_keys)
            don_pairs = list(map(tuple, map(sorted, don_bnd_keys)))

            for att_pair, don_pair in itertools.product(att_pairs, don_pairs):
                if not (
                    are_equivalent_atoms(rct1_gra, *att_pair)
                    or are_equivalent_atoms(rct2_gra, *don_pair)
                ):
                    don_pairs_ = list(itertools.permutations(don_pair))
                else:
                    don_pairs_ = [don_pair]

                for don_pair_ in don_pairs_:
                    att1_key, att2_key = att_pair
                    don1_key, don2_key = don_pair_

                    prds_gra = rcts_gra
                    prds_gra = add_bonds(
                        prds_gra, [(att1_key, don1_key), (att2_key, don2_key)]
                    )
                    prds_gra = remove_bonds(prds_gra, [(don1_key, don2_key)])

                    prd_gras = connected_components(prds_gra)

                    if len(prd_gras) == 1:
                        forw_tsg = ts.graph(
                            rcts_gra,
                            frm_bnd_keys=[(att1_key, don1_key), (att2_key, don2_key)],
                            brk_bnd_keys=[(don1_key, don2_key)],
                        )
                        back_tsg = ts.graph(
                            prds_gra,
                            frm_bnd_keys=[(don1_key, don2_key)],
                            brk_bnd_keys=[(att1_key, don1_key), (att2_key, don2_key)],
                        )

                        # Create the reaction object
                        rcts_keys = list(map(atom_keys, [rct1_gra, rct2_gra]))
                        prds_keys = list(map(atom_keys, prd_gras))
                        rxns.append(
                            from_forward_reverse(
                                cla=ReactionClass.INSERTION,
                                ftsg=forw_tsg,
                                rtsg=back_tsg,
                                rcts_keys=rcts_keys,
                                prds_keys=prds_keys,
                            )
                        )

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return unique(rxns)


# Cycle through the different finders and gather all possible reactions
FINDERS = {
    # unimolecular reactions
    ReactionClass.HYDROGEN_MIGRATION: hydrogen_migrations,
    ReactionClass.HOMOLYT_SCISSION: homolytic_scissions,
    ReactionClass.BETA_SCISSION: beta_scissions,
    ReactionClass.RING_FORM_SCISSION: ring_forming_scissions,
    ReactionClass.ELIMINATION: eliminations,
    # bimolecular reactions
    ReactionClass.HYDROGEN_ABSTRACTION: hydrogen_abstractions,
    ReactionClass.ADDITION: additions,
    ReactionClass.INSERTION: insertions,
    # par.ReactionClass.Typ.SUBSTITUTION: substitutions,
}


def enumerate_reactions(rct_gras, rxn_type=None, viable_only=True):
    """enumerate all possible reactions that a given set of reactants might
        undergo

    For single reactants, all unimolecular reactions will be enumerated.
    For pairs of reactants, all bimolecular reactions between them will be
    enumerated.

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]
    """

    # Check if the reactants are unstable and should yield no reactants
    unstable = any(bool(instability_product_graphs(gra)) for gra in rct_gras)

    # Determine products of the reactions if those are possible
    if not unstable:
        if rxn_type is not None:
            rxns = FINDERS[rxn_type](rct_gras, viable_only=viable_only)
        else:
            # Cycle through the all finders and gather all possible reactions
            finders_ = tuple(FINDERS.values())
            rxns = tuple(
                itertools.chain(
                    *(f_(rct_gras, viable_only=viable_only) for f_ in finders_)
                )
            )
    else:
        rxns = ()

    return rxns


# initialize ReactionInfo from a string
def reaction_info_from_string(rxn_str):
    """
    :param rxn_str: string representation of reaction information
    :type rxn_str: string
    """
    class_ = None
    spin_ = None
    is_rad_rad_ = False
    is_isc_ = False

    if any(
            value in rxn_str
            for value in ReactionSpin):
        spin_ = [value for value in ReactionSpin if value in rxn_str][0]
    else:
        spin_ = ReactionSpin.NONE

    if any(
            value in rxn_str
            for value in ReactionClass):
        class_ = ReactionClass(
            [value for value in ReactionClass if value in rxn_str][0])
    else:
        print(f"Reaction {rxn_str} is not an option")

    is_rad_rad_ = 'radical-radical' in rxn_str
    is_isc_ = 'intersystem-crossing' in rxn_str
    return ReactionInfo(
        class_, spin_, is_rad_rad_, is_isc_)
