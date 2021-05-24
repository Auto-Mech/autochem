"""
    Enumerate possible reactions
"""

import itertools
from automol import par
from automol.util import dict_
from automol.graph import ts
from automol.graph import relabel
from automol.graph import union
from automol.graph import implicit
from automol.graph import atom_keys
from automol.graph import bond_keys
from automol.graph import remove_atoms
from automol.graph import remove_bonds
from automol.graph import backbone_keys
from automol.graph import add_bonds
from automol.graph import add_bonded_atom
from automol.graph import rings_atom_keys
from automol.graph import connected_components
from automol.graph import unsaturated_atom_keys
from automol.graph import atom_neighbor_atom_key
from automol.graph import atoms_neighbor_atom_keys
from automol.graph import resonance_avg_bond_orders
from automol.reac._reac import Reaction
from automol.reac._reac import ts_unique
from automol.reac._reac import filter_viable_reactions
from automol.reac._util import argsort_reagents
from automol.reac._util import assert_is_valid_reagent_graph_list


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migrations(rct_gras, viable_only=True):
    """ find all possible hydrogen migration products for these reactants

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
        rct_gra, = rct_gras

        # Identify unsaturated sites
        rct_add_key = max(atom_keys(rct_gra)) + 1
        rct_rad_keys = unsaturated_atom_keys(rct_gra)
        rct_hyd_keys = atom_keys(rct_gra, sym='H')

        for rct_rad_key in rct_rad_keys:
            # Add a hydrogen to the radical/unsaturated site
            rct_h_gra = add_bonded_atom(
                rct_gra, 'H', rct_rad_key, bnd_atm_key=rct_add_key)

            # Identify donor sites
            rct_don_keys = backbone_keys(rct_h_gra) - {rct_rad_key}
            for rct_don_key in rct_don_keys:
                rct_hyd_key = atom_neighbor_atom_key(
                    rct_gra, rct_don_key, symbs_first=['H'], symbs_last=[])
                if rct_hyd_key in rct_hyd_keys:
                    prd_gra = remove_atoms(rct_h_gra, {rct_hyd_key})
                    prd_gra = relabel(prd_gra, {rct_add_key: rct_hyd_key})

                    forw_tsg = ts.graph(
                        rct_gra,
                        frm_bnd_keys=[(rct_rad_key, rct_hyd_key)],
                        brk_bnd_keys=[(rct_don_key, rct_hyd_key)])

                    back_tsg = ts.graph(
                        prd_gra,
                        frm_bnd_keys=[(rct_don_key, rct_hyd_key)],
                        brk_bnd_keys=[(rct_rad_key, rct_hyd_key)])

                    rxns.append(Reaction(
                        rxn_cls=par.ReactionClass.HYDROGEN_MIGRATION,
                        forw_tsg=forw_tsg,
                        back_tsg=back_tsg,
                        rcts_keys=[atom_keys(rct_gra)],
                        prds_keys=[atom_keys(prd_gra)],
                    ))

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return ts_unique(rxns)


# 2. Beta scissions
def beta_scissions(rct_gras, viable_only=True):
    """ find all possible beta scission products for these reactants

    :param rct_gras: graphs for the reactants, without stereo and without
        overlapping keys
    :param viable_only: Filter out reactions with non-viable products?
    :type viable_only: bool
    :returns: a list of Reaction objects
    :rtype: tuple[Reaction]

    Beta scissions are enumerated by identifying all pure single bonds (single
    bonds with no resonances), and looping over the results of breaking each of
    them. If this gives rise to two distinct fragments, the reaction is added
    to the list.
    """
    assert_is_valid_reagent_graph_list(rct_gras)

    rxns = []

    if len(rct_gras) == 1:
        rct_gra, = rct_gras

        # Identify all pure single bonds
        avg_bnd_ord_dct = resonance_avg_bond_orders(implicit(rct_gra))
        brk_bnd_keys = dict_.keys_by_value(avg_bnd_ord_dct, lambda x: x == 1)

        for brk_bnd_key in brk_bnd_keys:
            prds_gra = remove_bonds(rct_gra, [brk_bnd_key])
            prd_gras = connected_components(prds_gra)

            if len(prd_gras) == 2:
                prd_gras = list(map(prd_gras.__getitem__,
                                    argsort_reagents(prd_gras)))

                forw_tsg = ts.graph(rct_gra,
                                    frm_bnd_keys=[],
                                    brk_bnd_keys=[brk_bnd_key])
                back_tsg = ts.graph(prds_gra,
                                    frm_bnd_keys=[brk_bnd_key],
                                    brk_bnd_keys=[])

                # Create the reaction object
                rxns.append(Reaction(
                    rxn_cls=par.ReactionClass.BETA_SCISSION,
                    forw_tsg=forw_tsg,
                    back_tsg=back_tsg,
                    rcts_keys=list(map(atom_keys, rct_gras)),
                    prds_keys=list(map(atom_keys, prd_gras)),
                ))

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return ts_unique(rxns)


# 3. Ring-forming scissions (skip for now)


# 4. Eliminations
def eliminations(rct_gras, viable_only=True):
    """ find all possible elimination products for these reactants

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
        rct_gra, = rct_gras

        ngb_keys_dct = atoms_neighbor_atom_keys(rct_gra)

        frm1_keys = atom_keys(rct_gra, excl_syms=('H',))
        frm2_keys = atom_keys(rct_gra)
        bnd_keys = bond_keys(rct_gra)

        frm_bnd_keys = [(frm1_key, frm2_key) for frm1_key, frm2_key
                        in itertools.product(frm1_keys, frm2_keys)
                        if frm1_key != frm2_key and
                        not frozenset({frm1_key, frm2_key}) in bnd_keys]

        for frm1_key, frm2_key in frm_bnd_keys:
            # Bond the radical atom to the hydrogen atom
            prds_gra = add_bonds(rct_gra, [(frm2_key, frm1_key)])

            # Get keys to the ring formed by this extra bond
            rng_keys = next((ks for ks in rings_atom_keys(prds_gra)
                             if frm2_key in ks and frm1_key in ks), None)
            # Eliminations (as far as I can tell) only happen through TSs with
            # 3- or 4-membered rings
            if rng_keys is not None and len(rng_keys) < 5:
                frm1_ngb_key, = ngb_keys_dct[frm1_key] & set(rng_keys)
                frm2_ngb_key, = ngb_keys_dct[frm2_key] & set(rng_keys)

                # Break the bonds on either side of the newly formed bond
                prds_gra = remove_bonds(prds_gra, [(frm1_key, frm1_ngb_key)])
                prds_gra = remove_bonds(prds_gra, [(frm2_key, frm2_ngb_key)])

                prd_gras = connected_components(prds_gra)

                if len(prd_gras) == 2:
                    forw_tsg = ts.graph(
                        rct_gra,
                        frm_bnd_keys=[(frm1_key, frm2_key)],
                        brk_bnd_keys=[(frm1_key, frm1_ngb_key),
                                      (frm2_key, frm2_ngb_key)])
                    back_tsg = ts.graph(
                        prds_gra,
                        frm_bnd_keys=[(frm1_key, frm1_ngb_key),
                                      (frm2_key, frm2_ngb_key)],
                        brk_bnd_keys=[(frm1_key, frm2_key)])

                    rcts_atm_keys = list(map(atom_keys, rct_gras))
                    prds_atm_keys = list(map(atom_keys, prd_gras))

                    if frm2_key not in prds_atm_keys[1]:
                        prds_atm_keys = list(reversed(prds_atm_keys))

                    # Create the reaction object
                    rxns.append(Reaction(
                        rxn_cls=par.ReactionClass.ELIMINATION,
                        forw_tsg=forw_tsg,
                        back_tsg=back_tsg,
                        rcts_keys=rcts_atm_keys,
                        prds_keys=prds_atm_keys,
                    ))

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return ts_unique(rxns)


# Bimolecular reactions
# 2. Additions
def additions(rct_gras, viable_only=True):
    """ find all possible addition products for these reactants

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
        rct1_gra, rct2_gra, = rct_gras

        rct1_atm_keys = unsaturated_atom_keys(rct1_gra)
        rct2_atm_keys = unsaturated_atom_keys(rct2_gra)

        for frm_bnd_key in itertools.product(rct1_atm_keys, rct2_atm_keys):
            rcts_gra = union(rct1_gra, rct2_gra)
            prd_gra = add_bonds(rcts_gra, [frm_bnd_key])
            prd_gras = [prd_gra]

            forw_tsg = ts.graph(rcts_gra,
                                frm_bnd_keys=[frm_bnd_key],
                                brk_bnd_keys=[])
            back_tsg = ts.graph(prd_gra,
                                frm_bnd_keys=[],
                                brk_bnd_keys=[frm_bnd_key])

            # Create the reaction object
            rxns.append(Reaction(
                rxn_cls=par.ReactionClass.ADDITION,
                forw_tsg=forw_tsg,
                back_tsg=back_tsg,
                rcts_keys=list(map(atom_keys, rct_gras)),
                prds_keys=list(map(atom_keys, prd_gras)),
            ))

    if viable_only:
        rxns = filter_viable_reactions(rxns)

    return ts_unique(rxns)


def enumerate_reactions(rct_gras, viable_only=True):
    """ enumerate all possible reactions that a given set of reactants might
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

    # Cycle through the different finders and gather all possible reactions
    finders_ = [
        # unimolecular reactions
        hydrogen_migrations,
        beta_scissions,
        # ring_forming_scissions,
        eliminations,
        # bimolecular reactions
        # hydrogen_abstractions,
        additions,
        # insertions,
        # substitutions,
    ]

    rxns = tuple(itertools.chain(*(f_(rct_gras, viable_only=viable_only)
                                   for f_ in finders_)))

    return rxns


if __name__ == "__main__":
    import automol

    # RCT_SMIS = ['C=CCC[CH2]']         # hydrogen migration
    # RCT_SMIS = ['C=C[CH]CC']          # beta scission
    # RCT_SMIS = ['C=CC=C', '[CH3]']    # addition
    RCT_SMIS = ['CCCO[O]']            # elimination
    RCT_ICHS = list(map(automol.smiles.inchi, RCT_SMIS))
    RCT_GEOS = list(map(automol.convert.geom.geometry, RCT_ICHS))
    RCT_GRAS = list(map(automol.convert.geom.connectivity_graph, RCT_GEOS))
    RCT_GRAS, _ = automol.graph.standard_keys_for_sequence(RCT_GRAS)

    # RXNS = enumerate_reactions(RCT_GRAS)
    RXNS = eliminations(RCT_GRAS)
    print(len(RXNS))

    for rxn in RXNS:
        print(automol.graph.string(rxn.forward_ts_graph))
        print()

    print(len(RXNS))
