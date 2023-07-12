""" stereo functionality for reaction objects
"""
import automol.graph
import automol.geom
import automol.chi
from automol.reac._core import Reaction
from automol.reac._core import reverse
from automol.reac._core import atom_mapping
from automol.reac._core import reactant_graphs
from automol.reac._core import product_graphs
from automol.reac._core import reactants_graph
from automol.reac._core import products_graph
from automol.reac._core import forming_bond_keys
from automol.reac._core import breaking_bond_keys
from automol.reac._core import relabel
from automol.reac._conv import amchi
from automol.graph import ts


def add_stereo_from_geometries(rxn, rct_geos, prd_geos):
    """ Add stereo assignments to this reaction object from geometries.

    :param rxn: a reaction object
    :type rxn: Reaction
    :param rct_geos: the reactant geometries
    :param prd_geos: the product geometries
    :returns: a reaction object with stereo assignments
    :rtype: Reaction
    """
    rct_gras = list(reactant_graphs(rxn))
    prd_gras = list(product_graphs(rxn))

    def _index_dict(gra):
        return dict(map(reversed,
                        enumerate(sorted(automol.graph.atom_keys(gra)))))

    rct_gras = [
        automol.graph.set_stereo_from_geometry(
            rct_gra, rct_geo, geo_idx_dct=_index_dict(rct_gra))
        for rct_gra, rct_geo in zip(rct_gras, rct_geos)]

    prd_gras = [
        automol.graph.set_stereo_from_geometry(
            prd_gra, prd_geo, geo_idx_dct=_index_dict(prd_gra))
        for prd_gra, prd_geo in zip(prd_gras, prd_geos)]

    rcts_gra = automol.graph.union_from_sequence(rct_gras)
    prds_gra = automol.graph.union_from_sequence(prd_gras)

    rev_rxn = reverse(rxn)

    forw_tsg = automol.graph.ts.graph(
        rcts_gra, forming_bond_keys(rxn), breaking_bond_keys(rxn))
    back_tsg = automol.graph.ts.graph(
        prds_gra, forming_bond_keys(rev_rxn), breaking_bond_keys(rev_rxn))

    srxn = Reaction(rxn.class_, forw_tsg, back_tsg,
                    rxn.reactants_keys, rxn.products_keys)
    return srxn


def add_stereo_from_inchis(rxn, rct_ichs, prd_ichs, check=True):
    """ Add stereo assignments to this reaction object from inchis.

    :param rxn: a reaction object
    :type rxn: Reaction
    :param rct_ichs: the reactant inchis
    :param prd_ichs: the product inchis
    :param check: check the inchis for stereo completeness? if false, then
        arbitrary stereo will be assigned for inchis missing stereo assignments
    :param check: bool
    :returns: a reaction object with stereo assignments
    :rtype: Reaction
    """
    if check:
        assert all(map(automol.chi.is_complete, rct_ichs)), (
            f"Some inchis are not complete: {rct_ichs}")
        assert all(map(automol.chi.is_complete, prd_ichs)), (
            f"Some inchis are not complete: {prd_ichs}")

    rct_geos = list(map(automol.chi.geometry, rct_ichs))
    prd_geos = list(map(automol.chi.geometry, prd_ichs))
    srxn, _ = add_stereo_from_unordered_geometries(rxn, rct_geos, prd_geos,
                                                   reorder_atoms=False)
    return srxn


def add_stereo_from_unordered_geometries(rxn, rct_geos, prd_geos,
                                         reorder_atoms=True):
    """ Add stereo assignments to this reaction object from geometries that may
    not match the ordering of reactants, products, and atoms in the reaction
    object.

    :param rxn: a reaction object
    :type rxn: Reaction
    :param rct_geos: the reactant geometries
    :param prd_geos: the product geometries
    :param reorder_atoms: reorder the atoms in the reaction object to match
        those in the reactants or products? if false, stereo will be added to
        `rxn` without changing anything else about it.
    :type reorder_atoms: bool
    :returns: a reaction object with stereo assignments
    :rtype: Reaction
    """
    rct_gras = list(map(automol.geom.graph, rct_geos))
    prd_gras = list(map(automol.geom.graph, prd_geos))
    print('add streo')
    for rgra in rct_gras + prd_gras:
        print(automol.graph.string(rgra))
    found_srxn = None
    order = None

    for srxn in expand_stereo(rxn):
        comp_rct_gras = reactant_graphs(srxn)
        comp_prd_gras = product_graphs(srxn)
        print(comp_rct_gras)
        print(comp_prd_gras)
        rct_order, _ = automol.graph.sequence_isomorphism(
            rct_gras, comp_rct_gras)
        prd_order, _ = automol.graph.sequence_isomorphism(
            prd_gras, comp_prd_gras)
        print('orders', rct_order, prd_order)
        if rct_order is not None and prd_order is not None:
            rct_gras = [rct_gras[i] for i in rct_order]
            prd_gras = [prd_gras[i] for i in prd_order]
            rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
            prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

            comp_rcts_gra = reactants_graph(srxn)
            comp_prds_gra = products_graph(srxn)
            print('again')
            print(comp_rcts_gra)
            print(comp_prds_gra)
            rcts_gra = automol.graph.union_from_sequence(rct_gras)
            prds_gra = automol.graph.union_from_sequence(prd_gras)
            print('again2')
            print(rcts_gra)
            print(prds_gra)
            rct_iso_dct = automol.graph.isomorphism(comp_rcts_gra, rcts_gra)
            prd_iso_dct = automol.graph.isomorphism(comp_prds_gra, prds_gra)

            # if requested, reorder atoms in the reaction object to match those
            # in the geometries
            if reorder_atoms:
                srxn = relabel(srxn, rct_iso_dct)
                srxn = relabel(srxn, prd_iso_dct, product=True)

            found_srxn = srxn
            order = (rct_order, prd_order)
            break

    return found_srxn, order


def expand_stereo(rxn, enant=True):
    """ Expand all possible stereo assignments for the reactants and products
    of this reaction. Only includes possibilities that are mutually consistent
    with each other.

        :param rxn: a reaction object
        :type rxn: Reaction
        :param enant: Include all enantiomers? Otherwise, includes only
            canonical enantiomer species and reactions.
        :type enant: bool
        :returns: a sequence reaction objects with stereo assignments
        :rtype: Reaction
    """
    rxn_cls = rxn.class_
    forw_tsg = automol.graph.without_stereo(rxn.forward_ts_graph)
    back_tsg = rxn.backward_ts_graph
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys

    key_dct = atom_mapping(rxn)

    srxns = []
    for forw_ste_tsg, back_ste_tsg in (
            ts.expand_reaction_stereo(forw_tsg, enant=enant)):
        back_ste_tsg = automol.graph.relabel(back_ste_tsg, key_dct)

        # But for dummy atoms, we could just do the conversion directly,
        # but this avoids loss of dummy atoms from the products graph.
        tsg_ = back_tsg
        tsg_ = automol.graph.set_atom_stereo_parities(
            tsg_, automol.graph.atom_stereo_parities(back_ste_tsg))
        tsg_ = automol.graph.set_bond_stereo_parities(
            tsg_, automol.graph.bond_stereo_parities(back_ste_tsg))
        back_ste_tsg = tsg_

        srxn = Reaction(rxn_cls, forw_ste_tsg, back_ste_tsg,
                        rcts_keys, prds_keys)
        srxns.append(srxn)

    srxns = tuple(srxns)
    return srxns


def stereo_is_physical(srxn):
    """ Does this reaction have consistent stereo assignments for reactants and
    products?

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: True if stereo assignments are consisent, False if not
    :rtype: bool
    """
    ftsg_loc = automol.graph.to_local_stereo(srxn.forward_ts_graph)
    rtsg_loc = automol.graph.to_local_stereo(srxn.backward_ts_graph)
    key_dct = atom_mapping(srxn, rev=True)
    assert key_dct, (f"Forward and backward TS graphs don't match:\n"
                     f"{automol.graph.string(srxn.forward_ts_graph)}\n"
                     f"{automol.graph.string(srxn.backward_ts_graph)}\n")
    rtsg_loc = automol.graph.relabel(rtsg_loc, key_dct)
    return automol.graph.ts.reaction_stereo_is_physical(ftsg_loc, rtsg_loc)


def is_canonical_enantiomer(srxn):
    """ Does this reaction have a canonical combination of enantiomers?

        :param srxn: a reaction object with stereo assignments
        :type srxn: Reaction
        :returns: Whether or not the reaction is canonical
        :rtype: bool
    """
    rct_chis, prd_chis = amchi(srxn)
    return automol.amchi.is_canonical_enantiomer_reaction(rct_chis, prd_chis)


def canonical_enantiomer(srxn):
    """ Convert this reaction into a canonical combination of enantiomers

        :param srxn: a reaction object with stereo assignments
        :type srxn: Reaction
        :returns: Whether or not the reaction is canonical
        :rtype: bool
    """
    if not is_canonical_enantiomer(srxn):
        srxn = reflect(srxn)
    return srxn


def reflect(srxn):
    """ Reflect all graphs in this reaction, to obtain their mirror images

        :param srxn: a reaction object with stereo assignments
        :type srxn: Reaction
        :returns: a reflected reaction object
    """
    forw_tsg = automol.graph.reflect(srxn.forward_ts_graph)
    back_tsg = automol.graph.reflect(srxn.backward_ts_graph)
    srxn = Reaction(srxn.class_, forw_tsg, back_tsg,
                    srxn.reactants_keys, srxn.products_keys)
    return srxn
