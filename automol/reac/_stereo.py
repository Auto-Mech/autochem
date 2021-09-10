""" stereo functionality for reaction objects
"""
import automol.graph
from automol.reac._reac import Reaction
from automol.reac._reac import reverse
from automol.reac._reac import atom_mapping
from automol.reac._reac import reactant_graphs
from automol.reac._reac import product_graphs
from automol.reac._reac import forming_bond_keys
from automol.reac._reac import breaking_bond_keys
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


def expand_stereo(rxn):
    """ Expand all possible stereo assignments for the reactants and products
    of this reaction. Only includes possibilities that are mutually consistent
    with each other.

    :param rxn: a reaction object
    :type rxn: Reaction
    :returns: a sequence reaction objects with stereo assignments
    :rtype: Reaction
    """
    rxn_cls = rxn.class_
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys

    key_dct = atom_mapping(rxn)

    forw_ste_tsgs = ts.stereomers(forw_tsg)

    srxns = []
    for forw_ste_tsg in forw_ste_tsgs:
        for back_ste_tsg in ts.compatible_reverse_stereomers(forw_ste_tsg):
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


def expand_product_stereo(srxn):
    """ Expand all possible stereo assignments for the products of this
    reaction, given a set of stereo assignments for the reactants. Stereo
    assignments for the products will be ignored.

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: a sequence reaction objects with stereo assignments
    :rtype: Reaction
    """
    rxn_cls = srxn.class_
    back_tsg = srxn.backward_ts_graph
    rcts_keys = srxn.reactants_keys
    prds_keys = srxn.products_keys

    key_dct = atom_mapping(srxn)

    forw_ste_tsg = srxn.forward_ts_graph

    srxns = []
    for back_ste_tsg in ts.compatible_reverse_stereomers(forw_ste_tsg):
        back_ste_tsg = automol.graph.relabel(back_ste_tsg, key_dct)

        # But for dummy atoms, we could just do the conversion directly, but
        # this avoids loss of dummy atoms from the products graph.
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


def is_stereo_consistent(srxn):
    """ Does this reaction have consistent stereo assignments for reactants and
    products?

    :param srxn: a reaction object with stereo assignments
    :type srxn: Reaction
    :returns: True if stereo assignments are consisent, False if not
    :rtype: bool
    """
    back_tsg = srxn.backward_ts_graph

    srxns = expand_product_stereo(srxn)
    return any(
        automol.graph.isomorphism(back_tsg, s.backward_ts_graph, stereo=True)
        for s in srxns)
