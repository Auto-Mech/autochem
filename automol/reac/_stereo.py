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
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys

    key_dct = atom_mapping(rxn)

    forw_tsg = rxn.forward_ts_graph
    forw_ste_tsgs = ts.stereomers(forw_tsg)

    ste_rxns = []
    for forw_ste_tsg in forw_ste_tsgs:
        for back_ste_tsg in ts.compatible_reverse_stereomers(forw_ste_tsg):
            back_ste_tsg = automol.graph.relabel(back_ste_tsg, key_dct)

            ste_rxn = Reaction(rxn_cls, forw_ste_tsg, back_ste_tsg,
                               rcts_keys, prds_keys)
            ste_rxns.append(ste_rxn)

    ste_rxns = tuple(ste_rxns)
    return ste_rxns


def expand_product_stereo(rxn):
    """ Expand all possible stereo assignments for the products of this
    reaction, given a set of stereo assignments for the reactants. Stereo
    assignments for the products will be ignored.

    :param rxn: a reaction object
    :type rxn: Reaction
    :returns: a sequence reaction objects with stereo assignments
    :rtype: Reaction
    """
    rxn_cls = rxn.class_
    rcts_keys = rxn.reactants_keys
    prds_keys = rxn.products_keys

    key_dct = atom_mapping(rxn)

    forw_ste_tsg = rxn.forward_ts_graph

    ste_rxns = []
    for back_ste_tsg in ts.compatible_reverse_stereomers(forw_ste_tsg):
        back_ste_tsg = automol.graph.relabel(back_ste_tsg, key_dct)

        ste_rxn = Reaction(rxn_cls, forw_ste_tsg, back_ste_tsg,
                           rcts_keys, prds_keys)
        ste_rxns.append(ste_rxn)

    ste_rxns = tuple(ste_rxns)
    return ste_rxns
