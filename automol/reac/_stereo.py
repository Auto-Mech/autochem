""" stereo functionality for reaction objects
"""
import automol.graph
from automol.reac._reac import Reaction
from automol.reac._reac import string
from automol.reac._reac import reverse
from automol.reac._reac import reactant_graphs
from automol.reac._reac import product_graphs
from automol.reac._reac import forming_bond_keys
from automol.reac._reac import breaking_bond_keys
from automol.graph import ts


def add_stereo_from_geometries(rxn, rct_geos, prd_geos):
    """ Add stereo assignments to this reaction object from geometries.

    :param rxn: a Reaction object
    :param rct_geos: the reactant geometries
    :param prd_geos: the product geometries
    :returns: a Reaction object with stereo assignments
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


def is_stereo_consistent(rxn):
    """ Determine if the reaction has mutually consistent stereo assignments
    for reactants and products.

    :param rxn: a Reaction object
    :returns: True if the stereo assignments are consistent; otherwise False
    :rtype: bool
    """
    # 1. Determine the isomorphism taking reactant atoms into product atoms
    iso_dct = automol.graph.full_isomorphism(
        rxn.forward_ts_graph, ts.reverse(rxn.backward_ts_graph))
    print(iso_dct)


def substereomers(rxn):
    """ Obtain all stereomer combinations compatible with the current
    assignments (if any)

    :param rxn: a Reaction object
    :returns: a list of stereo-specific reaction objects
    :rtype: tuple[Reaction]
    """
    print(string(rxn))
    print('Hi')
