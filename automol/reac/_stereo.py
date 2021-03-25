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


def conserved_atom_stereo_keys(rxn):
    """ Determine the atom stereo keys which are conserved by the reaction (are
    stereo centers for both reactants and products)

    :param rxn: a Reaction object
    :returns: The conserved stereo atom keys as a dictionary. Dictionary keys
        are atom keys for the reactants; dictionary values are atom keys for
        the products.
    :rtype: dict
    """
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph

    # Relabel the backward TSG to correspond to the atom ordering of the
    # forward TSG. At this point, if there is symmetry in the connectivity
    # graph, there is no guarantee that symmetric stereo centers won't be
    # swapped.
    iso_dct = automol.graph.isomorphism(
        ts.reverse(back_tsg), forw_tsg, stereo=False, dummy=False)
    back_tsg = automol.graph.relabel(back_tsg, iso_dct)

    ste_atm_keys = sorted(automol.graph.atom_stereo_keys(forw_tsg) &
                          automol.graph.atom_stereo_keys(back_tsg))
    print(iso_dct)
    print(ste_atm_keys)


def is_stereo_consistent(rxn):
    """ Determine if the reaction has mutually consistent stereo assignments
    for reactants and products.

    :param rxn: a Reaction object
    :returns: True if the stereo assignments are consistent; otherwise False
    :rtype: bool
    """


def substereomers(rxn):
    """ Obtain all stereomer combinations compatible with the current
    assignments (if any)

    :param rxn: a Reaction object
    :returns: a list of stereo-specific reaction objects
    :rtype: tuple[Reaction]
    """
    print(string(rxn))
    print('Hi')
