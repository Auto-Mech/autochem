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


def created_atom_stereo_keys(rxn):
    """ Determine atoms which become stereogenic for the products.

    :param rxn: a Reaction object
    :returns: Atoms which become stereogenic for the products.
    :rtype: frozenset
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

    keys1 = automol.graph.atom_stereo_keys(forw_tsg)
    keys2 = automol.graph.atom_stereo_keys(back_tsg)

    cre_ste_atm_keys = frozenset(keys2 - keys1)

    return cre_ste_atm_keys


def destroyed_atom_stereo_keys(rxn):
    """ Determine atoms which are no longer stereogenic for the products.

    :param rxn: a Reaction object
    :returns: Atoms which are no longer stereogenic for the products.
    :rtype: frozenset
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

    keys1 = automol.graph.atom_stereo_keys(forw_tsg)
    keys2 = automol.graph.atom_stereo_keys(back_tsg)

    des_ste_atm_keys = frozenset(keys1 - keys2)

    return des_ste_atm_keys


def created_bond_stereo_keys(rxn):
    """ Determine bonds which become stereogenic for the products.

    :param rxn: a Reaction object
    :returns: Atoms which become stereogenic for the products.
    :rtype: frozenset
    """
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph

    # Relabel the backward TSG to correspond to the bond ordering of the
    # forward TSG. At this point, if there is symmetry in the connectivity
    # graph, there is no guarantee that symmetric stereo centers won't be
    # swapped.
    iso_dct = automol.graph.isomorphism(
        ts.reverse(back_tsg), forw_tsg, stereo=False, dummy=False)
    back_tsg = automol.graph.relabel(back_tsg, iso_dct)

    keys1 = automol.graph.bond_stereo_keys(forw_tsg)
    keys2 = automol.graph.bond_stereo_keys(back_tsg)

    cre_ste_bnd_keys = frozenset(keys2 - keys1)

    return cre_ste_bnd_keys


def destroyed_bond_stereo_keys(rxn):
    """ Determine bonds which are no longer stereogenic for the products.

    :param rxn: a Reaction object
    :returns: Atoms which are no longer stereogenic for the products.
    :rtype: frozenset
    """
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph

    # Relabel the backward TSG to correspond to the bond ordering of the
    # forward TSG. At this point, if there is symmetry in the connectivity
    # graph, there is no guarantee that symmetric stereo centers won't be
    # swapped.
    iso_dct = automol.graph.isomorphism(
        ts.reverse(back_tsg), forw_tsg, stereo=False, dummy=False)
    back_tsg = automol.graph.relabel(back_tsg, iso_dct)

    keys1 = automol.graph.bond_stereo_keys(forw_tsg)
    keys2 = automol.graph.bond_stereo_keys(back_tsg)

    des_ste_bnd_keys = frozenset(keys1 - keys2)

    return des_ste_bnd_keys


def expand_stereo_from_reactants(rxn):
    """ Determine possible products that are consistent with the stereo
    assignments of the reactants.
    """


def reaction_isomorphism(rxn, stereo=True, dummy=False):
    """ Determine the isomorphism taking mapping reactant atoms to product
    atoms.

    If the `stereo` flag is set to `True`, incompatible stereo parities on the
    reactants and products will result in a `None` for the isomorphism. This is
    a way of testing for stereo compatibility.

    :param rxn: a Reaction object
    :param stereo: Consider stereo?
    :type stereo: bool
    :param dummy: Consider dummy atoms?
    :type dummy: bool
    :returns: The isomorphism mapping reactant atoms to product atoms.
    :rtype: dict
    """
    forw_tsg = rxn.forward_ts_graph
    back_tsg = rxn.backward_ts_graph

    nca_keys1, nca_keys2 = nonconserved_atom_stereo_keys(rxn)
    ncb_keys1, ncb_keys2 = nonconserved_bond_stereo_keys(rxn)

    tsg1 = forw_tsg
    tsg2 = back_tsg

    if stereo:
        # Remove nonconserved stereo parities and relabel, in case there are
        # higher-order stereo assignments which depend on these.
        tsg1 = automol.graph.to_index_based_stereo(tsg1)
        tsg1 = automol.graph.remove_atom_stereo_parities(tsg1, nca_keys1)
        tsg1 = automol.graph.remove_bond_stereo_parities(tsg1, ncb_keys1)
        tsg1 = automol.graph.from_index_based_stereo(tsg1)

        # Remove nonconserved stereo parities, reverse the TS graph, and
        # relabel. This results in a TS graph whose stereo parities can be
        # directly compared to the one above.
        tsg2 = automol.graph.to_index_based_stereo(tsg2)
        tsg2 = automol.graph.remove_atom_stereo_parities(tsg2, nca_keys2)
        tsg2 = automol.graph.remove_bond_stereo_parities(tsg2, ncb_keys2)
        tsg2 = ts.reverse(tsg2)
        tsg2 = automol.graph.from_index_based_stereo(tsg2)

    iso_dct = automol.graph.isomorphism(tsg1, tsg2, stereo=stereo, dummy=dummy)
    return iso_dct


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
