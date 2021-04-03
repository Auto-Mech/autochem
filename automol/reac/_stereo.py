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


def nonconserved_atom_stereo_keys(rxn):
    """ Determine the atom stereo keys which are not conserved by the reaction.

    That is, atoms which are stereo centers for the reactants but not for the
    products, or vice versa.

    :param rxn: a Reaction object
    :returns: Nonconserved stereo atoms for the reactants and the products.
    :rtype: (frozenset, frozenset)
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

    inv_iso_dct = dict(map(reversed, iso_dct.items()))

    keys1 = automol.graph.atom_stereo_keys(forw_tsg)
    keys2 = automol.graph.atom_stereo_keys(back_tsg)
    forw_keys = frozenset(keys1 - keys2)
    back_keys = frozenset(map(inv_iso_dct.__getitem__, keys2 - keys1))

    return forw_keys, back_keys


def nonconserved_bond_stereo_keys(rxn):
    """ Determine the bond stereo keys which are not conserved by the reaction.

    That is, bonds which are stereo centers for the reactants but not for the
    products, or vice versa.

    :param rxn: a Reaction object
    :returns: Nonconserved stereo bonds for the reactants and the products.
    :rtype: (frozenset, frozenset)
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

    inv_iso_dct = dict(map(reversed, iso_dct.items()))

    keys1 = automol.graph.bond_stereo_keys(forw_tsg)
    keys2 = automol.graph.bond_stereo_keys(back_tsg)
    forw_keys = frozenset(keys1 - keys2)
    back_keys = frozenset(map(inv_iso_dct.__getitem__, keys2 - keys1))

    return forw_keys, back_keys


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

    tsg1 = forw_tsg
    tsg2 = back_tsg

    if not stereo:
        tsg1 = automol.graph.without_stereo_parities(tsg1)
        tsg2 = automol.graph.without_stereo_parities(tsg2)

    if not dummy:
        tsg1 = automol.graph.without_dummy_atoms(tsg1)
        tsg2 = automol.graph.without_dummy_atoms(tsg2)

    tsg2 = automol.graph.to_index_based_stereo(tsg2)
    # Here, we need to eliminate all but the conserved stereo parities.
    tsg2 = ts.reverse(tsg2)
    tsg2 = automol.graph.from_index_based_stereo(tsg2)


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
