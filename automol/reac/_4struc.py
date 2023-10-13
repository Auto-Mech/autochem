""" TS geometries for specific reaction classes
"""
from typing import List

from automol import geom, graph, zmat
from automol.reac._0core import (
    Reaction,
    class_,
    from_data,
    mapping,
    product_graphs,
    product_structures,
    products_conversion_info,
    reactant_graphs,
    reactant_structures,
    reactants_conversion_info,
    set_structures,
    structure_type,
    ts_conversion_info,
    ts_graph,
    ts_structure,
    without_stereo,
)
from automol.util import DummyConv, dict_, dummy_conv


def with_structures(
    rxn: Reaction,
    struc_typ: str,
    rct_strucs=None,
    prd_strucs=None,
    rct_dcs=None,
    prd_dcs=None,
) -> Reaction:
    """Convert 'geom' structures to 'zmat', or convert 'zmat' structures to 'geom',
    updating graphs and reagent keys accordingly

    :param rxn: The reaction object
    :type rxn: Reaction
    :param struc_type: The new structure type ('zmat' or 'geom')
    :type struc_type: str
    :param rct_strucs: Specify the reactant structures
    :type rct_strucs: List[automol geom or zmat data structure]
    :param prd_strucs: Specify the product structures
    :type prd_strucs: List[automol geom or zmat data structure]
    :param rct_dcs: Specify the reactant z-matrix conversions, defaults to None
    :type rct_dcs: List[DummyConv], optional
    :param prd_dcs: Specify the product z-matrix conversions, defaults to None
    :type prd_dcs: List[DummyConv], optional
    :return: The new reaction object
    :rtype: Reaction
    """
    orig_struc_typ = structure_type(rxn)
    no_strucs = rct_strucs is None and prd_strucs is None
    new_strucs = orig_struc_typ == struc_typ and not no_strucs

    # 1. Check that we are adding a valid structural type
    if struc_typ not in ("geom", "zmat"):
        raise ValueError(f"Requesting invalid structure type {struc_typ}")

    # If this Reaction object doesn't have structures, or we are updating a Reaction
    # object with new reagent structures, add the requested structure type
    if struc_typ == "geom" and (orig_struc_typ is None or new_strucs):
        print("Adding geom structures")
        rxn_ = _with_geom_structures(rxn, rct_geos=rct_strucs, prd_geos=prd_strucs)
    elif struc_typ == "zmat" and (orig_struc_typ is None or new_strucs):
        print("Adding zmat structures")
        rxn_ = _with_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_dcs=rct_dcs,
            prd_dcs=prd_dcs,
        )
    # If the Reaction already has this structure type, and no structures were passed in,
    # return as-is
    elif struc_typ == orig_struc_typ and no_strucs:
        print("No change")
        rxn_ = rxn
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "zmat" and orig_struc_typ == "geom":
        print("Converting geom structures to zmat structures")
        rxn_ = _convert_geom_to_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_dcs=rct_dcs,
            prd_dcs=prd_dcs,
        )
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "geom" and orig_struc_typ == "zmat":
        print("Converting zmat structures to geom structures")
        rxn_ = _convert_zmat_to_geom_structures(
            rxn,
            rct_geos=rct_strucs,
            prd_geos=prd_strucs,
            rct_dcs=rct_dcs,
            prd_dcs=prd_dcs,
        )

    return rxn_


def _with_geom_structures(rxn: Reaction, rct_geos=None, prd_geos=None) -> Reaction:
    """Add geometry structures to a Reaction object

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_geos: Specify the reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param prd_geos: Specify the product geometries
    :type prd_geos: List[automol geom data structure]
    :returns: A new reaction object
    :rtype: Reaction
    """
    if rct_geos is None:
        rct_geos = tuple(map(graph.geometry, reactant_graphs(rxn)))

    if prd_geos is None:
        prd_geos = tuple(map(graph.geometry, product_graphs(rxn)))

    tsg = ts_graph(rxn)
    geo_idx_dct = mapping(rxn, "T", "R")
    ts_geo = graph.ts_geometry_from_reactants(tsg, rct_geos, geo_idx_dct=geo_idx_dct)

    grxn = set_structures(
        rxn,
        ts_geo,
        rct_geos,
        prd_geos,
        struc_typ="geom",
    )
    return grxn


def _with_zmat_structures(
    rxn: Reaction,
    rct_zmas=None,
    prd_zmas=None,
    rct_dcs: List[DummyConv] = None,
    prd_dcs: List[DummyConv] = None,
) -> Reaction:
    """Add structures to a Reaction object

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_zmas: Optionally, specify the reactant z-matrices
    :type rct_zmas: List[automol zmat data structure]
    :param prd_zmas: Optionally, specify the product z-matrices
    :type prd_zmas: List[automol zmat data structure]
    :param rct_dcs: Optionally, specify dummy conversions for the z-matrices relative to
        the reactant keys, defaults to None
    :type rct_dcs: List[DummyConv], optional
    :param prd_dcs: Optionally, specify dummy conversions for the z-matrices relative to
        the product keys, defaults to None
    :type prd_dcs: List[DummyConv], optional
    :returns: A new reaction object
    :rtype: Reaction
    """

    def _reagent_geometries(zmas, dcs):
        if zmas is None:
            return None

        dcs = [None] * len(zmas) if dcs is None else dcs
        geos = [zmat.geometry(z, dc_=d) for z, d in zip(zmas, dcs)]
        return geos

    rct_geos = _reagent_geometries(rct_zmas, rct_dcs)
    prd_geos = _reagent_geometries(prd_zmas, prd_dcs)

    grxn = _with_geom_structures(rxn, rct_geos=rct_geos, prd_geos=prd_geos)
    zrxn = _convert_geom_to_zmat_structures(
        grxn, rct_zmas=rct_zmas, prd_zmas=prd_zmas, rct_dcs=rct_dcs, prd_dcs=prd_dcs
    )
    return zrxn


def _convert_geom_to_zmat_structures(
    rxn: Reaction,
    rct_zmas=None,
    prd_zmas=None,
    rct_dcs: List[DummyConv] = None,
    prd_dcs: List[DummyConv] = None,
) -> Reaction:
    """Convert a reaction with 'geom' structures to 'zmat' structures

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_zmas: Optionally, specify the reactant z-matrices, defaults to None
    :type rct_zmas: List[automol zmat data structure], optional
    :param prd_zmas: Optionally, specify the product z-matrices, defaults to None
    :type prd_zmas: List[automol zmat data structure], optional
    :param rct_dcs: Optionally, specify dummy conversions for the z-matrices relative to
        the reactant geometries, defaults to None
    :type rct_dcs: List[DummyConv], optional
    :param prd_dcs: Optionally, specify dummy conversions for the z-matrices relative to
        the product geometries, defaults to None
    :type prd_dcs: List[DummyConv], optional
    :return: The new reaction object
    :rtype: Reaction
    """
    assert structure_type(rxn) == "geom"

    tsg0 = ts_graph(rxn)
    ts_geo = ts_structure(rxn)

    # Convert the TS geometry
    ts_zma, ts_dc = geom.zmatrix_with_conversion_info(ts_geo, gra=tsg0)
    tsg = graph.apply_dummy_conversion(tsg0, ts_dc)

    # Get the z-matrices and updated keys for each set of reagents
    def _reagent_zmats_conversions_and_keys(geos, gras, zmas, dcs, typ):
        # 1. Convert the reagent geometries and get the conversion info
        if zmas is None:
            zmas, dcs = zip(
                *(geom.zmatrix_with_conversion_info(*gg) for gg in zip(geos, gras))
            )
        elif dcs is None:
            # Otherwise, get the dummy conversions for the provided reagent z-matrices
            dcs = tuple(map(zmat.dummy_conversion, zmas))

        # 2. Use the conversion info to determine the updated reagent keys
        keys_lst = []
        start = 0
        for dc_ in dcs:
            # a. Shift the reagent keys
            dc_ = dummy_conv.shift(dc_, start, original=True)

            # b. Update the starting index for the next reagent
            start = max(dummy_conv.real_keys(dc_, original=True)) + 1

            # c. Map the reagent keys onto the TS keys
            dc_ = dummy_conv.relabel(dc_, mapping(rxn, typ, "T"), original=True)

            # d. Find the subgraph isomorphism for this reagent's final keys onto the
            # final TS keys
            iso_dct = dummy_conv.subgraph_isomorphism(ts_dc, dc_)
            keys = dict_.keys_sorted_by_value(iso_dct)

            # e. Save the results
            keys_lst.append(keys)

        return zmas, dcs, keys_lst

    rct_zmas, rct_dcs, rcts_keys = _reagent_zmats_conversions_and_keys(
        geos=reactant_structures(rxn),
        gras=reactant_graphs(rxn),
        zmas=rct_zmas,
        dcs=rct_dcs,
        typ="R",
    )
    prd_zmas, prd_dcs, prds_keys = _reagent_zmats_conversions_and_keys(
        geos=product_structures(rxn),
        gras=product_graphs(rxn),
        zmas=prd_zmas,
        dcs=prd_dcs,
        typ="P",
    )

    return from_data(
        tsg=tsg,
        rcts_keys=rcts_keys,
        prds_keys=prds_keys,
        cla=class_(rxn),
        ts_struc=ts_zma,
        rct_strucs=rct_zmas,
        prd_strucs=prd_zmas,
        struc_typ="zmat",
        ts_dc=ts_dc,
        rct_dcs=rct_dcs,
        prd_dcs=prd_dcs,
    )


def _convert_zmat_to_geom_structures(
    rxn: Reaction,
    rct_geos=None,
    prd_geos=None,
    rct_dcs: List[DummyConv] = None,
    prd_dcs: List[DummyConv] = None,
) -> Reaction:
    """Convert a reaction with 'zmat' structures to 'geom' structures

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_geos: Optionally, specify the reactant geometries, defaults to None
    :type rct_geos: List[automol zmat data structure], optional
    :param prd_geos: Optionally, specify the product geometries, defaults to None
    :type prd_geos: List[automol zmat data structure], optional
    :param rct_dcs: Optionally, specify dummy conversions for the geometries relative to
        the reactant z-matrices, defaults to None
    :type rct_dcs: List[DummyConv], optional
    :param prd_dcs: Optionally, specify dummy conversions for the geometries relative to
        the product z-matrices, defaults to None
    :type prd_dcs: List[DummyConv], optional
    :return: The new reaction object
    :rtype: Reaction
    """
    assert structure_type(rxn) == "zmat"

    tsg0 = ts_graph(rxn)
    ts_zma = ts_structure(rxn)

    # Convert the TS z-matrix
    ts_dc = ts_conversion_info(rxn)
    ts_geo, ts_dc = zmat.geometry_with_conversion_info(ts_zma, dc_=ts_dc)
    tsg = graph.reverse_dummy_conversion(tsg0, ts_dc)

    def _reagent_geoms_conversions_and_keys(zmas, gras, geos, dcs, typ):
        # 1. Convert the reagent z-matrices and get the conversion info
        if geos is None:
            geos, dcs = zip(
                *(
                    zmat.geometry_with_conversion_info(z, dc_=d)
                    for z, d in zip(zmas, dcs)
                )
            )
        elif dcs is None:
            dcs = tuple(map(zmat.geometry_with_conversion_info, zmas))

        # Need to account for extra shifts due to additional dummy atoms in the TS

        # # 2. Use the conversion info to determine the updated reagent keys
        # keys_lst = []
        # start = 0
        # for dc_ in dcs:
        #     # a. Shift the reagent keys
        #     dc_ = dummy_conv.shift(dc_, start, original=False)

        #     # b. Update the starting index for the next reagent
        #     start0 = max(dummy_conv.real_keys(dc_, original=True)) + 1
        #     start = max(dc_)

        #     # c. Map the current reagent keys onto the current TS keys
        #     dc_ = dummy_conv.relabel(dc_, mapping(rxn, typ, "T"), original=False)

        #     # d. Get the relabelling from z-matrix TS keys back to geometry TS keys
        #     rel_dct = dummy_conv.relabel_dict(dc_)


def ts_geometry_from_reactants(
    rxn: Reaction,
    rct_geos,
    stereo=True,
    max_dist_err=2e-1,
    log=False,
):
    """Generate a TS geometry for this reaction object

    DEPRECATE THIS FUNCTION -- if someone wants to get the TS geometry directly, apart
    from a Reaction object, they should call graph.ts_geometry_from_reactants

    :param rxn: a Reaction object
    :type rxn: Reaction
    :param rct_geos: the reactant geometries
    :param stereo: Enforce correct stereochemistry?, default True
    :type stereo: bool, optional
    :param max_dist_err: The distance convergence threshold, in angstroms
    :type max_dist_err: float, optional
    :param log: Print optimization log?, defaults to False
    :type log: bool, optional
    :returns: the TS geometry
    """
    if not stereo:
        rxn = without_stereo(rxn)

    tsg = ts_graph(rxn)
    geo_idx_dct = mapping(rxn, "T", "R")
    ts_geo = graph.ts_geometry_from_reactants(
        tsg,
        rct_geos,
        geo_idx_dct=geo_idx_dct,
        max_dist_err=max_dist_err,
        log=log,
    )
    return ts_geo
