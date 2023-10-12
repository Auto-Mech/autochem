""" TS geometries for specific reaction classes
"""
# from typing import List

from automol import geom, graph, zmat
from automol.reac._0core import (
    Reaction,
    class_,
    from_data,
    mapping,
    product_graphs,
    product_structures,
    reactant_graphs,
    reactant_structures,
    set_structures,
    structure_type,
    ts_graph,
    ts_structure,
    without_stereo,
)
from automol.util import dict_, dummy_conv


def with_geom_structures(rxn: Reaction, rct_geos=None, prd_geos=None) -> Reaction:
    """Add structures to a Reaction object

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_geos: Optionally, specify the reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param prd_geos: Optionally, specify the product geometries
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

    grxn = set_structures(rxn, ts_geo, rct_geos, prd_geos, "geom")
    return grxn


def convert_structures(rxn: Reaction, struc_typ: str = "zmat") -> Reaction:
    """Convert 'geom' structures to 'zmat', or convert 'zmat' structures to 'geom',
    updating graphs and reagent keys accordingly

    :param rxn: The reaction object
    :type rxn: Reaction
    :param struc_type: The new structure type ('zmat' or 'geom')
    :type struc_type: str
    :return: The new reaction object
    :rtype: Reaction
    """
    orig_struc_typ = structure_type(rxn)

    # Check that it has structures to begin with
    if orig_struc_typ not in ("geom", "zmat"):
        raise ValueError(
            f"Reaction with this structure type cannot be converted:\n{rxn}"
        )

    # Check that we are converting to a valid structure
    if struc_typ not in ("geom", "zmat"):
        raise ValueError(
            f"Cannot convert this reaction to structure type {struc_typ}:\n{rxn}"
        )

    # If this is already the structure type, return as-is
    if struc_typ == orig_struc_typ:
        return rxn

    # If we are converting geom => z-matrix, handle that case
    if orig_struc_typ == "geom" and struc_typ == "zmat":
        return convert_geom_to_zmat_structures(rxn)

    # If we are converting geom => z-matrix, handle that case
    if orig_struc_typ == "zmat" and struc_typ == "geom":
        raise NotImplementedError("zmat => geom conversion not yet implemented")


def convert_geom_to_zmat_structures(
    rxn: Reaction,
    rct_zmas=None,
    prd_zmas=None,
) -> Reaction:
    """Convert a reaction with 'geom' structures to 'zmat' structures

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_zmas: Optionally, specify the reactant z-matrices, defaults to None
    :type rct_zmas: List[automol zmat data structure], optional
    :param prd_zmas: Optionally, specify the product z-matrices, defaults to None
    :type prd_zmas: List[automol zmat data structure], optional
    :return: The new reaction object
    :rtype: Reaction
    """
    assert structure_type(rxn) == "geom"

    tsg0 = ts_graph(rxn)
    ts_geo = ts_structure(rxn)

    # Build the z-matrix for the TS
    ts_zma, ts_dc = geom.zmatrix_with_conversion_info(ts_geo, gra=tsg0)
    tsg = graph.apply_dummy_conversion(tsg0, ts_dc)

    # Loop over reactants and products to get z-matrices and updated keys for each
    rgt_zmas_lst = []
    rgts_keys_lst = []
    for typ in ("R", "P"):
        geos = reactant_structures(rxn) if typ == "R" else product_structures(rxn)
        gras = reactant_graphs(rxn) if typ == "R" else product_graphs(rxn)
        zmas = rct_zmas if typ == "R" else prd_zmas

        # 1. Build the z-matrices and get the individual dummy conversions
        if zmas is None:
            # If if no reagent z-matrices were provided, build them
            zmas, dcs = zip(
                *(geom.zmatrix_with_conversion_info(*gg) for gg in zip(geos, gras))
            )
        else:
            # Otherwise, get the dummy conversions for the provided reagent z-matrices
            dcs = tuple(map(zmat.dummy_conversion, zmas))

        # 2. Use the dummy conversion info to determine the updated reagent keys
        keys_lst = []
        start = 0
        for dc_ in dcs:
            # a. Shift the original reagent keys
            dc_ = dummy_conv.shift(dc_, start, original=True)

            # b. Update the starting index for the next reagent
            start = max(dummy_conv.true_atom_keys(dc_, original=True)) + 1

            # c. Map the original reagent keys onto the original TS keys
            dc_ = dummy_conv.relabel(dc_, mapping(rxn, typ, "T"), original=True)

            # d. Find the subgraph isomorphism for this reagent's final keys onto the
            # final TS keys
            iso_dct = dummy_conv.subgraph_isomorphism(ts_dc, dc_)
            keys = dict_.keys_sorted_by_value(iso_dct)

            # e. Save the result
            keys_lst.append(keys)

        # 3. Save the final results
        rgt_zmas_lst.append(zmas)
        rgts_keys_lst.append(keys_lst)

    rct_zmas, prd_zmas = rgt_zmas_lst
    rcts_keys, prds_keys = rgts_keys_lst

    return from_data(
        tsg=tsg,
        rcts_keys=rcts_keys,
        prds_keys=prds_keys,
        cla=class_(rxn),
        ts_struc=ts_zma,
        rct_strucs=rct_zmas,
        prd_strucs=prd_zmas,
        struc_typ="zmat",
    )


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
