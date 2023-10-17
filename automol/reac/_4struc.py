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
    product_mappings,
    product_structures,
    products_keys,
    reactant_graphs,
    reactant_mappings,
    reactant_structures,
    reactants_keys,
    set_structures,
    string,
    structure_type,
    ts_graph,
    ts_structure,
    without_stereo,
)
from automol.util import ZmatConv, dict_, zmat_conv


def with_structures(
    rxn: Reaction,
    struc_typ: str,
    rct_strucs=None,
    prd_strucs=None,
    rct_zcs=None,
    prd_zcs=None,
    log: bool = False,
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
    :param rct_zcs: Specify the reactant z-matrix conversions, defaults to None
    :type rct_zcs: List[ZmatConv], optional
    :param prd_zcs: Specify the product z-matrix conversions, defaults to None
    :type prd_zcs: List[ZmatConv], optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :return: The new reaction object
    :rtype: Reaction
    """
    orig_struc_typ = structure_type(rxn)
    no_strucs = all(x is None for x in [rct_strucs, prd_strucs, rct_zcs, prd_zcs])
    new_strucs = orig_struc_typ == struc_typ and not no_strucs

    # 1. Check that we are adding a valid structural type
    if struc_typ not in ("geom", "zmat"):
        raise ValueError(f"Requesting invalid structure type {struc_typ}")

    # If this Reaction object doesn't have structures, or we are updating a Reaction
    # object with new reagent structures, add the requested structure type
    if struc_typ == "geom" and (orig_struc_typ is None or new_strucs):
        rxn_ = _with_geom_structures(
            rxn,
            rct_geos=rct_strucs,
            prd_geos=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    elif struc_typ == "zmat" and (orig_struc_typ is None or new_strucs):
        rxn_ = _with_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    # If the Reaction already has this structure type, and no structures were passed in,
    # return as-is
    elif struc_typ == orig_struc_typ and no_strucs:
        rxn_ = rxn
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "zmat" and orig_struc_typ == "geom":
        rxn_ = _convert_geom_to_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "geom" and orig_struc_typ == "zmat":
        rxn_ = _convert_zmat_to_geom_structures(
            rxn,
            rct_geos=rct_strucs,
            prd_geos=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )

    return rxn_


def _with_geom_structures(
    rxn: Reaction,
    rct_geos=None,
    prd_geos=None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
    log: bool = False,
) -> Reaction:
    """Add geometry structures to a Reaction object

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_geos: Specify the reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param prd_geos: Specify the product geometries
    :type prd_geos: List[automol geom data structure]
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: List[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: List[ZmatConv], optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :returns: A new reaction object
    :rtype: Reaction
    """
    if log:
        print(f"Adding geometry structures to this Reaction:\n{string(rxn)}")
        print(f"Supplied reactant geometries (may be None):\n{rct_geos}")
        print(f"Supplied product geometries (may be None):\n{prd_geos}")

    if rct_geos is None:
        rct_geos = tuple(map(graph.geometry, reactant_graphs(rxn)))

    if prd_geos is None:
        prd_geos = tuple(map(graph.geometry, product_graphs(rxn)))

    tsg = ts_graph(rxn)
    geo_idx_dct = mapping(rxn, "T", "R")

    if log:
        print(f"Building geometry for TS graph:\n{tsg}")
        print(f"... using these reactant geometries:\n{rct_geos}")

    ts_geo = graph.ts_geometry_from_reactants(
        tsg, rct_geos, geo_idx_dct=geo_idx_dct, log=log
    )

    grxn = set_structures(
        rxn,
        ts_geo,
        rct_geos,
        prd_geos,
        struc_typ="geom",
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )
    return grxn


def _with_zmat_structures(
    rxn: Reaction,
    rct_zmas=None,
    prd_zmas=None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
    log: bool = False,
) -> Reaction:
    """Add structures to a Reaction object

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_zmas: Reactant z-matrices, defaults to None
    :type rct_zmas: List[automol zmat data structure]
    :param prd_zmas: Product z-matrices, defaults to None
    :type prd_zmas: List[automol zmat data structure]
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: List[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: List[ZmatConv], optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :returns: A new reaction object
    :rtype: Reaction
    """
    if log:
        print(f"Adding z-matrix structures to this Reaction:\n{string(rxn)}")
        print(f"Supplied reactant z-matrices (may be None):\n{rct_zmas}")
        print(f"Supplied product z-matrices (may be None):\n{prd_zmas}")

    def _reagent_geometries(zmas, zcs):
        if zmas is None:
            return None

        zcs = [None] * len(zmas) if zcs is None else zcs
        geos = [zmat.geometry(zma, zc_=zc) for zma, zc in zip(zmas, zcs)]
        return geos

    rct_geos = _reagent_geometries(rct_zmas, rct_zcs)
    prd_geos = _reagent_geometries(prd_zmas, prd_zcs)

    grxn = _with_geom_structures(rxn, rct_geos=rct_geos, prd_geos=prd_geos)
    zrxn = _convert_geom_to_zmat_structures(
        grxn, rct_zmas=rct_zmas, prd_zmas=prd_zmas, rct_zcs=rct_zcs, prd_zcs=prd_zcs
    )
    return zrxn


def _convert_geom_to_zmat_structures(
    rxn: Reaction,
    rct_zmas=None,
    prd_zmas=None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
    log: bool = False,
) -> Reaction:
    """Convert a reaction with 'geom' structures to 'zmat' structures

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_zmas: Reactant z-matrices, defaults to None
    :type rct_zmas: List[automol zmat data structure]
    :param prd_zmas: Product z-matrices, defaults to None
    :type prd_zmas: List[automol zmat data structure]
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: List[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: List[ZmatConv], optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :return: The new reaction object
    :rtype: Reaction
    """
    if log:
        print(f"Converting geom => zmat structures to this Reaction:\n{string(rxn)}")
        print(f"Supplied reactant z-matrices (may be None):\n{rct_zmas}")
        print(f"Supplied product z-matrices (may be None):\n{prd_zmas}")

    assert structure_type(rxn) == "geom"

    gtsg = ts_graph(rxn)
    ts_geo = ts_structure(rxn)

    # Convert the TS geometry
    ts_zma, ts_zc = geom.zmatrix_with_conversion_info(ts_geo, gra=gtsg)
    ztsg = graph.apply_zmatrix_conversion(gtsg, ts_zc)

    # Get the z-matrix for each reagent, along with the conversion info
    def _reagent_zmats_with_conversion_info(geos, gras, zmas, zcs):
        if zmas is None:
            # If no z-matrices were given create them from the geometry
            # If `zcs` is not None, try to impose the requested z-matrix conversion
            # (will fail if it can't be replicated)
            zcs = [None] * len(geos) if zcs is None else zcs
            zmas, zcs = zip(
                *(
                    geom.zmatrix_with_conversion_info(geo, gra=gra, zc_=zc)
                    for geo, gra, zc in zip(geos, gras, zcs)
                )
            )
        elif zcs is None:
            # If there were z-matrices but no conversions given, get the conversion info
            # from the z-matrices
            zcs = tuple(map(zmat.conversion_info, zmas))
        else:
            # If both `zmas` and `zcs` are specified, pass them along as given
            pass

        return zmas, zcs

    def _get_reagent_zmat_keys(rgts_gkeys, rgt_zcs, rgt_gkey_dcts):
        """
        :param rgts_gkeys: Keys for each reagent realtive to the TS geometry
        :param rgt_zcs: Z-matrix conversions for each reagent
        :param rgt_gkey_dcts: Mappings from TS gkeys onto each reagent's gkeys
        """
        rgts_zkeys = []
        for rgt_gkeys, rgt_zc, rgt_gkey_dct in zip(rgts_gkeys, rgt_zcs, rgt_gkey_dcts):
            # 1. {TS zkeys: TS gkeys} (extract the subset matching reagent gkeys)
            zc_ = zmat_conv.subset(ts_zc, rgt_gkeys, "geom")

            # 2. => {TS zkeys: reagent gkeys} (map TS gkeys into reagent gkeys)
            zc_ = zmat_conv.relabel(zc_, rgt_gkey_dct, "geom")

            # 3. => {TS zkeys: reagent zkeys} (get TS => reagent zkey mapping)
            tz_to_rz = zmat_conv.subset_isomorphism(zc_, rgt_zc)

            # 4. => [reagent TS zkeys] (get sorted TS zkeys for the reagent)
            rgt_zkeys = dict_.keys_sorted_by_value(tz_to_rz)

            # 5. Save the result
            rgts_zkeys.append(rgt_zkeys)

        return rgts_zkeys

    rct_zmas, rct_zcs = _reagent_zmats_with_conversion_info(
        geos=reactant_structures(rxn),
        gras=reactant_graphs(rxn),
        zmas=rct_zmas,
        zcs=rct_zcs,
    )
    prd_zmas, prd_zcs = _reagent_zmats_with_conversion_info(
        geos=product_structures(rxn),
        gras=product_graphs(rxn),
        zmas=prd_zmas,
        zcs=prd_zcs,
    )

    rcts_zkeys = _get_reagent_zmat_keys(
        rgts_gkeys=reactants_keys(rxn),
        rgt_zcs=rct_zcs,
        rgt_gkey_dcts=reactant_mappings(rxn, rev=True),
    )

    prds_zkeys = _get_reagent_zmat_keys(
        rgts_gkeys=products_keys(rxn),
        rgt_zcs=prd_zcs,
        rgt_gkey_dcts=product_mappings(rxn, rev=True),
    )

    return from_data(
        tsg=ztsg,
        rcts_keys=rcts_zkeys,
        prds_keys=prds_zkeys,
        cla=class_(rxn),
        ts_struc=ts_zma,
        rct_strucs=rct_zmas,
        prd_strucs=prd_zmas,
        struc_typ="zmat",
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )


def _convert_zmat_to_geom_structures(
    rxn: Reaction,
    rct_geos=None,
    prd_geos=None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
    log: bool = False,
) -> Reaction:
    """Convert a reaction with 'geom' structures to 'zmat' structures

    :param rxn: The reaction object
    :type rxn: Reaction
    :param rct_geos: Reactant geometries, defaults to None
    :type rct_geos: List[automol geom data structure]
    :param prd_geos: Product geometries, defaults to None
    :type prd_geos: List[automol geom data structure]
    :param rct_zcs: Z-matrix conversions for reactant structures, defaults to None
    :type rct_zcs: List[ZmatConv], optional
    :param prd_zcs: Z-matrix conversions for product structures, defaults to None
    :type prd_zcs: List[ZmatConv], optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :return: The new reaction object
    :rtype: Reaction
    """
    if log:
        print(f"Converting zmat => geom structures to this Reaction:\n{string(rxn)}")
        print(f"Supplied reactant geometries (may be None):\n{rct_geos}")
        print(f"Supplied product geometries (may be None):\n{prd_geos}")

    assert structure_type(rxn) == "zmat"

    ztsg = ts_graph(rxn)
    ts_zma = ts_structure(rxn)

    # Convert the TS geometry
    ts_geo, ts_zc = zmat.geometry_with_conversion_info(ts_zma)
    gtsg = graph.reverse_zmatrix_conversion(ztsg, ts_zc)

    # Get the geometry for each reagent, along with the conversion info
    def _reagent_geoms_with_conversion_info(zmas, geos, zcs):
        if geos is None:
            # If no geometries were given, get them from the z-matrix
            # Passing in a z-matrix conversion allows us to retrieve geometries with a
            # different atom ordering
            zcs = [None] * len(zmas) if zcs is None else zcs
            geos, zcs = zip(
                *(
                    zmat.geometry_with_conversion_info(zma, zc_=zc)
                    for zma, zc in zip(zmas, zcs)
                )
            )
        elif zcs is None:
            # If geometries were provided but no conversion info was specified, get it
            # from the z-matrices
            zcs = tuple(map(zmat.conversion_info, zmas))
        else:
            # If both `zmas` and `zcs` are specified, pass them along as given
            pass

        return geos, zcs

    def _get_reagent_geom_keys(rgts_zkeys, rgt_zcs, rgt_zkey_dcts):
        """
        :param rgts_zkeys: Keys for each reagent realtive to the TS z-matrix
        :param rgt_zcs: Z-matrix conversions for each reagent
        :param rgt_zkey_dcts: Mappings from TS zkeys onto each reagent's zkeys
        """
        rgts_gkeys = []
        for rgt_zkeys, rgt_zc, rgt_zkey_dct in zip(rgts_zkeys, rgt_zcs, rgt_zkey_dcts):
            # 1. {TS zkeys: TS gkeys} (extract the subset matching reagent zkeys)
            zc_ = zmat_conv.subset(ts_zc, rgt_zkeys, "zmat")

            # 2. => {reagent zkeys: TS gkeys} (map TS zkeys into reagent zkeys)
            zc_ = zmat_conv.relabel(zc_, rgt_zkey_dct, "zmat")

            # 3. {TS gkeys: reagent gkeys} (get TS => reagent gkeys mapping)
            tg_to_rz = zmat_conv.relabel_dict(zc_)
            rz_to_rg = zmat_conv.relabel_dict(rgt_zc, rev=True)
            tg_to_rg = dict_.compose(rz_to_rg, tg_to_rz)

            # 4. [reagent TS gkeys] (get sorted TS gkeys for the reagent)
            rgt_gkeys = dict_.keys_sorted_by_value(tg_to_rg)

            # 5. Save the result
            rgts_gkeys.append(rgt_gkeys)

        return rgts_gkeys

    rct_geos, rct_zcs = _reagent_geoms_with_conversion_info(
        zmas=reactant_structures(rxn),
        geos=rct_geos,
        zcs=rct_zcs,
    )
    prd_geos, prd_zcs = _reagent_geoms_with_conversion_info(
        zmas=product_structures(rxn),
        geos=prd_geos,
        zcs=prd_zcs,
    )

    rcts_gkeys = _get_reagent_geom_keys(
        rgts_zkeys=reactants_keys(rxn),
        rgt_zkey_dcts=reactant_mappings(rxn, rev=True),
        rgt_zcs=rct_zcs,
    )

    prds_gkeys = _get_reagent_geom_keys(
        rgts_zkeys=products_keys(rxn),
        rgt_zkey_dcts=product_mappings(rxn, rev=True),
        rgt_zcs=prd_zcs,
    )

    return from_data(
        tsg=gtsg,
        rcts_keys=rcts_gkeys,
        prds_keys=prds_gkeys,
        cla=class_(rxn),
        ts_struc=ts_geo,
        rct_strucs=rct_geos,
        prd_strucs=prd_geos,
        struc_typ="geom",
        ts_zc=ts_zc,
        rct_zcs=rct_zcs,
        prd_zcs=prd_zcs,
    )


def ts_geometry_from_reactants(
    rxn: Reaction,
    rct_geos,
    stereo=True,
    max_dist_err=0.2,
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