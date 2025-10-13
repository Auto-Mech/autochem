"""TS geometries for specific reaction classes"""

from typing import List

from .. import geom, graph, zmat
from ..util import ZmatConv
from ._0core import (
    Reaction,
    apply_zmatrix_conversion,
    mapping,
    product_graphs,
    product_structures,
    products_conversion_info,
    reactant_graphs,
    reactant_structures,
    reactants_conversion_info,
    reset_conversion_info,
    reverse_without_recalculating,
    set_structures,
    string,
    structure_type,
    ts_conversion_info,
    ts_graph,
    ts_structure,
    undo_zmatrix_conversion,
    update_structures,
    without_dummy_atoms,
    without_structures,
)


def clean_ts_structure(rxn, ts_struc, struc_typ: str) -> object:
    """Clean a TS structure based on a Reaction object.

    :param rxn: The reaction object
    :param ts_struc: The TS structure to clean
    :param struc_typ: The structure type, "geom" or "zmat"
    :return: The cleaned TS structure
    """
    xrxn = with_structures(rxn, struc_typ=struc_typ)
    xrxn = update_structures(xrxn, ts_struc=ts_struc)
    grxn = with_structures(xrxn, "geom")
    ts_gra = ts_graph(grxn)
    ts_geo = ts_structure(grxn)
    ts_geo = graph.clean_ts_geometry(ts_gra, ts_geo=ts_geo)
    grxn = update_structures(grxn, ts_struc=ts_geo)
    xrxn = with_structures(grxn, struc_typ=struc_typ)
    return ts_structure(xrxn)


def with_structures(
    rxn: Reaction,
    struc_typ: str,
    ignore_zc: bool = False,
    rct_strucs=None,
    prd_strucs=None,
    rct_zcs=None,
    prd_zcs=None,
    log: bool = False,
) -> Reaction:
    """Convert 'geom' structures to 'zmat', or convert 'zmat' structures to 'geom',
    updating graphs and reagent keys accordingly.

    Dummy atoms will be added for 'zmat' structures, and removed for 'geom' structures,
    so there is no guarantee that the TS graph will stay the same, but all information
    will be kept internally consistent and consistent with the reactant/product inputs.

    :param rxn: The reaction object
    :type rxn: Reaction
    :param struc_type: The new structure type ('zmat' or 'geom')
    :type struc_type: str
    :param ignore_zc: Whether to ignore z-matrix conversion info.
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
    if struc_typ not in ("geom", "zmat", None):
        raise ValueError(f"Requesting invalid structure type {struc_typ}")

    # If there are no current structures, remove dummy atoms from the graph to avoid
    # problems
    if orig_struc_typ is None:
        rxn = without_dummy_atoms(rxn)

    # If requested to ignore z-matrix conversion info, re-set it
    if ignore_zc:
        rxn = reset_conversion_info(rxn)

    # If a structure type of `None` was requested, simply update with the structures
    # passed in
    if struc_typ is None:
        rxn = update_structures(
            rxn,
            rct_strucs=rct_strucs,
            prd_strucs=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
        )
    # If this Reaction object doesn't have structures, or we are updating a Reaction
    # object with new reagent structures, add the requested structure type
    elif struc_typ == "geom" and (orig_struc_typ is None or new_strucs):
        rxn = _with_geom_structures(
            rxn,
            rct_geos=rct_strucs,
            prd_geos=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    elif struc_typ == "zmat" and (orig_struc_typ is None or new_strucs):
        rxn = _with_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "zmat" and orig_struc_typ == "geom":
        rxn = _convert_geom_to_zmat_structures(
            rxn,
            rct_zmas=rct_strucs,
            prd_zmas=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    # If we are converting geom => z-matrix, handle that case
    elif struc_typ == "geom" and orig_struc_typ == "zmat":
        rxn = _convert_zmat_to_geom_structures(
            rxn,
            rct_geos=rct_strucs,
            prd_geos=prd_strucs,
            rct_zcs=rct_zcs,
            prd_zcs=prd_zcs,
            log=log,
        )
    # If the Reaction already has this structure type, and no structures were passed in,
    # return as-is

    return rxn


def reverse(rxn: Reaction, recalc: bool = True) -> Reaction:
    """Get the reaction object for the reverse reaction

    :param rxn: A reaction object
    :type rxn: Reaction
    :param recalc: Recalculate the TS structure while reversing? defaults to True
    :type recalc: bool, optional
    :returns: The reversed reaction object
    :rtype: Reaction
    """
    if not recalc or ts_structure(rxn) is None:
        return reverse_without_recalculating(rxn)

    # Convert to geometry structures before reversal, to avoid z-matrix issues
    rxn_ = without_structures(rxn)
    rxn_ = undo_zmatrix_conversion(rxn_)
    rev_rxn = reverse_without_recalculating(rxn_)

    return with_structures(
        rev_rxn,
        structure_type(rxn),
        rct_strucs=product_structures(rxn),
        prd_strucs=reactant_structures(rxn),
        rct_zcs=products_conversion_info(rxn),
        prd_zcs=reactants_conversion_info(rxn),
    )


# # helper functions
# # # structure addition/conversion
def _with_geom_structures(
    rxn: Reaction,
    rct_geos=None,
    prd_geos=None,
    rct_zcs: List[ZmatConv] = None,
    prd_zcs: List[ZmatConv] = None,
    log: bool = False,
) -> Reaction:
    """Add geometry structures to a Reaction object.

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

    ts_geo = graph.geometry(tsg, rct_geos=rct_geos, geo_idx_dct=geo_idx_dct, log=log)

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

    grxn = _with_geom_structures(rxn, rct_geos=rct_geos, prd_geos=prd_geos, log=log)
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
    """Convert a reaction with 'geom' structures to 'zmat' structures.

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

    # Read in pre-existing z-matrix conversion info, if there is some
    ts_zc = ts_conversion_info(rxn)
    rct_zcs = reactants_conversion_info(rxn) if rct_zcs is None else rct_zcs
    prd_zcs = products_conversion_info(rxn) if prd_zcs is None else prd_zcs

    gtsg = ts_graph(rxn)
    ts_geo = ts_structure(rxn)

    # Convert the TS geometry
    ts_zma, ts_zc = geom.zmatrix_with_conversion_info(ts_geo, gra=gtsg, zc_=ts_zc)

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

    rxn = without_structures(rxn)
    rxn = apply_zmatrix_conversion(rxn, ts_zc=ts_zc, rct_zcs=rct_zcs, prd_zcs=prd_zcs)
    return update_structures(
        rxn, ts_struc=ts_zma, rct_strucs=rct_zmas, prd_strucs=prd_zmas
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

    # Read in pre-existing z-matrix conversion info, if there is some
    ts_zc = ts_conversion_info(rxn)
    rct_zcs = reactants_conversion_info(rxn) if rct_zcs is None else rct_zcs
    prd_zcs = products_conversion_info(rxn) if prd_zcs is None else prd_zcs

    ts_zma = ts_structure(rxn)

    # Convert the TS geometry
    ts_geo, ts_zc = zmat.geometry_with_conversion_info(ts_zma, zc_=ts_zc)

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

    rxn = without_structures(rxn)
    rxn = undo_zmatrix_conversion(
        rxn, ts_zc=ts_zc, rct_zcs=rct_zcs, prd_zcs=prd_zcs, keep_info=True
    )
    return update_structures(
        rxn, ts_struc=ts_geo, rct_strucs=rct_geos, prd_strucs=prd_geos
    )
