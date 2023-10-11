""" TS geometries for specific reaction classes
"""
import automol.graph
from automol.reac._0core import (
    Reaction,
    mapping,
    product_graphs,
    reactant_graphs,
    set_structures,
    ts_graph,
    without_stereo,
)


def with_geom_structures(rxn, rct_geos=None, prd_geos=None) -> Reaction:
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
        rct_geos = tuple(map(automol.graph.geometry, reactant_graphs(rxn)))

    if prd_geos is None:
        prd_geos = tuple(map(automol.graph.geometry, product_graphs(rxn)))

    tsg = ts_graph(rxn)
    geo_idx_dct = mapping(rxn, "T", "R")
    ts_geo = automol.graph.ts_geometry_from_reactants(
        tsg, rct_geos, geo_idx_dct=geo_idx_dct
    )

    grxn = set_structures(rxn, ts_geo, rct_geos, prd_geos, "geom")
    return grxn


def ts_geometry_from_reactants(
    rxn: Reaction,
    rct_geos,
    stereo=True,
    max_dist_err=2e-1,
    log=False,
):
    """Generate a TS geometry for this reaction object

    DEPRECATE THIS FUNCTION -- if someone wants to get the TS geometry directly, apart
    from a Reaction object, they should call automol.graph.ts_geometry_from_reactants

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
    ts_geo = automol.graph.ts_geometry_from_reactants(
        tsg,
        rct_geos,
        geo_idx_dct=geo_idx_dct,
        max_dist_err=max_dist_err,
        log=log,
    )
    return ts_geo
