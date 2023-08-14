""" TS geometries for specific reaction classes
"""
import automol.geom.ts
from automol.reac._0core import Reaction, mapping, ts_graph, without_stereo


def ts_geometry(
    rxn: Reaction,
    rct_geos,
    stereo=True,
    max_dist_err=2e-1,
    debug_visualize=False,
    log=False,
):
    """Generate a TS geometry for this reaction object

    :param rxn: a Reaction object
    :type rxn: Reaction
    :param rct_geos: the reactant geometries
    :param stereo: Enforce correct stereochemistry?, default False
    :type stereo: bool, optional
    :param max_dist_err: The distance convergence threshold, in angstroms
    :type max_dist_err: float, optional
    :param debug_visualize: Prompt visualizations in a Jupyter notebook for debugging?
    :type debug_visualize: bool, optional
    :param log: Print optimization log?, defaults to False
    :type log: bool, optional
    :returns: the TS geometry
    """
    if not stereo:
        rxn = without_stereo(rxn)

    tsg = ts_graph(rxn)
    geo_idx_dct = mapping(rxn, "T", "R")
    ts_geo = automol.geom.ts.geometry_from_reactants(
        rct_geos,
        tsg,
        geo_idx_dct=geo_idx_dct,
        max_dist_err=max_dist_err,
        debug_visualize=debug_visualize,
        log=log,
    )
    return ts_geo
