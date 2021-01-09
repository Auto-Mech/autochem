""" generate ts geometries
"""
import automol.graph
import automol.convert.geom


def from_reactant_geometries(geos, bnd_dist_dct, bnd_ang_dct=None,
                             angstrom=True, degree=True):
    """ generate a ts geometry from reactant geometries

    :param geos: the reactant geometries
    :param bnd_dist_dct: a dictionary of desired bond distances for breaking
        and forming (or other) bonds; the keys are pairs of atoms, indexed in
        the order they appear in `geos`
    :param bnd_ang_dct: a dictionary of desired bond angles, presumably at the
        reaction site but could be anywhere
    """
    gras = list(map(automol.convert.geom.connectivity_graph, geos))
    gra, _ = automol.graph.standard_keys_for_sequence(gras)
