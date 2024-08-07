""" TS scanning functions
"""

import more_itertools as mit
import numpy

from .. import geom
from ..const import ReactionClass
from ..graph import ts
from ._0core import Reaction, class_, ts_graph, ts_structure
from ._1util import (
    hydrogen_migration_atom_keys,
    hydrogen_migration_might_dissociate,
    ring_forming_scission_chain,
)


def scan_coordinates(rxn: Reaction) -> tuple[int, ...]:
    """Determine the scan coordinates for a Reaction object

    :param rxn: A Reaction object
    :return: The scan coordinate
    """
    if class_(rxn) == ReactionClass.HYDROGEN_MIGRATION:
        (frm_bkey,) = ts.forming_bond_keys(ts_graph(rxn))
        scan_coo = tuple(sorted(frm_bkey))
        return (scan_coo,)
    if class_(rxn) in (ReactionClass.BETA_SCISSION, ReactionClass.RING_FORM_SCISSION):
        (frm_bkey,) = ts.breaking_bond_keys(ts_graph(rxn))
        scan_coo = tuple(sorted(frm_bkey))
        return (scan_coo,)

    return ()


def scan_values(rxn: Reaction) -> tuple[numpy.ndarray, ...]:
    """Determine the scan coordinate values for a Reaction object

    :param rxn: A Reaction object with structures assigned
    :return: The scan grid (one sequence for each coordinate)
    """
    geo = ts_structure(rxn)
    coos = scan_coordinates(rxn)
    dists = [geom.distance(geo, *c, angstrom=True) for c in coos]

    # 1. Form a grid for the first coordinate
    dist1 = dists[0]
    npoints1 = start1 = end1 = None
    if class_(rxn) == ReactionClass.HYDROGEN_MIGRATION:
        npoints1, start1, end1 = (16, dist1 + 0.05, dist1 - 0.55)
    if class_(rxn) == ReactionClass.BETA_SCISSION:
        npoints1, start1, end1 = (14, dist1 + 0.1, dist1 + 0.8)
    if class_(rxn) == ReactionClass.RING_FORM_SCISSION:
        npoints1, start1, end1 = (7, dist1 + 0.1, dist1 + 0.7)

    grid = numpy.linspace(start1, end1, npoints1)

    return (grid,)


def constraint_coordinates(rxn: Reaction) -> tuple[tuple[int, ...]]:
    """Determine the coordinates of a Reaction object to constrain during a scan

    :param rxn: A Reaction object
    :return: The constraint coordinates
    """
    if class_(rxn) == ReactionClass.HYDROGEN_MIGRATION:
        att_key, _, don_key, att_nkey = hydrogen_migration_atom_keys(rxn)
        diss_coo1 = tuple(sorted([att_key, att_nkey]))
        diss_coo2 = hydrogen_migration_might_dissociate(rxn, att_key, att_nkey, don_key)
        return (diss_coo1,) if diss_coo2 is None else (diss_coo1, diss_coo2)
    if class_(rxn) == ReactionClass.RING_FORM_SCISSION:
        chain_keys = ring_forming_scission_chain(rxn)
        ang_keys_lst = []
        ang_keys_lst.extend(sorted(mit.windowed(chain_keys[1:], 3)))
        ang_keys_lst.extend(sorted(mit.windowed(chain_keys, 4)))
        # Use standard coordinate keys
        rev_ang_keys_lst = list(map(tuple, map(reversed, ang_keys_lst)))
        ang_keys_lst = tuple(map(min, zip(ang_keys_lst, rev_ang_keys_lst)))
        return ang_keys_lst

    return ()
