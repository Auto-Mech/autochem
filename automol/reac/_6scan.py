""" TS scanning functions
"""

import numpy

from automol import geom
from automol.const import ReactionClass
from automol.graph import ts
from automol.reac._0core import Reaction, class_, ts_graph, ts_structure
from automol.reac._1util import (
    hydrogen_migration_atom_keys,
    hydrogen_migration_might_dissociate,
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
    npoints1, start1, end1 = {
        # ReactionClass.HYDROGEN_MIGRATION: (18, dist1 - 0.1, dist1 + 0.3),
        # ReactionClass.HYDROGEN_MIGRATION: (9, dist1 - 0.35, dist1 + 0.05),
        ReactionClass.HYDROGEN_MIGRATION: (16, dist1 + 0.05, dist1 - 0.35),
    }.get(class_(rxn), (None, None, None))

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

    return ()
