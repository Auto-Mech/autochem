""" random functions that I am gonna distribute to other parts of automol
"""

from automol.zmat._zmat import geometry
from automol.geom._base import distance as _distance

# def z_atom_closest_to(zma, idx1, idx2, chk_idxs):
#     """
#     """
#     geo = automol.zmat.geometry(zma)
#     return def atom_closest_to(geo, idx1, idx2, chk_idxs)
#
#
# def atom_closest_to(geo, idx1, idx2, chk_idxs):
#     """ Check to see whether the radical atom is still closest to the bond
#         formation site.
#     """
#
#     atom_closest = True
#
#     dist1 = automol.geom.distance(geo, idx1, idx2)
#     for idx in chk_idxs:
#         if dist1 < automol.geom.distance(geo, idx1, idx) - 0.01:
#                 atom_closest = False
#
#     return atom_closest


def is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.
    """
    geo = geometry(zma)
    atom_closest = True
    for idx, _ in enumerate(geo):
        if idx < idx_rad:
            distance = _distance(geo, idx, idx_rad)
            if distance < bond_dist-0.01:
                atom_closest = False
    return atom_closest
