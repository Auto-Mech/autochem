""" extra Level 4 Z-Matrix functions
"""

from automol import geom
from automol.zmat._conv import geometry


def is_atom_closest_to_bond_atom(zma, idx_rad, bond_dist):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.

    TODO: DEPRECATE THIS FUNCTION -- as far as I can tell, it doesn't do what it is
    supposed to do
    """
    geo = geometry(zma)
    atom_closest = True
    for idx, _ in enumerate(geo):
        if idx < idx_rad:
            dist = geom.distance(geo, idx, idx_rad)
            if dist < bond_dist-0.01:
                atom_closest = False
    return atom_closest
