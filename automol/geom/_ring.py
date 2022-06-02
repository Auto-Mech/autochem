""" get lvl 4
"""

from phydat import phycon
import automol.graph
from automol.geom.base import central_angle
from automol.geom.base import from_subset


ATHRESH = 94.0 * phycon.DEG2RAD


def all_rings_angles_reasonable(geo, rings_atoms, thresh=ATHRESH):
    """ Assess if all angles of ring are cool
    """

    condition = True
    for ring_atoms in rings_atoms:
        condition = ring_angles_reasonable(
            geo, ring_atoms, thresh=thresh)

    return condition


def ring_angles_reasonable(geo, ring_atoms, thresh=ATHRESH):
    """ Assess whether any of the angles of the prospective rings of a geometry
        are bent at unphysical angles.

        :param geo: molecular geometry
        :type geo: automol.geom object
        :param rng_atoms: idxs for atoms inside rings
        :type rng_atoms: list
        :param thresh: angular threshold for large angles
        :type thresh: float
        :rtype: bool
    """

    condition = True
    for i, ring_atom in enumerate(ring_atoms):
        _atoms = [ring_atom, ring_atoms[i-1], ring_atoms[i-2]]
        cangle = central_angle(geo, *_atoms, degree=False)
        if cangle < thresh:
            condition = False
            break

    return condition


def ring_fragments_geometry(geo, rings_atoms=None, ngbs=None):
    """ Fragment out the ring and its neighbors in a geometry.(?)

        :param geo: molecular geometry
        :type geo: automol.geom object
    """

    gra = automol.geom.graph(geo)
    if rings_atoms is None:
        rings_atoms = automol.graph.rings_atom_keys(gra)
    if ngbs is None:
        ngbs = automol.graph.atoms_sorted_neighbor_atom_keys(gra)
    # print('geo in ring_frag_geo:', automol.geom.string(geo))
    # print('rings_atoms:', rings_atoms)
    # print('ngbs:', ngbs)

    ring_idxs = []
    ret = None
    for ring_atoms in rings_atoms:
        for ring_atom in ring_atoms:
            ring_ngbs = ngbs[ring_atom]
            if ring_atom not in ring_idxs:
                ring_idxs.append(ring_atom)
            for ngb in ring_ngbs:
                if ngb not in ring_idxs:
                    ring_idxs.append(ngb)
    # print('ring idxs:', ring_idxs)
    if ring_idxs:
        ret = from_subset(geo, ring_idxs)

    return ret
