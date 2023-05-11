""" working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import numpy
from phydat import phycon
from automol import util
import automol.geom.base
from automol.graph.base._core import atom_keys
from automol.graph.base._core import bonds_neighbor_atom_keys
from automol.graph.base._core import bonds_neighbor_bond_keys
from automol.graph.base._core import explicit
from automol.graph.base._algo import branch
from automol.graph.base._kekule import rigid_planar_bond_keys


# stereo parity evaluations
def atom_parity_from_geometry(gra, geo, key, nkeys, geo_idx_dct=None):
    r""" Calculate the atom parity of a geometry directly, based on a sorted
        list of neighbor keys

        Atom parity is defined as follows:

        The four keys passed in are apices of a tetrahedron. Looking at 2, 3,
        and 4 from 1, they will either ascend in clockwise or counterclockwise
        order.

        If ascending in counterclockwise order, the parity is False ('-').
        If ascending in clockwise order, the parity is True ('+').

              2                   2
             /1\                 /1\
            3---4               4---3

            counterclockwise    clockwise
            False               True
            '-'                 '+'

        (Viewed looking down from 1)

        If only three keys are passed in, they will be treated as keys 2, 3,
        and 4 above and it will be assumed that there is a lone pair at 1.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :param key: the atom key whose parity is being evaluated
        :type key: int
        :param nkeys: the neighboring atom keys, pre-sorted by priority
        :type nkeys: list[int]
    """
    assert gra == explicit(gra), (
        "Explicit graph should be used when getting parities from geometry.")

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(keys))})
    xyzs = automol.geom.base.coordinates(geo)
    xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in keys}

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    if len(nkeys) == 4:
        keys = nkeys
    else:
        assert len(nkeys) == 3
        keys = [key] + list(nkeys)

    xyzs = list(map(list, map(xyz_dct.__getitem__, keys)))
    det_mat = numpy.ones((4, 4))
    det_mat[:, 1:] = xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.  # for now, assume no four-atom planes
    par = det_val > 0.
    return par


def bond_parity_from_geometry(gra, geo, key1, key2, nkey1s, nkey2s,
                              geo_idx_dct=None):
    r""" Calculate the bond parity of a geometry directly, based on a sorted
         list of neighbor keys

        Bond parity is defined as follows:

        For each atom in the double bond, find the heavy-atom neighbor with the
        higher canonical number. Although hydrogen atoms have higher canonical
        numbers, they are always given lowest priority.

        If the neighbors are cis to each other, the parity is False ('-').
        If the neighbors are trans to each other, the parity is True ('+').

            max    max      max    min
              \   /           \   /
               A=B             A=B
              /   \           /   \
            min    min      min    max

            cis             trans
            False           True
            '-'             '+'

        If one side only has a single neighbor, then it is compared with the
        maximum neighbor on the other side.

            max    nei      max
              \   /           \
               A=B             A=B
              /               /   \
            min             min    nei

            cis             trans
            False           True
            '-'             '+'

        If both sides have only single neighbors, then they are compared to
        each other.

        :param gra: molecular graph
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :param key1: the first atom key in the bond being evaluated
        :type key1: int
        :param key2: the second atom key in the bond being evaluated
        :type key2: int
        :param nkey1s: neighbors of the first atom key, pre-sorted by priority
        :type nkey1s: list[int]
        :param nkey2s: neighbors of the second atom key, pre-sorted by priority
        :type nkey2s: list[int]
    """
    assert gra == explicit(gra), (
        "Explicit graph should be used when getting parities from geometry.")

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(keys))})
    xyzs = automol.geom.base.coordinates(geo)
    xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in keys}

    xyz1 = xyz_dct[key1]
    xyz2 = xyz_dct[key2]
    nxyz1 = xyz_dct[nkey1s[-1]]
    nxyz2 = xyz_dct[nkey2s[-1]]

    bnd1_vec = numpy.subtract(nxyz1, xyz1)
    bnd2_vec = numpy.subtract(nxyz2, xyz2)

    dot_val = numpy.vdot(bnd1_vec, bnd2_vec)
    assert dot_val != 0.    # for now, assume not collinear
    par = dot_val < 0.
    return par


# corrections
def linear_vinyl_corrected_geometry(gra, geo, geo_idx_dct=None,
                                    tol=2.*phycon.DEG2RAD):
    """ correct a geometry for linear vinyl groups

        :param gra: molecular graph with stereo parities
        :type gra: automol graph data structure
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param geo_idx_dct: If they don't already match, specify which graph
            keys correspond to which geometry indices.
        :type geo_idx_dct: dict[int: int]
        :param tol: tolerance of bond angle(s) for determing linearity
        :type tol: float
    """
    atm_keys = atom_keys(gra)
    bnakeys_dct = bonds_neighbor_atom_keys(gra)
    bnbkeys_dct = bonds_neighbor_bond_keys(gra)

    geo_idx_dct = (geo_idx_dct if geo_idx_dct is not None
                   else {k: i for i, k in enumerate(sorted(atm_keys))})

    bnd_keys = rigid_planar_bond_keys(gra)

    for bnd1_key in bnd_keys:
        for bnd2_key in bnbkeys_dct[bnd1_key]:
            atm2_key, = bnd1_key & bnd2_key
            atm1_key, = bnd1_key - {atm2_key}
            atm3_key, = bnd2_key - {atm2_key}

            atm1_idx = geo_idx_dct[atm1_key]
            atm2_idx = geo_idx_dct[atm2_key]
            atm3_idx = geo_idx_dct[atm3_key]

            ang = automol.geom.base.central_angle(
                geo, atm1_idx, atm2_idx, atm3_idx)

            if numpy.abs(ang - numpy.pi) < tol:
                atm0_key = next(iter(bnakeys_dct[bnd1_key] - {atm3_key}), None)
                atm0_key = atm3_key if atm0_key is None else atm0_key
                atm0_idx = geo_idx_dct[atm0_key]

                xyzs = automol.geom.base.coordinates(geo)

                atm0_xyz = xyzs[atm0_idx]
                atm1_xyz = xyzs[atm1_idx]
                atm2_xyz = xyzs[atm2_idx]

                rot_axis = util.vec.unit_perpendicular(atm0_xyz, atm1_xyz,
                                                       orig_xyz=atm2_xyz)

                rot_atm_keys = atom_keys(
                    branch(gra, atm2_key, {atm2_key, atm3_key}))

                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = automol.geom.rotate(
                    geo, rot_axis, numpy.pi/3,
                    orig_xyz=atm2_xyz, idxs=rot_idxs)

    return geo
