""" working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numbers
from collections import abc
from typing import Dict, List

import more_itertools as mit
import numpy
from phydat import phycon

import automol.geom.base as geom_base
from automol import util
from automol.graph.base._0core import (
    atom_keys,
    atom_neighbor_atom_keys,
    atom_stereo_sorted_neighbor_keys,
    atoms_neighbor_atom_keys,
    atoms_sorted_neighbor_atom_keys,
    backbone_bond_keys,
    bond_stereo_sorted_neighbor_keys,
    bonds_neighbor_atom_keys,
    bonds_neighbor_bond_keys,
    explicit,
    relabel,
    ts_reacting_atom_keys,
    ts_transferring_atoms,
    without_dummy_atoms,
)
from automol.graph.base._2algo import (
    branch_atom_keys,
    ring_systems_atom_keys,
    rings_bond_keys,
)
from automol.graph.base._3kekule import rigid_planar_bond_keys
from automol.util.vec import Vector


# stereo parity evaluations
def geometry_atom_parity(gra, geo, atm_key, nkeys=None, geo_idx_dct=None):
    r""" Calculate an atom parity directly from a geometry

    Neighboring atom keys (`nkeys`) must be passed in as a priority-sorted
    list. If `None`, a local parity calculation will occur based on the
    atom keys in the molecular graph.

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
    :param atm_key: the atom key whose parity is being evaluated
    :type atm_key: int
    :param nkeys: the neighboring atom keys, pre-sorted by priority
    :type nkeys: list[int]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    nkeys = atom_stereo_sorted_neighbor_keys(gra, atm_key) if nkeys is None else nkeys

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    if len(nkeys) == 4:
        keys = nkeys
    else:
        assert len(nkeys) == 3
        keys = [atm_key] + list(nkeys)

    idxs = list(map(geo_idx_dct.__getitem__, keys))
    xyzs = geom_base.coordinates(geo, idxs=idxs)
    det_mat = numpy.ones((4, 4))
    det_mat[:, 1:] = xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.0  # for now, assume no four-atom planes
    par = bool(det_val > 0.0)
    return par


def geometry_bond_parity(gra, geo, bnd_key, bnd_nkeys=None, geo_idx_dct=None):
    r""" Calculate a bond parity directly from a geometry

    Neighboring bond keys (`bnd_nkeys`) must be passed in as a pair of
    priority-sorted lists corresponding to the first and second atoms in
    `bnd_key`. Note that the latter must be an *ordered list* in this case!
    If `None`, a local parity calculation will occur based on the atom keys
    in the molecular graph.

    Bond parity is defined as follows:

    For each atom in the double bond, find the heavy-atom neighbor with the
    higher canonical number. Although hydrogen atoms have higher canonical
    numbers, they are always given lowest priority.

    If the neighbors are cis to each other, the parity is False ('-').
    If the neighbors are trans to each other, the parity is True ('+').

        max     max      max     min
           \   /            \   /
            A=B              A=B
           /   \            /   \
        min     min      min     max

        cis              trans
        False            True
        '-'              '+'

    If one side only has a single neighbor, then it is compared with the
    maximum neighbor on the other side.

        max     nei      max
           \   /            \
            A=B              A=B
           /                /   \
        min              min     nei

        cis              trans
        False            True
        '-'              '+'

    If both sides have only single neighbors, then they are compared to
    each other.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param bnd_key: the bond key. If using `bnd_nkeys`, this must be an
        ordered list!
    :type bnd_key: list[int]
    :param bnd_nkeys: a pair of lists of neighboring keys for the first and
        second atoms in `bnd_key`, respectively.
    :type bnd_nkeys: list[list[int]]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."

    assert (
        isinstance(bnd_key, abc.Collection) and len(bnd_key) == 2
    ), f"{bnd_key} is not a valid bond key."
    key1, key2 = bnd_key

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    if bnd_nkeys is None:
        nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(gra, key1, key2)
    else:
        assert (
            isinstance(bnd_nkeys, abc.Collection)
            and len(bnd_nkeys) == 2
            and all(isinstance(nk, abc.Collection) for nk in bnd_nkeys)
        ), f"Bond neighbor keys should be a pair of lists: {bnd_nkeys}"
        nkey1s, nkey2s = bnd_nkeys

    idx1 = geo_idx_dct[key1]
    idx2 = geo_idx_dct[key2]
    nidx1 = geo_idx_dct[nkey1s[-1]]
    nidx2 = geo_idx_dct[nkey2s[-1]]

    (xyz1,) = geom_base.coordinates(geo, idxs=(idx1,))
    (xyz2,) = geom_base.coordinates(geo, idxs=(idx2,))
    (nxyz1,) = geom_base.coordinates(geo, idxs=(nidx1,))
    (nxyz2,) = geom_base.coordinates(geo, idxs=(nidx2,))

    bnd1_vec = numpy.subtract(nxyz1, xyz1)
    bnd2_vec = numpy.subtract(nxyz2, xyz2)

    dot_val = numpy.vdot(bnd1_vec, bnd2_vec)
    assert dot_val != 0.0  # for now, assume not collinear
    par = bool(dot_val < 0.0)
    return par


def geometry_local_parity(gra, geo, key, geo_idx_dct=None):
    """Calculate the local parity of an atom or bond

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param key: the atom or bond key whose parity is being evaluated
    :type key: int
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    if isinstance(key, numbers.Number):
        par = geometry_atom_parity(gra, geo, key, geo_idx_dct=geo_idx_dct)
    else:
        par = geometry_bond_parity(gra, geo, key, geo_idx_dct=geo_idx_dct)
    return par


def geometries_have_matching_parities(gra, geo1, geo2, keys, geo_idx_dct=None):
    """Check whether two geometries have matching parities at a list of sites

    Keys in list may be atom or bond keys.  Any stereo in the graph object
    gets ignored.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo1: the first molecular geometry
    :type geo1: automol geometry data structure
    :param geo2: the second molecular geometry
    :type geo2: automol geometry data structure
    :param keys: list of atom or bond keys for comparison sites
    :type keys: list
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :returns: true if they match, false if not
    """
    return all(
        (
            geometry_local_parity(gra, geo1, key, geo_idx_dct=geo_idx_dct)
            == geometry_local_parity(gra, geo2, key, geo_idx_dct=geo_idx_dct)
        )
        for key in keys
    )


def geometries_parity_mismatches(gra, geo1, geo2, keys, geo_idx_dct=None):
    """Check where two geometries have mismatched parities and return keys to
    those sites

    Keys in list may be atom or bond keys.  Any stereo in the graph object
    gets ignored.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo1: the first molecular geometry
    :type geo1: automol geometry data structure
    :param geo2: the second molecular geometry
    :type geo2: automol geometry data structure
    :param keys: list of atom or bond keys for comparison sites
    :type keys: list
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :returns: keys to sites at which they don't match
    """
    return tuple(
        key
        for key in keys
        if geometry_local_parity(gra, geo1, key, geo_idx_dct=geo_idx_dct)
        != geometry_local_parity(gra, geo2, key, geo_idx_dct=geo_idx_dct)
    )


def linear_segment_dummy_direction(
    gra, geo, seg_keys: List[int], geo_idx_dct: Dict[int, int] = None
) -> Vector:
    """Get a good direction for placing dummy atoms over a linear segment

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param seg_keys: Keys for the linear segment, sorted by connectivity so that the
        first and last keys are at the ends of th segment
    :type seg_keys: List[int]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: Dict[int, int], optional
    :return: The linear direction
    :rtype: Vector
    """
    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )
    gra = relabel(gra, geo_idx_dct)
    gra = without_dummy_atoms(gra)

    # 1. Find atoms book-ending either side of the linear segment
    nkeys_dct = atoms_sorted_neighbor_atom_keys(gra)
    tra_dct = ts_transferring_atoms(gra)
    end_keys = []
    for key in (seg_keys[0], seg_keys[-1]):
        nkeys = tra_dct[key] if key in tra_dct else nkeys_dct[key]
        end_key, *_ = (set(nkeys) - set(seg_keys)) - set(end_keys)
        end_keys.append(end_key)

    assert len(end_keys) == 2, "Sanity check"

    # 2. Determine the linear segment direction
    end_xyzs = geom_base.coordinates(geo, idxs=end_keys)
    seg_vec = util.vec.unit_norm(numpy.subtract(*end_xyzs))

    # 3. Find an auxiliary, non-parallel direction, if possible
    aux_vec = None
    tra_key = next((k for k in seg_keys if k in tra_dct), None)
    if tra_key is not None and len(nkeys_dct[tra_key]) > 2:
        # For a transferring atom, with spectator neighbors, find a direction that
        # avoids them as much as possible
        nkeys = set(nkeys_dct[tra_key]) - set(tra_dct[tra_key])
        nkey1, nkey2, *_ = nkeys
        nxyz1, nxyz2 = geom_base.coordinates(geo, idxs=(nkey1, nkey2))
        (tra_xyz,) = geom_base.coordinates(geo, idxs=(tra_key,))
        aux_vec = util.vec.unit_bisector(nxyz1, nxyz2, tra_xyz, outer=len(nkeys) < 3)
    else:
        # Otherwise, look for a non-parallel neighbor to one of the book-ending keys
        for key in end_keys:
            nkeys = set(nkeys_dct[key]) - set(seg_keys)
            for nkey in nkeys:
                xyz, nxyz = geom_base.coordinates(geo, idxs=(key, nkey))
                end_nvec = numpy.subtract(nxyz, xyz)
                if not util.vec.are_parallel(seg_vec, end_nvec, anti=True):
                    aux_vec = util.vec.unit_norm(end_nvec)
                    break

    # 4. Find the dummy direction
    if aux_vec is None:
        # If we don't have an auxiliary vector, choose an arbitrary vector perpendicular
        dummy_vec = util.vec.arbitrary_unit_perpendicular(seg_vec)
    else:
        # If we have an auxiliary vector, ortogonalize it against the segment direction
        # vector to get a nice perpendicular direction
        dummy_vec = util.vec.orthogonalize(seg_vec, aux_vec, normalize=True)

    return dummy_vec


# corrections
def geometry_correct_nonplanar_pi_bonds(
    gra, geo, geo_idx_dct=None, pert=5.0 * phycon.DEG2RAD
):
    """correct a geometry for non-planar pi-bonds

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param pert: Perturbation angle, in radians, to prevent perfect planarity
    :type pert: float
    """
    akeys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )
    gra = relabel(gra, geo_idx_dct)

    bnd_keys = rigid_planar_bond_keys(gra)

    for bkey in bnd_keys:
        key1, key2 = bkey

        # Figure out the correct angle to flip the stereo parity
        nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(gra, key1, key2)
        if nkey1s and nkey2s:
            nkey1 = nkey1s[-1]
            nkey2 = nkey2s[-1]
            dih_ang = geom_base.dihedral_angle(geo, nkey1, key1, key2, nkey2)
            if abs(dih_ang) < numpy.pi / 2.0:
                ang = pert - dih_ang
            else:
                ang = numpy.pi - pert - dih_ang

            geo = geometry_rotate_bond(gra, geo, [key1, key2], ang)

    return geo


def geometry_correct_linear_vinyls(
    gra, geo, geo_idx_dct=None, tol=2.0 * phycon.DEG2RAD
):
    """correct a geometry for linear vinyl groups

    :param gra: molecular graph
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

    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(atm_keys))}
    )

    bnd_keys = rigid_planar_bond_keys(gra)

    for bnd1_key in bnd_keys:
        for bnd2_key in bnbkeys_dct[bnd1_key]:
            (atm2_key,) = bnd1_key & bnd2_key
            (atm1_key,) = bnd1_key - {atm2_key}
            (atm3_key,) = bnd2_key - {atm2_key}

            atm1_idx = geo_idx_dct[atm1_key]
            atm2_idx = geo_idx_dct[atm2_key]
            atm3_idx = geo_idx_dct[atm3_key]

            ang = geom_base.central_angle(geo, atm1_idx, atm2_idx, atm3_idx)

            if numpy.abs(ang - numpy.pi) < tol:
                atm0_key = next(iter(bnakeys_dct[bnd1_key] - {atm3_key}), None)
                atm0_key = atm3_key if atm0_key is None else atm0_key
                atm0_idx = geo_idx_dct[atm0_key]

                xyzs = geom_base.coordinates(geo)

                atm0_xyz = xyzs[atm0_idx]
                atm1_xyz = xyzs[atm1_idx]
                atm2_xyz = xyzs[atm2_idx]

                rot_axis = util.vec.unit_perpendicular(
                    atm0_xyz, atm1_xyz, orig_xyz=atm2_xyz
                )

                rot_atm_keys = branch_atom_keys(gra, atm2_key, atm3_key)

                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = geom_base.rotate(
                    geo, rot_axis, numpy.pi / 3, orig_xyz=atm2_xyz, idxs=rot_idxs
                )

    return geo


def geometry_pseudorotate_atom(
    gra, geo, key, ang=numpy.pi, degree=False, geo_idx_dct=None
):
    r"""Pseudorotate an atom in a molecular geometry by a certain amount

    'Pseudorotate' here means to rotate all but two of the atom's neighbors, which can
    be used to invert/correct stereochemistry at an atom:

        1   2                                     1   2
         \ /                                       \ /
          C--3   = 1,4 pseudorotation by pi =>   3--C
          |                                         |
          4                                         4

    The two fixed atoms will be chosen to prevent the structural 'damage' from the
    rotation as much as possible. For example, atoms in rings will be favored to be
    fixed.

    If such a choice is not possible -- for example, if three or more neighbors are
    locked into connected rings -- then no geometry will be returned.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param key: The graph key of the atom to be rotated
    :type key: frozenset[int]
    :param ang: The angle of rotation (in radians, unless `degree = True`)
    :type ang: float
    :param degree: Is the angle of rotation in degrees?, default False
    :type degree: bool
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    ang = ang * phycon.DEG2RAD if degree else ang
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(atm_keys)} if geo_idx_dct is None else geo_idx_dct
    )

    # For simplicity, relabel the graph to match the geometry
    gra = relabel(gra, geo_idx_dct)

    rxn_keys = ts_reacting_atom_keys(gra)
    rsy_keys_lst = ring_systems_atom_keys(gra, lump_spiro=False)
    nkeys = atom_neighbor_atom_keys(gra, key)
    # Group together neighbors connected in a ring system
    nkey_sets = [nkeys & ks for ks in rsy_keys_lst if nkeys & ks]
    # Add the other neighbors as singletons
    nkey_sets.extend({k} for k in nkeys if not any(k in ks for ks in nkey_sets))
    # Put the biggest groups first
    nkey_sets = sorted(nkey_sets, key=len, reverse=True)
    # Move groups containing reacting atom keys last
    nkey_sets = sorted(nkey_sets, key=lambda g: bool(g & rxn_keys))

    # Now, find a pair of atoms to keep fixed
    found_pair = False
    for nkeys1, nkeys2 in mit.pairwise(nkey_sets + [set()]):
        if len(nkeys1) == 2 or len(nkeys1 | nkeys2) == 2:
            found_pair = True
            nkey1, nkey2, *_ = list(nkeys1) + list(nkeys2)
            break

    if found_pair:
        # Determine the rotational axis as the unit bisector between the fixed pair
        xyz, nxyz1, nxyz2 = geom_base.coordinates(geo, idxs=(key, nkey1, nkey2))
        rot_axis = util.vec.unit_bisector(nxyz1, nxyz2, orig_xyz=xyz)

        # Identify the remaining keys to be rotated
        rot_nkeys = nkeys - {nkey1, nkey2}
        rot_keys = set(
            itertools.chain(*(branch_atom_keys(gra, key, k) for k in rot_nkeys))
        )

        geo = geom_base.rotate(geo, rot_axis, ang, orig_xyz=xyz, idxs=rot_keys)
    else:
        geo = None

    return geo


def geometry_rotate_bond(gra, geo, key, ang, degree=False, geo_idx_dct=None):
    """Rotate a bond in a molecular geometry by a certain amount

    If no angle is passed in, the bond will be rotated to flip stereochemistry

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param key: The graph key of the bond to be rotated
    :type key: frozenset[int]
    :param ang: The angle of rotation (in radians, unless `degree = True`)
    :type ang: float
    :param degree: Is the angle of rotation in degrees?, default False
    :type degree: bool
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    """
    ang = ang * phycon.DEG2RAD if degree else ang
    akeys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )
    gra = relabel(gra, geo_idx_dct)

    key1, key2 = key
    xyzs = geom_base.coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]

    rot_axis = numpy.subtract(xyz2, xyz1)
    rot_keys = branch_atom_keys(gra, key1, key2)

    geo = geom_base.rotate(geo, rot_axis, ang, orig_xyz=xyz1, idxs=rot_keys)
    return geo


def geometry_dihedrals_near_value(
    gra,
    geo,
    ang,
    geo_idx_dct=None,
    tol=None,
    abs_=True,
    degree=False,
    rings=False,
):
    """Identify dihedrals of a certain value

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ang: The angle to check for
    :type ang: float
    :param tol: Tolerance for comparison (in radians, unless `degree = True`).
        Default is 5 degrees.
    :type tol: float
    :param abs_: Compare absolute values?
    :type abs_: bool
    :param rings: Include ring diherals?, defaults to False
    :type rings: bool, optional
    :returns: Quartets of dihedral keys matching this value
    :rtype: frozenset[tuple[int]]
    """
    ref_ang = ang * phycon.DEG2RAD if degree else ang
    if abs_:
        ref_ang = numpy.abs(ref_ang)
    tol = (
        5.0 * phycon.DEG2RAD
        if tol is None
        else (tol * phycon.DEG2RAD if degree else tol)
    )
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(atom_keys(gra)))}
    )

    nkeys_dct = atoms_neighbor_atom_keys(gra)
    bnd_keys = backbone_bond_keys(gra)
    if not rings:
        bnd_keys -= set(itertools.chain(*rings_bond_keys(gra)))

    dih_keys = []
    for atm2_key, atm3_key in bnd_keys:
        atm1_keys = nkeys_dct[atm2_key] - {atm3_key}
        atm4_keys = nkeys_dct[atm3_key] - {atm2_key}
        for atm1_key, atm4_key in itertools.product(atm1_keys, atm4_keys):
            if atm1_key != atm4_key:
                atm1_idx = geo_idx_dct[atm1_key]
                atm2_idx = geo_idx_dct[atm2_key]
                atm3_idx = geo_idx_dct[atm3_key]
                atm4_idx = geo_idx_dct[atm4_key]
                ang = geom_base.dihedral_angle(
                    geo, atm1_idx, atm2_idx, atm3_idx, atm4_idx
                )
                if abs_:
                    ang = numpy.abs(ang)
                if numpy.abs(ang - ref_ang) < tol:
                    dih_keys.append((atm1_key, atm2_key, atm3_key, atm4_key))
    return frozenset(dih_keys)
