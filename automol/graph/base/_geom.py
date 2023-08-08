""" working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numbers
from collections import abc

import more_itertools as mit
import numpy
from phydat import phycon

import automol.geom.base
from automol import util
from automol.graph.base._algo import branch_atom_keys, ring_systems_atom_keys
from automol.graph.base._core import (
    atom_keys,
    atom_stereo_sorted_neighbor_keys,
    atom_neighbor_atom_keys,
    atoms_neighbor_atom_keys,
    backbone_bond_keys,
    bond_stereo_sorted_neighbor_keys,
    bonds_neighbor_atom_keys,
    bonds_neighbor_bond_keys,
    explicit,
    relabel,
    ts_reacting_atom_keys,
    without_dummy_atoms,
)
from automol.graph.base._kekule import rigid_planar_bond_keys


# stereo parity evaluations
def geometry_atom_parity(
    gra, geo, atm_key, nkeys=None, geo_idx_dct=None, neg_hkeys=True
):
    r""" Calculate an atom parity directly from a geometry

    Neighboring atom keys (`nkeys`) must be passed in as a priority-sorted
    list. If `None`, a local parity calculation will occur based on the
    atom keys in the molecular graph. In this case, `neg_hkeys` can be used
    to determine whether to negate hydrogen keys, giving them lowest
    priority.

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
    :neg_hkeys: negate hydrogen keys, to match with InChI parities?
        Only has an effect if `nkeys` is `None`.
    :neg_hkeys: bool
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."
    gra = without_dummy_atoms(gra)

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(keys))}
    )
    xyzs = automol.geom.base.coordinates(geo)
    xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in keys}

    if nkeys is None:
        hkeys = atom_keys(gra, symb="H")
        sgn = -1 if neg_hkeys else +1

        pri_dct = {k: (k if k not in hkeys else sgn * k) for k in keys}

        nkeys = atom_stereo_sorted_neighbor_keys(gra, atm_key, pri_dct=pri_dct)

    # If there are only three groups, use the stereo atom itself as
    # the top apex of the tetrahedron.
    if len(nkeys) == 4:
        keys = nkeys
    else:
        assert len(nkeys) == 3
        keys = [atm_key] + list(nkeys)

    xyzs = list(map(list, map(xyz_dct.__getitem__, keys)))
    det_mat = numpy.ones((4, 4))
    det_mat[:, 1:] = xyzs
    det_val = numpy.linalg.det(det_mat)
    assert det_val != 0.0  # for now, assume no four-atom planes
    par = bool(det_val > 0.0)
    return par


def geometry_bond_parity(
    gra, geo, bnd_key, bnd_nkeys=None, geo_idx_dct=None, neg_hkeys=True
):
    r""" Calculate a bond parity directly from a geometry

    Neighboring bond keys (`bnd_nkeys`) must be passed in as a pair of
    priority-sorted lists corresponding to the first and second atoms in
    `bnd_key`. Note that the latter must be an *ordered list* in this case!
    If `None`, a local parity calculation will occur based on the atom keys
    in the molecular graph. In this case, `neg_hkeys` can be used to
    determine whether to negate hydrogen keys, giving them lowest priority.

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
    :neg_hkeys: negate hydrogen keys, to match with InChI parities?
        Only has an effect if `nkeys` is `None`.
    :neg_hkeys: bool
    """
    assert gra == explicit(
        gra
    ), "Explicit graph should be used when getting parities from geometry."
    gra = without_dummy_atoms(gra)

    assert (
        isinstance(bnd_key, abc.Collection) and len(bnd_key) == 2
    ), f"{bnd_key} is not a valid bond key."
    key1, key2 = bnd_key

    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(keys))}
    )
    xyzs = automol.geom.base.coordinates(geo)
    xyz_dct = {k: xyzs[geo_idx_dct[k]] for k in keys}

    if bnd_nkeys is None:
        hkeys = atom_keys(gra, symb="H")
        sgn = -1 if neg_hkeys else +1

        pri_dct = {k: (k if k not in hkeys else sgn * k) for k in keys}

        nkey1s, nkey2s = bond_stereo_sorted_neighbor_keys(
            gra, key1, key2, pri_dct=pri_dct
        )
    else:
        assert (
            isinstance(bnd_nkeys, abc.Collection)
            and len(bnd_nkeys) == 2
            and all(isinstance(nk, abc.Collection) for nk in bnd_nkeys)
        ), f"Bond neighbor keys should be a pair of lists: {bnd_nkeys}"
        nkey1s, nkey2s = bnd_nkeys

    xyz1 = xyz_dct[key1]
    xyz2 = xyz_dct[key2]
    nxyz1 = xyz_dct[nkey1s[-1]]
    nxyz2 = xyz_dct[nkey2s[-1]]

    bnd1_vec = numpy.subtract(nxyz1, xyz1)
    bnd2_vec = numpy.subtract(nxyz2, xyz2)

    dot_val = numpy.vdot(bnd1_vec, bnd2_vec)
    assert dot_val != 0.0  # for now, assume not collinear
    par = bool(dot_val < 0.0)
    return par


def geometry_local_parity(gra, geo, key, geo_idx_dct=None, neg_hkeys=True):
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
    :neg_hkeys: negate hydrogen keys, to match with InChI parities?
    :neg_hkeys: bool
    """
    if isinstance(key, numbers.Number):
        par = geometry_atom_parity(
            gra, geo, key, geo_idx_dct=geo_idx_dct, neg_hkeys=neg_hkeys
        )
    else:
        par = geometry_bond_parity(
            gra, geo, key, geo_idx_dct=geo_idx_dct, neg_hkeys=neg_hkeys
        )
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


# corrections
def linear_vinyl_corrected_geometry(
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

            ang = automol.geom.base.central_angle(geo, atm1_idx, atm2_idx, atm3_idx)

            if numpy.abs(ang - numpy.pi) < tol:
                atm0_key = next(iter(bnakeys_dct[bnd1_key] - {atm3_key}), None)
                atm0_key = atm3_key if atm0_key is None else atm0_key
                atm0_idx = geo_idx_dct[atm0_key]

                xyzs = automol.geom.base.coordinates(geo)

                atm0_xyz = xyzs[atm0_idx]
                atm1_xyz = xyzs[atm1_idx]
                atm2_xyz = xyzs[atm2_idx]

                rot_axis = util.vec.unit_perpendicular(
                    atm0_xyz, atm1_xyz, orig_xyz=atm2_xyz
                )

                rot_atm_keys = branch_atom_keys(gra, atm2_key, atm3_key)

                rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

                geo = automol.geom.rotate(
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
    nkey_sets.extend({k} for k in nkeys if not any(k in ks for ks in nkey_sets))
    # Put the biggest groups first
    nkey_sets = sorted(nkey_sets, key=len, reverse=True)
    # Move groups containing reacting atom keys last
    nkey_sets = sorted(nkey_sets, key=lambda g: bool(g & rxn_keys))

    # Now, find a pair of atoms to keep fixed
    found_pair = False
    for nkeys1, nkeys2 in mit.pairwise(nkey_sets):
        if len(nkeys1) == 2 or len(nkeys1 | nkeys2) == 2:
            found_pair = True
            nkey1, nkey2, *_ = nkeys1 | nkeys2
            break

    if found_pair:
        # Determine the rotational axis as the unit bisector between the fixed pair
        xyz, nxyz1, nxyz2 = automol.geom.base.coordinates(geo, idxs=(key, nkey1, nkey2))
        rot_axis = util.vec.unit_bisector(nxyz1, nxyz2, orig_xyz=xyz)

        # Identify the remaining keys to be rotated
        rot_nkeys = nkeys - {nkey1, nkey2}
        rot_keys = set(
            itertools.chain(*(branch_atom_keys(gra, key, k) for k in rot_nkeys))
        )

        geo = automol.geom.rotate(geo, rot_axis, ang, orig_xyz=xyz, idxs=rot_keys)
    else:
        geo = None

    return geo


def geometry_rotate_bond(gra, geo, key, ang, degree=False, geo_idx_dct=None):
    """Rotate a bond in a molecular geometry by a certain amount

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

    atm1_key, atm2_key = key
    xyzs = automol.geom.base.coordinates(geo)
    atm1_xyz = xyzs[atm1_key]
    atm2_xyz = xyzs[atm2_key]

    rot_axis = numpy.subtract(atm2_xyz, atm1_xyz)
    rot_atm_keys = branch_atom_keys(gra, atm1_key, atm2_key)
    rot_idxs = list(map(geo_idx_dct.__getitem__, rot_atm_keys))

    geo = automol.geom.rotate(geo, rot_axis, ang, orig_xyz=atm1_xyz, idxs=rot_idxs)
    return geo


def geometry_dihedrals_near_value(
    gra, geo, ang, geo_idx_dct=None, tol=None, abs_=True, degree=False
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
    bnd_keys = list(map(sorted, backbone_bond_keys(gra)))

    dih_keys = []
    for atm2_key, atm3_key in bnd_keys:
        atm1_keys = nkeys_dct[atm2_key] - {atm3_key}
        atm4_keys = nkeys_dct[atm3_key] - {atm2_key}
        for atm1_key, atm4_key in itertools.product(atm1_keys, atm4_keys):
            atm1_idx = geo_idx_dct[atm1_key]
            atm2_idx = geo_idx_dct[atm2_key]
            atm3_idx = geo_idx_dct[atm3_key]
            atm4_idx = geo_idx_dct[atm4_key]
            ang = automol.geom.base.dihedral_angle(
                geo, atm1_idx, atm2_idx, atm3_idx, atm4_idx
            )
            if abs_:
                ang = numpy.abs(ang)
            if numpy.abs(ang - ref_ang) < tol:
                dih_keys.append((atm1_key, atm2_key, atm3_key, atm4_key))
    return frozenset(dih_keys)


def perturb_geometry_planar_dihedrals(
    gra, geo, geo_idx_dct=None, ang=5.0 * phycon.DEG2RAD, degree=False
):
    """Remove symmetry from a geometry by perturbing planar dihedrals?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param ang: The angle of rotation
    :type ang: float
    :param degree: Is the angle of rotation in degrees?, default True
    :type degree: bool
    :returns: a geometry in which all planar dihedrals have been shifted by the
        given amount
    :rtype: automol geometry data structure
    """
    ang = ang * phycon.DEG2RAD if degree else ang
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(atom_keys(gra)))}
    )

    # Keep checking for planar dihedrals and rotating the corresponding bonds
    # until none are left
    cis_keys_lst = geometry_dihedrals_near_value(
        gra, geo, 0.0, geo_idx_dct=geo_idx_dct, tol=ang
    )
    trans_keys_lst = geometry_dihedrals_near_value(
        gra, geo, numpy.pi, geo_idx_dct=geo_idx_dct, tol=ang
    )

    # Otherwise, adjust the dihedral to give it a non-planar value
    dih_dct = {}
    dih_dct.update({k: 0.0 + ang for k in cis_keys_lst})
    dih_dct.update({k: numpy.pi - ang for k in trans_keys_lst})
    for keys, dih in dih_dct.items():
        idx, idx1, idx2, idx3 = list(map(geo_idx_dct.__getitem__, keys))
        geo = automol.geom.change_zmatrix_row_values(
            geo, idx=idx, idx1=idx1, idx2=idx2, idx3=idx3, dih=dih, degree=False
        )
    return geo
