"""working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""

import itertools
import numbers

import numpy
from phydat import phycon

from ... import util
from ...geom import base as geom_base
from ._00core import (
    align_with_geometry,
    atoms_neighbor_atom_keys,
    backbone_bond_keys,
    bond_neighbor_atom_keys,
    ts_reactants_graph_without_stereo,
)
from ._02algo import (
    branch_atom_keys,
    rings_atom_keys,
    rings_bond_keys,
)
from ._03kekule import (
    rigid_planar_bonds,
    vinyl_radical_atom_bond_keys,
)
from ._05stereo import geometry_atom_parity, geometry_bond_parity


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
def geometry_correct_nonplanar_pi_bonds(
    gra,
    geo,
    geo_idx_dct=None,
    pert: float = 5.0 * phycon.DEG2RAD,
    excl_keys: frozenset[int] = frozenset(),
):
    """Correct a geometry for non-planar pi-bonds.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param pert: Perturbation angle, in radians, to prevent perfect planarity
    :param excl_keys: Atom keys whose bonds should be excluded (for linear atoms)
    """
    # Align the graph and the geometry keys/indices
    gra, geo, excl_keys, _, idx_dct = align_with_geometry(
        gra, geo, excl_keys, geo_idx_dct
    )

    rp_dct = rigid_planar_bonds(gra, min_ring_size=numpy.inf, excl_keys=excl_keys)

    for bkey, bnkeys in rp_dct.items():
        key1, key2 = sorted(bkey)
        nkey1, nkey2 = (nks[-1] for nks in bnkeys)
        dih_ang = geom_base.dihedral_angle(geo, nkey1, key1, key2, nkey2)

        # Rotate bonds that are closer to trans to 175 degrees
        if numpy.pi / 2 < abs(dih_ang) < 3 * numpy.pi / 2:
            ang = numpy.pi - pert - dih_ang
        # Rotate bonds that are closer to cis to 5 degrees
        else:
            ang = pert - dih_ang

        geo = geometry_rotate_bond(gra, geo, [key1, key2], ang)

    # Restore the original atom ordering of the geometry
    return geom_base.reorder(geo, idx_dct)


def geometry_correct_nonplanar_ring_pi_bond(gra, geo, bkey):
    """Correct a geometry for a non-planar ring pi-bond.

    :param gra: molecular graph
    :param geo: molecular geometry
    :param bkey: The bond key of the pi-bond to be corrected
    """
    rkeys_lst = list(map(set, rings_atom_keys(gra)))
    key1, key2 = sorted(bkey)
    rkeys = next(ks for ks in rkeys_lst if {key1, key2} <= ks)
    nkey1s, nkey2s = bond_neighbor_atom_keys(gra, key1, key2)
    rnkey = next(iter(nkey1s & rkeys), None)
    nkey1 = next(iter(nkey1s - rkeys), None)
    nkey2 = next(iter(nkey2s - rkeys), None)

    assert rnkey is not None, "No ring neighbor key found for pi-bond"
    xyz1, xyz2 = geom_base.coordinates(geo, idxs=(key1, key2))
    rot_axis = numpy.subtract(xyz2, xyz1)

    for key, nkey, xyz in [(key1, nkey1, xyz1), (key2, nkey2, xyz2)]:
        if nkey is not None:
            dih_ang = geom_base.dihedral_angle(geo, rnkey, key1, key2, nkey)
            rot_ang = numpy.pi - dih_ang
            rot_keys = branch_atom_keys(gra, key, nkey)
            geo = geom_base.rotate(geo, rot_axis, rot_ang, orig_xyz=xyz, idxs=rot_keys)

    return geo


def geometry_correct_linear_vinyls(
    gra,
    geo,
    geo_idx_dct=None,
    tol: float = 2.0 * phycon.DEG2RAD,
    excl_keys: frozenset[int] = frozenset(),
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
    :param excl_keys: Atom keys whose bonds should be excluded (linear atoms)
    """
    if geo is None:
        return None

    # Align the graph and the geometry keys/indices
    gra, geo, excl_keys, _, idx_dct = align_with_geometry(
        gra, geo, excl_keys, geo_idx_dct
    )

    gra = ts_reactants_graph_without_stereo(gra)
    rng_bkeys = set(itertools.chain(*rings_bond_keys(gra)))
    nkeys_dct = atoms_neighbor_atom_keys(gra, ts_=False)

    vin_dct = vinyl_radical_atom_bond_keys(gra)

    for key, bkey in vin_dct.items():
        if bkey not in rng_bkeys and not bkey & excl_keys:
            key2 = key
            (key1,) = bkey - {key}
            key3s = nkeys_dct[key] - {key, key1}
            if key3s:
                (key3,) = key3s

                ang = geom_base.central_angle(geo, key1, key2, key3)

                if numpy.abs(ang - numpy.pi) < tol:
                    key0 = next(iter(nkeys_dct[key1] - {key2}), None)
                    key0 = key3 if key0 is None else key0

                    xyz0, xyz1, xyz2 = geom_base.coordinates(
                        geo, idxs=(key0, key1, key2)
                    )

                    rot_axis = util.vector.unit_perpendicular(xyz0, xyz1, orig_xyz=xyz2)

                    rot_keys = branch_atom_keys(gra, key2, key3)

                    geo = geom_base.rotate(
                        geo, rot_axis, numpy.pi / 3, orig_xyz=xyz2, idxs=rot_keys
                    )

    # Restore the original atom ordering of the geometry
    return geom_base.reorder(geo, idx_dct)


def geometry_rotate_bond(gra, geo, key, ang, degree=False, geo_idx_dct=None):
    """Rotate a bond in a molecular geometry by a certain amount.

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
    # Align the graph and the geometry keys/indices
    gra, geo, (key,), _, idx_dct = align_with_geometry(gra, geo, (key,), geo_idx_dct)

    key1, key2 = key
    xyzs = geom_base.coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]

    rot_axis = numpy.subtract(xyz2, xyz1)
    rot_keys = branch_atom_keys(gra, key1, key2)

    geo = geom_base.rotate(geo, rot_axis, ang, orig_xyz=xyz1, idxs=rot_keys)

    # Restore the original atom ordering of the geometry
    return geom_base.reorder(geo, idx_dct)


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
    ref_ang = numpy.abs(ref_ang) if abs_ else ref_ang
    tol = (
        5.0 * phycon.DEG2RAD
        if tol is None
        else (tol * phycon.DEG2RAD if degree else tol)
    )
    # Align the graph and the geometry keys/indices
    gra, geo, _, key_dct, _ = align_with_geometry(gra, geo, (), geo_idx_dct)

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
                ang = geom_base.dihedral_angle(
                    geo, atm1_key, atm2_key, atm3_key, atm4_key
                )
                if abs_:
                    ang = numpy.abs(ang)
                if numpy.abs(ang - ref_ang) < tol:
                    dih_keys.append((atm1_key, atm2_key, atm3_key, atm4_key))

    return frozenset(util.translate(dih_keys, key_dct))
