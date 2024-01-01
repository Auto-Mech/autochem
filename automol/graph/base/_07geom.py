""" working with geometries

BEFORE ADDING ANYTHING, SEE IMPORT HIERARCHY IN __init__.py!!!!
"""
import itertools
import numbers

import more_itertools as mit
import numpy
from automol import util
from automol.geom import base as geom_base
from automol.graph.base._00core import (
    atom_keys,
    atom_neighbor_atom_keys,
    atoms_neighbor_atom_keys,
    backbone_bond_keys,
    relabel,
    ts_reactants_graph_without_stereo,
    ts_reacting_atom_keys,
    vinyl_radical_candidates,
)
from automol.graph.base._02algo import (
    branch_atom_keys,
    ring_systems_atom_keys,
    rings_bond_keys,
)
from automol.graph.base._03kekule import rigid_planar_bonds
from automol.graph.base._05stereo import geometry_atom_parity, geometry_bond_parity
from phydat import phycon


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

    rp_dct = rigid_planar_bonds(gra, min_ring_size=numpy.inf)

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
    keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    gra = ts_reactants_graph_without_stereo(gra)
    rng_bkeys = set(itertools.chain(*rings_bond_keys(gra)))
    nkeys_dct = atoms_neighbor_atom_keys(gra, ts_=False)

    vin_dct = vinyl_radical_candidates(gra, min_ncount=0)

    for key, bkey in vin_dct.items():
        if bkey not in rng_bkeys:
            key2 = key
            (key1,) = bkey - {key}
            key3s = nkeys_dct[key] - {key, key1}
            if key3s:
                (key3,) = key3s
                idx1, idx2, idx3 = tuple(map(geo_idx_dct.get, (key1, key2, key3)))

                ang = geom_base.central_angle(geo, idx1, idx2, idx3)

                if numpy.abs(ang - numpy.pi) < tol:
                    key0 = next(iter(nkeys_dct[key1] - {key2}), None)
                    key0 = key3 if key0 is None else key0
                    idx0 = geo_idx_dct[key0]

                    xyz0, xyz1, xyz2 = geom_base.coordinates(
                        geo, idxs=(idx0, idx1, idx2)
                    )

                    rot_axis = util.vector.unit_perpendicular(xyz0, xyz1, orig_xyz=xyz2)

                    rot_keys = branch_atom_keys(gra, key2, key3)
                    rot_idxs = list(map(geo_idx_dct.get, rot_keys))

                    geo = geom_base.rotate(
                        geo, rot_axis, numpy.pi / 3, orig_xyz=xyz2, idxs=rot_idxs
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
        rot_axis = util.vector.unit_bisector(nxyz1, nxyz2, orig_xyz=xyz)

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
    atm_keys = sorted(atom_keys(gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(atm_keys)} if geo_idx_dct is None else geo_idx_dct
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
