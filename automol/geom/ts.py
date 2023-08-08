""" generate ts geometries
"""
import numpy
import scipy
from phydat import phycon

import automol.graph.base
from automol.geom._conv import connectivity_graph
from automol.geom.base import (
    coordinates,
    count,
    distance,
    from_data,
    rotate,
    set_coordinates,
    symbols,
    translate,
)
from automol.util import dict_, vec


def geometry_from_reactants(
    geos, tsg, geo_idx_dct=None, fdist_factor=1.1, bdist_factor=0.9
):
    """Generate a TS geometry from reactants

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :return: TS geometry
    :rtype: automol geom data structure
    """
    # 1. If there are multiple reactants, combine the geometries
    if len(geos) == 1:
        (geo,) = geos
    else:
        assert (
            len(geos) == 2 and len(automol.graph.base.ts.forming_bond_keys(tsg)) == 1
        ), "Generating a TS geometry for this case is not implemented:\n{tsg}"
        geo = join_at_forming_bond(
            geos, tsg, geo_idx_dct=geo_idx_dct, fdist_factor=fdist_factor
        )

    # 2. Correct the stereochemistry against the TS graph, so it is consistent with both
    # reactants and products
    ts_geo = geo
    # ts_geo = automol.graph.stereo_corrected_geometry(tsg, geo,
    # geo_idx_dct=geo_idx_dct)

    print(bdist_factor)

    return ts_geo
    # This
    # # 1. Evaluate the distances, updating them for the TS graph
    # dist_dct = distances(
    #     geos,
    #     tsg,
    #     geo_idx_dct=geo_idx_dct,
    #     fdist_factor=fdist_factor,
    #     bdist_factor=bdist_factor,
    # )


def join_at_forming_bond(geos, tsg, geo_idx_dct=None, fdist_factor=1.1):
    """Join two reactant geometries at a single forming bond

    If stereochemistry is present in TS graph, the join will take stereochemistry into
    account.

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :returns: The joined geometry
    :rtype: automol geom data structure
    """
    akeys = sorted(automol.graph.base.atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    tsg = automol.graph.base.relabel(tsg, geo_idx_dct)
    frm_keys = automol.graph.base.ts.forming_bond_keys(tsg)
    assert len(frm_keys) == 1, f"This only works for one forming bond\n{tsg}"
    (frm_key,) = frm_keys

    # To begin, assume both ends are pi-forming
    assert len(geos) == 2, f"This requires two reactants, but {len(geos)} were given"
    len1, len2 = map(count, geos)
    idxs1 = tuple(range(len1))
    idxs2 = tuple(range(len1, len1 + len2))
    (fidx1,) = frm_key & set(idxs1)
    (fidx2,) = frm_key & set(idxs2)
    geo = sum(geos, ())
    # 1. Find the atoms on either side of the forming bond
    fxyz1 = coordinates(geo)[fidx1]
    fxyz2 = coordinates(geo)[fidx2]
    # 2. Translate to place them right on top of each other
    geo = translate(geo, numpy.subtract(fxyz1, fxyz2), idxs=idxs2)
    fxyz2 = coordinates(geo)[fidx2]
    # 3. Find the reacting electron directions
    rvec1, rvec2 = reacting_electron_directions(geo, tsg, fidx1, fidx2)
    # 4. Rotate to align them antiparallel
    rot_ang = vec.angle(rvec1, numpy.negative(rvec2))
    rot_axis = vec.unit_perpendicular(rvec1, rvec2)
    geo = rotate(geo, rot_axis, rot_ang, orig_xyz=fxyz2, idxs=idxs2)
    # 5. Translate the second reagent to the appropriate distance away
    fdist = automol.graph.base.ts.heuristic_bond_distance(
        tsg, fidx1, fidx2, fdist_factor=fdist_factor, angstrom=True
    )
    fvec = numpy.multiply(rvec1, fdist)
    geo = translate(geo, fvec, idxs=idxs2, angstrom=True)
    return geo


def reacting_electron_directions(geo, tsg, key1, key2) -> (vec.Vector, vec.Vector):
    """Identify the direction of a reacting electron in a forming bond

    If stereochemistry is present in the TS graph, the directions will be consistent
    with the correct stereochemistry.

    :param geo: Reactants geometry, aligned to the TS graph
    :type geo: automol geom data structure
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key1: Key for the first atom in the forming bond
    :type key1: int
    :param key2: Key for the second atom in the forming bond
    :type key2: int
    :returns: A vector for first atom, and a vector for the second atom
    :rtype: (vec.Vector, vec.Vector)
    """
    frm_key = frozenset({key1, key2})
    assert frm_key in automol.graph.base.ts.forming_bond_keys(tsg)

    if automol.graph.base.has_stereo(tsg):
        tsg = automol.graph.base.to_local_stereo(tsg)

    no_rxn_bnd_gra = automol.graph.base.without_reacting_bonds(tsg)
    rcts_gra = automol.graph.base.ts.reactants_graph(tsg)
    vin_dct = automol.graph.base.vinyl_radical_atom_bond_keys(rcts_gra)
    sig_dct = automol.graph.base.sigma_radical_atom_bond_keys(rcts_gra)
    apar_dct = dict_.filter_by_value(
        automol.graph.base.atom_stereo_parities(tsg), lambda x: x is not None
    )
    bpar_dct = dict_.filter_by_value(
        automol.graph.base.bond_stereo_parities(tsg), lambda x: x is not None
    )

    xyzs = coordinates(geo)

    def _direction_vector(key):
        xyz = xyzs[key]
        # If this is a vinyl radical, the reacting electron is pointed in-plane at a 120
        # degree angle to the double bond
        if key in vin_dct:
            nkeys = automol.graph.base.atom_neighbor_atom_keys(rcts_gra, key)
            if len(nkeys) == 2:
                nxyz1, nxyz2 = map(xyzs.__getitem__, nkeys)
                rvec = vec.unit_bisector(nxyz1, nxyz2, orig_xyz=xyz, outer=True)
            elif len(nkeys) == 1:
                # Get a unit vector perpendicular to the vinyl plane
                vin_bkey = vin_dct[key]
                vin_nbh = automol.graph.base.bond_neighborhood(rcts_gra, vin_bkey)
                vin_pkeys = automol.graph.base.atom_keys(vin_nbh)
                vin_pxyzs = coordinates(geo, idxs=vin_pkeys)
                vin_nvec = vec.best_unit_perpendicular(xyzs=vin_pxyzs)
                # Get a unit vector along the vinyl bond
                nxyz = xyzs[next(iter(nkeys))]
                bvec = vec.unit_norm(numpy.subtract(nxyz, xyz))
                # Rotate the unit bond vector by 120 degrees in the vinyl plane
                rot_ = vec.rotator(vin_nvec, 2.0 * numpy.pi / 3.0)
                rvec = rot_(bvec)
                # Adjust direction to match stereochemistry
                bkey = next(
                    (bk for bk in bpar_dct if key in bk and not frm_key == bk), None
                )
                if bkey is not None:
                    # Check what the parity would be if the other atom were in this
                    # position
                    (key_,) = frm_key - {key}
                    xyz_ = numpy.add(xyz, rvec)
                    geo_ = set_coordinates(geo, {key_: xyz_})
                    par = automol.graph.base.geometry_bond_parity(tsg, geo_, bkey)
                    # If it isn't right, rotate by another 2 pi / 3
                    if par != bpar_dct[bkey]:
                        rvec = rot_(rvec)
        # If this is a sigma radical, the reacting electron is pointed outward along the
        # triple bond
        elif key in sig_dct:
            (nkey,) = sig_dct[key] - {key}
            nxyz = xyzs[nkey]
            rvec = vec.unit_norm(numpy.subtract(xyz, nxyz))
        # Otherwise, assume this is a reacting pi-electron, which is pointed
        # perpendicular to the plane of the neighboring atoms
        else:
            vin_nbh = automol.graph.base.atom_neighborhood(no_rxn_bnd_gra, key)
            # Use this to find the
            vin_pkeys = automol.graph.base.atom_keys(vin_nbh)
            vin_pxyzs = coordinates(geo, idxs=vin_pkeys)
            rvec = vec.best_unit_perpendicular(xyzs=vin_pxyzs)

            # Adjust direction to match stereochemistry
            if key in apar_dct:
                # Check what the parity would be if the other atom were in this position
                (key_,) = frm_key - {key}
                xyz_ = numpy.add(xyz, rvec)
                geo_ = set_coordinates(geo, {key_: xyz_})
                par = automol.graph.base.geometry_atom_parity(tsg, geo_, key)
                # If it isn't right, negate the direction vector
                if par != apar_dct[key]:
                    rvec = numpy.negative(rvec)

        return tuple(map(float, rvec))

    rvec1 = _direction_vector(key1)
    rvec2 = _direction_vector(key2)

    return rvec1, rvec2


def distances(
    geos,
    tsg,
    geo_idx_dct=None,
    fdist_factor=1.1,
    bdist_factor=0.9,
    angles=True,
    angstrom=True,
):
    """return a dictionary of distances for certain atoms in a sequence of
    reactant geometries, shifting the keys as needed

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param angles: include distances between the ends of bond angles?
    :type angles: bool
    :param angstrom: return the distances in angstroms?
    :type angstrom: bool
    """
    dist_dct = bond_distances(
        geos,
        tsg,
        geo_idx_dct=geo_idx_dct,
        fdist_factor=fdist_factor,
        bdist_factor=bdist_factor,
        angstrom=angstrom,
    )
    if angles:
        dist_dct.update(
            angle_distances(geos, tsg, geo_idx_dct=geo_idx_dct, angstrom=angstrom)
        )
    return dist_dct


def bond_distances(
    geos, tsg, geo_idx_dct=None, fdist_factor=1.1, bdist_factor=0.9, angstrom=True
):
    """return a dictionary of bond distances for a sequence of reactant
    geometries, shifting the keys as needed

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param angstrom: return the distances in angstroms?
    :type angstrom: bool
    """
    akeys = sorted(automol.graph.base.atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(geos, ())
    tsg = automol.graph.base.relabel(tsg, geo_idx_dct)
    gra = automol.graph.base.ts.reactants_graph(tsg)
    keys = automol.graph.base.bond_keys(gra)
    dists = [distance(geo, *k, angstrom=angstrom) for k in keys]
    dist_dct = dict(zip(keys, dists))
    for key in automol.graph.base.ts.reacting_bond_keys(tsg):
        dist_dct[key] = automol.graph.base.ts.heuristic_bond_distance(
            tsg,
            *key,
            fdist_factor=fdist_factor,
            bdist_factor=bdist_factor,
            angstrom=angstrom,
        )
    return dist_dct


def angle_distances(geos, tsg, geo_idx_dct=None, angstrom=True):
    """return a dictionary of distances between ends of a bond angle for a
    sequence of reactant geometries, shifting the keys as needed

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param angstrom: return the distances in angstroms?
    """
    akeys = sorted(automol.graph.base.atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(geos, ())
    gra = automol.graph.base.ts.reactants_graph(tsg)
    gra = automol.graph.base.relabel(gra, geo_idx_dct)
    keys = [frozenset({k1, k3}) for k1, _, k3 in automol.graph.base.angle_keys(gra)]
    dists = [distance(geo, *k, angstrom=angstrom) for k in keys]
    dist_dct = dict(zip(keys, dists))
    return dist_dct


def join(
    geo1,
    geo2,
    key2,
    key3,
    r23,
    a123=85.0,
    a234=85.0,
    d1234=85.0,
    key1=None,
    key4=None,
    angstrom=True,
    degree=True,
):
    """join two geometries based on four of their atoms, two on the first
    and two on the second

    Variables set the coordinates for 1-2...3-4 where 1-2 are bonded atoms in
    geo1 and 3-4 are bonded atoms in geo2.
    """
    key3 = key3 - count(geo1)
    a123 *= phycon.DEG2RAD if degree else 1
    a234 *= phycon.DEG2RAD if degree else 1
    d1234 *= phycon.DEG2RAD if degree else 1

    gra1, gra2 = map(connectivity_graph, (geo1, geo2))
    key1 = (
        automol.graph.base.atom_neighbor_atom_key(gra1, key2) if key1 is None else key1
    )
    key4 = (
        automol.graph.base.atom_neighbor_atom_key(gra2, key3) if key4 is None else key4
    )

    syms1 = symbols(geo1)
    syms2 = symbols(geo2)
    xyzs1 = coordinates(geo1, angstrom=angstrom)
    xyzs2 = coordinates(geo2, angstrom=angstrom)

    xyz1 = xyzs1[key1] if key1 is not None else None
    xyz2 = xyzs1[key2]
    orig_xyz3 = xyzs2[key3]
    orig_xyz4 = xyzs2[key4] if key4 is not None else [1.0, 1.0, 1.0]

    if key1 is None:
        # If the first fragment is monatomic, we can place the other one
        # anywhere we want (direction doesn't matter)
        xyz3 = numpy.add(xyz2, [0.0, 0.0, r23])
    else:
        # If the first fragment isn't monatomic, we need to take some care in
        # where we place the second one.
        # r23 and a123 are fixed, so the only degree of freedom we have to do
        # this is to optimize a dihedral angle relative to an arbitrary point
        # xyz0.
        # This dihedral angle is optimized to maximize the distance of xyz3
        # from all of the atoms in fragment 1 (xyzs1).
        xyz1 = xyzs1[key1]

        # Place xyz3 as far away from the atoms in geo1 as possible by
        # optimizing the undetermined dihedral angle
        xyz0 = vec.arbitrary_unit_perpendicular(xyz2, orig_xyz=xyz1)

        def _distance_norm(dih):  # objective function for minimization
            (dih,) = dih
            xyz3 = vec.from_internals(
                dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1, dih=dih, xyz3=xyz0
            )
            dist_norm = numpy.linalg.norm(numpy.subtract(xyzs1, xyz3))
            # return the negative norm so that minimum value gives maximum
            # distance
            return -dist_norm

        res = scipy.optimize.basinhopping(_distance_norm, 0.0)
        dih = res.x[0]

        # Now, get the next position with the optimized dihedral angle
        xyz3 = vec.from_internals(
            dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1, dih=dih, xyz3=xyz0
        )

    r34 = vec.distance(orig_xyz3, orig_xyz4)

    # If 2 doen't have neighbors or 1-2-3 are linear, ignore the dihedral angle
    if key1 is None or numpy.abs(a123 * phycon.RAD2DEG - 180.0) < 5.0:
        xyz4 = vec.from_internals(dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2)
    else:
        xyz4 = vec.from_internals(
            dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2, dih=d1234, xyz3=xyz1
        )

    align_ = vec.aligner(orig_xyz3, orig_xyz4, xyz3, xyz4)
    xyzs2 = tuple(map(align_, xyzs2))

    syms = syms1 + syms2
    xyzs = xyzs1 + xyzs2

    geo = from_data(syms, xyzs, angstrom=angstrom)
    return geo
