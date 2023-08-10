""" generate ts geometries
"""
import numpy
import scipy
from phydat import phycon

import automol.graph.base
from automol.extern import py3dmol_
from automol.geom._conv import graph_without_stereo, py3dmol_view
from automol.geom.base import (
    coordinates,
    count,
    distance,
    from_data,
    rotate,
    symbols,
    translate,
)
from automol.util import vec


def geometry_from_reactants(
    geos,
    tsg,
    geo_idx_dct=None,
    fdist_factor=1.1,
    bdist_factor=0.9,
    max_dist_err=2e-1,
    debug_visualize=False,
    log=False,
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
    :param max_dist_err: The distance convergence threshold, in angstroms
    :type max_dist_err: float, optional
    :param debug_visualize: Prompt visualizations in a Jupyter notebook for debugging?
    :type debug_visualize: bool, optional
    :param log: Print optimization log?, defaults to False
    :type log: bool, optional
    :return: TS geometry
    :rtype: automol geom data structure
    """
    keys = sorted(automol.graph.base.atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )
    tsg = automol.graph.base.relabel(tsg, geo_idx_dct)

    # 1. Join geometries for bimolecular reactions, yielding a single starting structure
    if len(geos) > 1:
        assert (
            len(geos) == 2 and len(automol.graph.base.ts.forming_bond_keys(tsg)) == 1
        ), "Generating a TS geometry for this case is not implemented:\n{tsg}"
        geo = join_at_forming_bond(
            geos,
            tsg,
            fdist_factor=fdist_factor,
            debug_visualize=debug_visualize,
        )
    else:
        (geo,) = geos

    # 2. Correct the stereochemistry against the TS graph, so it is consistent with both
    # reactants and products
    ts_geo = automol.graph.stereo_corrected_geometry(tsg, geo, geo_idx_dct=geo_idx_dct)

    # 3. Embed the TS structure, using distances from the *original* reactant geometries
    # along with forming/breaking bond distances
    dist_dct = distances(
        geos,
        tsg,
        geo_idx_dct=geo_idx_dct,
        fdist_factor=fdist_factor,
        bdist_factor=bdist_factor,
    )
    dist_range_dct = automol.graph.embed.distance_ranges_from_coordinates(
        tsg, dist_dct, angstrom=True, degree=True
    )
    ts_geo = automol.graph.embed.clean_geometry(
        tsg,
        ts_geo,
        geos=geos,
        dist_range_dct=dist_range_dct,
        relax_angles=automol.graph.base.ts.has_reacting_ring(tsg),
        max_dist_err=max_dist_err,
        log=log,
    )
    return ts_geo


def join_at_forming_bond(
    geos, tsg, geo_idx_dct=None, fdist_factor=1.1, debug_visualize=False
):
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
    :param debug_visualize: Prompt visualizations in a Jupyter notebook for debugging?
    :type debug_visualize: bool, optional
    :returns: The joined geometry
    :rtype: automol geom data structure
    """
    rcts_gra = automol.graph.ts.reactants_graph(tsg, stereo=False)
    geos_gra = automol.graph.base.union_from_sequence(
        list(map(automol.geom.graph_without_stereo, geos)), shift_keys=True
    )
    assert geos_gra == rcts_gra, (
        f"The geometries don't match the TS graph. Have they been reordered?"
        f"\ngeos:\n{geos}\ntsg:\n{tsg}"
    )
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
    rvec1 = reacting_electron_direction(geo, tsg, fidx1)
    rvec2 = reacting_electron_direction(geo, tsg, fidx2)
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

    if debug_visualize:
        view = py3dmol_.create_view()
        # Visualize the geomtry
        rcts_gra = automol.graph.ts.reactants_graph(tsg)
        view = py3dmol_view(geo, gra=rcts_gra, view=view)
        # Visualize the first direction vector
        view = py3dmol_.view_vector(rvec1, orig_xyz=fxyz1, color="blue", view=view)
        # Visualize the rotational axis
        fxyz2 = coordinates(geo)[fidx2]
        view = py3dmol_.view_vector(rot_axis, orig_xyz=fxyz2, color="red", view=view)
        # Visualize the second direction vector
        rot_ = vec.rotator(rot_axis, rot_ang)
        rvec2 = rot_(rvec2)
        view = py3dmol_.view_vector(rvec2, orig_xyz=fxyz2, color="green", view=view)
        view.show()

    return geo


def reacting_electron_direction(geo, tsg, key) -> vec.Vector:
    """Identify the direction of a reacting electron on a bond-forming atom

    Does not correct for stereochemistry.

    :param geo: Reactants geometry, aligned to the TS graph
    :type geo: automol geom data structure
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key: Key for the atom, which must be part of a forming bond
    :type key: int
    :returns: A vector indicating the direction
    :rtype: vec.Vector
    """
    frm_key = next(
        (k for k in automol.graph.base.ts.forming_bond_keys(tsg) if key in k), None
    )
    assert frm_key is not None, f"Atom {key} is not forming a bond in this graph:{tsg}"

    # Get the normal vector
    pkeys = automol.graph.base.ts.plane_keys(tsg, key)
    pxyzs = coordinates(geo, idxs=pkeys)
    zvec = vec.best_unit_perpendicular(pxyzs)

    # Get the direction information:
    #   1. xkey: a bond key giving an x direction
    #   2. ykey: a bond key giving an y direction
    #   3. phi: a rotational angle
    # The electron direction is obtained by rotating the x direction by phi around a z
    # axis of a right-handed coordinate system (happens in the `else` below)
    xkey, ykey, phi = automol.graph.base.ts.reacting_electron_direction(tsg, key)
    # `None` indicates that the electron is perpendicular to the plane, along the normal
    # vector
    if xkey is None:
        rvec = zvec
    # Otherwise, the electron is in-plane, given by rotating an existing bond direction
    # by `phi`
    else:
        xxyz1, xxyz2 = coordinates(geo, idxs=xkey)
        xvec = numpy.subtract(xxyz2, xxyz1)
        if ykey is not None:
            yxyz1, yxyz2 = coordinates(geo, idxs=ykey)
            yvec = numpy.subtract(yxyz2, yxyz1)
            zvec = vec.flip_if_left_handed(xvec, yvec, zvec)
        xvec = vec.orthogonalize(zvec, xvec, normalize=True)
        rot_ = vec.rotator(zvec, phi)
        rvec = rot_(xvec)

    return rvec


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
    keys = sorted(automol.graph.base.atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(geos, ())
    tsg = automol.graph.base.relabel(tsg, geo_idx_dct)
    gra = automol.graph.base.ts.reactants_graph(tsg)

    # Add measured bond distances from the reactants geometry
    bnd_keys = automol.graph.base.bond_keys(gra)
    dists = [distance(geo, *k, angstrom=angstrom) for k in bnd_keys]
    dist_dct = dict(zip(bnd_keys, dists))
    # Add measured angle distances from the reactants geometry, if requested
    if angles:
        ang_keys = [
            frozenset({k1, k3}) for k1, _, k3 in automol.graph.base.angle_keys(gra)
        ]
        dists = [distance(geo, *k, angstrom=angstrom) for k in ang_keys]
        dist_dct = dict(zip(ang_keys, dists))

    # Add heuristic bond distances from the TS graph
    for key in automol.graph.base.ts.reacting_bond_keys(tsg):
        dist_dct[key] = automol.graph.base.ts.heuristic_bond_distance(
            tsg,
            *key,
            fdist_factor=fdist_factor,
            bdist_factor=bdist_factor,
            angstrom=angstrom,
        )
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

    gra1, gra2 = map(graph_without_stereo, (geo1, geo2))
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
