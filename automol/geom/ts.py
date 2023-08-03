""" generate ts geometries
"""

import numpy
import qcelemental
import scipy
from phydat import phycon

import automol.graph.base
from automol.geom._conv import connectivity_graph
from automol.geom.base import coordinates, count, distance, from_data, symbols
from automol.util import vec


def geometry_from_reactants(
    geos, ts_gra, geo_idx_dct=None, fdist_factor=1.1, bdist_factor=0.9
):
    """Generate a TS geometry from reactants

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :return: TS geometry
    :rtype: automol geom data structure
    """
    # 1. Evaluate the distances, updating them for the TS graph
    dist_dct = distances(
        geos,
        ts_gra,
        geo_idx_dct=geo_idx_dct,
        fdist_factor=fdist_factor,
        bdist_factor=bdist_factor,
    )
    print("distances:")
    print(dist_dct)

    # 2. Join geometries together, if bimolecular
    if len(geos) > 1:
        raise NotImplementedError

    (geo,) = geos

    ts_geo = automol.graph.stereo_corrected_geometry(
        ts_gra, geo, geo_idx_dct=geo_idx_dct
    )
    return ts_geo


def distances(
    geos,
    ts_gra,
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
    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
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
    dist_dct = bond_distances(geos, ts_gra, geo_idx_dct=geo_idx_dct, angstrom=angstrom)
    if angles:
        dist_dct.update(
            angle_distances(geos, ts_gra, geo_idx_dct=geo_idx_dct, angstrom=angstrom)
        )
    if ts_gra is not None:
        dist_dct.update(
            reacting_bond_distances(
                ts_gra,
                geo_idx_dct=geo_idx_dct,
                fdist_factor=fdist_factor,
                bdist_factor=bdist_factor,
                angstrom=angstrom,
            )
        )
    return dist_dct


def bond_distances(geos, ts_gra, geo_idx_dct=None, angstrom=True):
    """return a dictionary of bond distances for a sequence of reactant
    geometries, shifting the keys as needed

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param angstrom: return the distances in angstroms?
    """
    akeys = sorted(automol.graph.base.atom_keys(ts_gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(geos, ())
    gra = automol.graph.base.ts.reactants_graph(ts_gra)
    gra = automol.graph.base.relabel(gra, geo_idx_dct)
    keys = automol.graph.base.bond_keys(gra)
    dists = [distance(geo, *k, angstrom=angstrom) for k in keys]
    dist_dct = dict(zip(keys, dists))
    return dist_dct


def angle_distances(geos, ts_gra, geo_idx_dct=None, angstrom=True):
    """return a dictionary of distances between ends of a bond angle for a
    sequence of reactant geometries, shifting the keys as needed

    :param geos: Reactant geometries
    :type geos: List[automol geom data structure]
    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param angstrom: return the distances in angstroms?
    """
    akeys = sorted(automol.graph.base.atom_keys(ts_gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(geos, ())
    gra = automol.graph.base.ts.reactants_graph(ts_gra)
    gra = automol.graph.base.relabel(gra, geo_idx_dct)
    keys = [frozenset({k1, k3}) for k1, _, k3 in automol.graph.base.angle_keys(gra)]
    dists = [distance(geo, *k, angstrom=angstrom) for k in keys]
    dist_dct = dict(zip(keys, dists))
    return dist_dct


def reacting_bond_distances(
    ts_gra,
    geo_idx_dct=None,
    fdist_factor=1.1,
    bdist_factor=0.9,
    angstrom=True,
):
    """Get a dictionary of reacting bond distances for a TS geometry, based on a TS graph

    :param ts_gra: TS graph
    :type ts_gra: automol graph data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param angstrom: return the distances in angstroms?
    :type angstrom: True
    """
    units = "angstrom" if angstrom else "bohr"
    akeys = sorted(automol.graph.base.atom_keys(ts_gra))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )
    frm_bkeys = automol.graph.base.ts.forming_bond_keys(ts_gra)

    symb_dct = automol.graph.base.atom_symbols(ts_gra)
    dist_dct = {}
    for bkey in automol.graph.base.ts.reacting_bond_keys(ts_gra):
        # Calculate the average van der Waals radius
        dist = sum(qcelemental.vdwradii.get(symb_dct[k], units=units) for k in bkey) / 2

        # Multiply by the chosen factors for forming/breaking bonds
        if bkey in frm_bkeys:
            dist *= fdist_factor
        else:
            dist *= bdist_factor

        # Get the geometry indices for this bond key
        bidxs = frozenset(map(geo_idx_dct.__getitem__, bkey))
        dist_dct[bidxs] = dist

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
