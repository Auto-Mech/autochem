""" generate ts geometries
"""

import numpy
import scipy
from phydat import phycon
from automol.util import vec
import automol.graph
from automol.geom._conv import connectivity_graph
from automol.geom.base import from_data
from automol.geom.base import symbols
from automol.geom.base import coordinates
from automol.geom.base import count
from automol.geom.base import distance


def distances(geos, bonds=True, angles=True, angstrom=True):
    """ return a dictionary of distances for certain atoms in a sequence of
    reactant geometries, shifting the keys as needed

    :param geos: the reactant geometries
    :param bonds: include bond distances?
    :param angles: include distances between the ends of bond angles?
    :param angstrom: return the distances in angstroms?
    """
    dist_dct = {}
    if bonds:
        dist_dct.update(bond_distances(geos, angstrom=angstrom))
    if angles:
        dist_dct.update(angle_distances(geos, angstrom=angstrom))
    return dist_dct


def bond_distances(geos, angstrom=True):
    """ return a dictionary of bond distances for a sequence of reactant
    geometries, shifting the keys as needed

    :param geos: the reactant geometries
    :param angstrom: return the distances in angstroms?
    """
    dist_dct = {}

    shift = 0
    for geo in geos:
        gra = connectivity_graph(geo)
        pairs = list(automol.graph.bond_keys(gra))
        keys = [frozenset({k1+shift, k2+shift}) for k1, k2 in pairs]
        dists = [distance(geo, *p, angstrom=angstrom) for p in pairs]
        dist_dct.update(dict(zip(keys, dists)))

        shift += count(geo)

    return dist_dct


def angle_distances(geos, angstrom=True):
    """ return a dictionary of distances between ends of a bond angle for a
    sequence of reactant geometries, shifting the keys as needed

    :param geos: the reactant geometries
    :param angstrom: return the distances in angstroms?
    """
    dist_dct = {}

    shift = 0
    for geo in geos:
        gra = connectivity_graph(geo)
        pairs = [(k1, k3) for k1, _, k3 in automol.graph.angle_keys(gra)]
        keys = [frozenset({k1+shift, k2+shift}) for k1, k2 in pairs]
        dists = [distance(geo, *p, angstrom=angstrom) for p in pairs]
        dist_dct.update(dict(zip(keys, dists)))

        shift += count(geo)

    return dist_dct


def join(geo1, geo2, key2, key3, r23, a123=85., a234=85., d1234=85.,
         key1=None, key4=None, angstrom=True, degree=True):
    """ join two geometries based on four of their atoms, two on the first
    and two on the second

    Variables set the coordinates for 1-2...3-4 where 1-2 are bonded atoms in
    geo1 and 3-4 are bonded atoms in geo2.
    """
    key3 = key3 - count(geo1)
    a123 *= phycon.DEG2RAD if degree else 1
    a234 *= phycon.DEG2RAD if degree else 1
    d1234 *= phycon.DEG2RAD if degree else 1

    gra1, gra2 = map(connectivity_graph, (geo1, geo2))
    key1 = (automol.graph.atom_neighbor_atom_key(gra1, key2) if key1 is None
            else key1)
    key4 = (automol.graph.atom_neighbor_atom_key(gra2, key3) if key4 is None
            else key4)

    syms1 = symbols(geo1)
    syms2 = symbols(geo2)
    xyzs1 = coordinates(geo1, angstrom=angstrom)
    xyzs2 = coordinates(geo2, angstrom=angstrom)

    xyz1 = xyzs1[key1] if key1 is not None else None
    xyz2 = xyzs1[key2]
    orig_xyz3 = xyzs2[key3]
    orig_xyz4 = xyzs2[key4] if key4 is not None else [1., 1., 1.]

    if key1 is None:
        # If the first fragment is monatomic, we can place the other one
        # anywhere we want (direction doesn't matter)
        xyz3 = numpy.add(xyz2, [0., 0., r23])
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
            dih, = dih
            xyz3 = vec.from_internals(dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1,
                                      dih=dih, xyz3=xyz0)
            dist_norm = numpy.linalg.norm(numpy.subtract(xyzs1, xyz3))
            # return the negative norm so that minimum value gives maximum
            # distance
            return -dist_norm

        res = scipy.optimize.basinhopping(_distance_norm, 0.)
        dih = res.x[0]

        # Now, get the next position with the optimized dihedral angle
        xyz3 = vec.from_internals(dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1,
                                  dih=dih, xyz3=xyz0)

    r34 = vec.distance(orig_xyz3, orig_xyz4)

    # If 2 doen't have neighbors or 1-2-3 are linear, ignore the dihedral angle
    if key1 is None or numpy.abs(a123 * phycon.RAD2DEG - 180.) < 5.:
        xyz4 = vec.from_internals(dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2)
    else:
        xyz4 = vec.from_internals(dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2,
                                  dih=d1234, xyz3=xyz1)

    align_ = vec.aligner(orig_xyz3, orig_xyz4, xyz3, xyz4)
    xyzs2 = tuple(map(align_, xyzs2))

    syms = syms1 + syms2
    xyzs = xyzs1 + xyzs2

    geo = from_data(syms, xyzs, angstrom=angstrom)
    return geo


def _maximize_xyz3_distance(dist, xyz1, ang, xyz2, against_xyzs=()):
    """ place another atom, maximizing its distance from another set of atoms
    """
    xyz0 = vec.arbitrary_unit_perpendicular(xyz2, orig_xyz=xyz1)

    against_xyzs = numpy.array(against_xyzs)

    def _distance_norm(dih):
        xyz3 = vec.from_internals(dist=dist, xyz1=xyz1, ang=ang, xyz2=xyz2,
                                  dih=dih, xyz3=xyz0)
        dist_norm = numpy.linalg.norm(against_xyzs - numpy.array(xyz3))
        return dist_norm
