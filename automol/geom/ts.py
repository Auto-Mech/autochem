""" generate ts geometries
"""
import numpy
import scipy
from qcelemental import constants as qcc
import automol.graph
import automol.convert.geom
from automol.util import vec
from automol.geom._geom import count
from automol.geom._geom import symbols
from automol.geom._geom import distance
from automol.geom._geom import coordinates

DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


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
        gra = automol.convert.geom.connectivity_graph(geo)
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
        gra = automol.convert.geom.connectivity_graph(geo)
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
    a123 *= DEG2RAD if degree else 1
    a234 *= DEG2RAD if degree else 1
    d1234 *= DEG2RAD if degree else 1

    gra1, gra2 = map(automol.convert.geom.connectivity_graph, (geo1, geo2))
    key1 = (automol.graph.atom_neighbor_key(gra1, key2) if key1 is None
            else key1)
    key4 = (automol.graph.atom_neighbor_key(gra2, key3) if key4 is None
            else key4)

    syms1 = symbols(geo1)
    syms2 = symbols(geo2)
    xyzs1 = coordinates(geo1, angstrom=angstrom)
    xyzs2 = coordinates(geo2, angstrom=angstrom)

    xyz1, xyz2 = map(xyzs1.__getitem__, (key1, key2))
    orig_xyz3, orig_xyz4 = map(xyzs2.__getitem__, (key3, key4))

    r34 = vec.distance(orig_xyz3, orig_xyz4)

    # Place xyz3 as far away from the atoms in geo1 as possible by optimizing
    # the undetermined dihedral angle
    xyz0 = vec.arbitrary_unit_perpendicular(xyz2, orig_xyz=xyz1)

    def _distance_norm(dih):  # objective function for minimization
        xyz3 = vec.from_internals(dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1,
                                  dih=dih, xyz3=xyz0)
        dist_norm = numpy.linalg.norm(numpy.subtract(xyzs1, xyz3))
        # return the negative norm so that minimum value gives maximum distance
        return -dist_norm

    res = scipy.optimize.basinhopping(_distance_norm, 0.)
    dih = res.x[0]

    # Now, get the next position with the optimized dihedral angle
    xyz3 = vec.from_internals(dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1,
                              dih=dih, xyz3=xyz0)

    # Don't use the dihedral angle if 1-2-3 are linear
    if numpy.abs(a123 * RAD2DEG - 180.) > 5.:
        xyz4 = vec.from_internals(dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2,
                                  dih=d1234, xyz3=xyz1)
    else:
        xyz4 = vec.from_internals(dist=r34, xyz1=xyz3, ang=a234, xyz2=xyz2)

    align_ = vec.aligner(orig_xyz3, orig_xyz4, xyz3, xyz4)
    xyzs2 = tuple(map(align_, xyzs2))

    syms = syms1 + syms2
    xyzs = xyzs1 + xyzs2

    geo = automol.create.geom.from_data(syms, xyzs, angstrom=angstrom)
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
