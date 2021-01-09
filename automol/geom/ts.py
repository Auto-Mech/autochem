""" generate ts geometries
"""
import numpy
from qcelemental import constants as qcc
import automol.graph
import automol.convert.geom
from automol.util import vec
from automol.geom._geom import count
from automol.geom._geom import symbols
from automol.geom._geom import coordinates

DEG2RAD = qcc.conversion_factor('degree', 'radian')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


def join(geo1, geo2, key2, key3, r23, a123, a234, d1234=85., key1=None,
         key4=None, angstrom=True, degree=True):
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

    xyz3 = vec.from_internals(dist=r23, xyz1=xyz2, ang=a123, xyz2=xyz1)

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
