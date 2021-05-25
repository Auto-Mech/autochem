""" independent geom tools for graph
"""
import itertools
import numpy
from phydat import phycon
from phydat import ptab
import automol.util.vec
import automol.create as _create


def coordinates(geo, idxs=None, angstrom=False):
    """ Obtain the Cartesian coordinates of atoms in the molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indexs of atoms to obtain information for
        :type idxs: tuple(int)
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :rtype: tuple(tuple(float))
    """

    idxs = list(range(count(geo))) if idxs is None else idxs
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    xyzs = xyzs if not angstrom else numpy.multiply(xyzs, phycon.BOHR2ANG)
    xyzs = tuple(tuple(xyz) for idx, xyz in enumerate(xyzs) if idx in idxs)

    return xyzs


def count(geo):
    """ Obtain the number of rows of the molecular geometry, which corresponds to
        the number of atoms in the geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: int
    """
    return len(geo)


def geometry_join(
        geo1, geo2, dist_cutoff=3.0*phycon.ANG2BOHR, theta=0.0, phi=0.0):
    """ Join two molecular geometries together where the intermolecular
        separation and orientation can be specified.

        :param geo1: molecular geometry 1
        :type geo1: automol molecular geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol molecular geometry data structure
        :param dist_cutoff: threshhold for center-of-mass distance
        :type: dist_cutoff: float
        :param theta: theta angle for intermolecular orientation
        :type theta: float
        :param phi: phi angle for intermolecular orientation
        :type phi: float
        :rtype: automol molecular geometry data structure
    """

    if not geo1:
        symbs = symbols(geo2)
        xyzs = coordinates(geo2)
    elif not geo2:
        symbs = symbols(geo1)
        xyzs = coordinates(geo1)
    else:
        orient_vec = numpy.array([numpy.sin(theta) * numpy.cos(phi),
                                  numpy.sin(theta) * numpy.sin(phi),
                                  numpy.cos(theta)])
        neg_orient_vec = -1.0 * orient_vec

        # get the correct distance apart
        geo1 = mass_centered(geo1)
        geo2 = mass_centered(geo2)
        ext1 = max(numpy.vdot(orient_vec, xyz)
                   for xyz in coordinates(geo1))
        ext2 = max(numpy.vdot(neg_orient_vec, xyz)
                   for xyz in coordinates(geo2))

        cm_dist = ext1 + dist_cutoff + ext2
        dist_grid = numpy.arange(cm_dist, 0., -0.1)
        for dist in dist_grid:
            trans_geo2 = translate(geo2, orient_vec * dist)
            min_dist = minimum_distance(geo1, trans_geo2)
            if numpy.abs(min_dist - dist_cutoff) < 0.1:
                break

        geo2 = trans_geo2

        # now, join them together
        symbs = symbols(geo1) + symbols(geo2)
        xyzs = coordinates(geo1) + coordinates(geo2)

    return _create.geom.from_data(symbs, xyzs)


def mass_centered(geo):
    """ Generate a new geometry where the coordinates of the input geometry
        have been translated to the center-of-mass.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: tuple(float)
    """
    return translate(geo, numpy.negative(center_of_mass(geo)))


def translate(geo, xyz):
    """ Translate the coordinates of a molecular geometry along
        a three-dimensiona vector.
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param xyz: vector to translate along
        :type xyz: tuple(float)
        :rtype: automol molecular geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.add(xyzs, xyz)
    return _create.geom.from_data(symbs, xyzs)


def symbols(geo, idxs=None):
    """ Obtain the atomic symbols atoms in the molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indexs of atoms to obtain information for
        :type idxs: tuple(int)
        :rtype: tuple(str)
    """

    idxs = list(range(count(geo))) if idxs is None else idxs

    if geo:
        symbs, _ = zip(*geo)
    else:
        symbs = ()

    symbs = tuple(symb for idx, symb in enumerate(symbs) if idx in idxs)
    return symbs


def minimum_distance(geo1, geo2):
    """ get the minimum distance between atoms in geo1 and those in geo2

        :param geo1: molecular geometry 1
        :type geo1: automol molecular geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol molecular geometry data structure
        :rtype: float
    """

    xyzs1 = coordinates(geo1)
    xyzs2 = coordinates(geo2)
    return min(automol.util.vec.distance(xyz1, xyz2)
               for xyz1, xyz2 in itertools.product(xyzs1, xyzs2))


def center_of_mass(geo):
    """ Determine the center-of-mass for a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: tuple(float)
    """

    xyzs = coordinates(geo)
    amas = masses(geo)
    cm_xyz = tuple(
        sum(numpy.multiply(xyz, ama) for xyz, ama in zip(xyzs, amas)) /
        sum(amas))

    return cm_xyz


def masses(geo, amu=True):
    """ Build a list of the atomic masses that corresponds to the list
        of atomic sybmols of a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(float)
    """

    symbs = symbols(geo)
    amas = list(map(ptab.to_mass, symbs))

    if not amu:
        amas = numpy.multiply(amas, phycon.AMU2EMASS)

    amas = tuple(amas)

    return amas


def is_atom(geo):
    """ Determine if the molecular geometry corresponds to an atom.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: bool
    """

    return len(symbols(geo)) == 1
