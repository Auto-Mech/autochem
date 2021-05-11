"""
  Chemical properties
"""

import more_itertools as mit
import numpy
from phydat import phycon, ptab
from automol import util
from automol.geom import _base as geom_base
from automol.graph.geom import mass_centered
from automol.graph.geom import is_atom as _is_atom


def is_atom(geo):
    """ Determine if the molecular geometry corresponds to an atom.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: bool
    """

    return _is_atom(geo)


def is_linear(geo, tol=2.*phycon.DEG2RAD):
    """ Determine if the molecular geometry is linear.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param tol: tolerance of bond angle(s) for determing linearity
        :type tol: float
        :rtype: bool
    """

    ret = True

    if len(geo) == 1:
        ret = False
    elif len(geo) == 2:
        ret = True
    else:
        keys = range(len(geom_base.symbols(geo)))
        for key1, key2, key3 in mit.windowed(keys, 3):
            cangle = numpy.abs(geom_base.central_angle(geo, key1, key2, key3))
            if not (numpy.abs(cangle) < tol or
                    numpy.abs(cangle - numpy.pi) < tol):
                ret = False
    return ret


def masses(geo, amu=True):
    """ Build a list of the atomic masses that corresponds to the list
        of atomic sybmols of a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(float)
    """

    symbs = geom_base.symbols(geo)
    amas = list(map(ptab.to_mass, symbs))

    if not amu:
        amas = numpy.multiply(amas, phycon.AMU2EMASS)

    amas = tuple(amas)

    return amas


def total_mass(geo, amu=True):
    """ Calculate the total mass of a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(float)
    """
    return sum(masses(geo, amu=amu))


def center_of_mass(geo):
    """ Determine the center-of-mass for a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: tuple(float)
    """

    xyzs = geom_base.coordinates(geo)
    amas = masses(geo)
    cm_xyz = tuple(
        sum(numpy.multiply(xyz, ama) for xyz, ama in zip(xyzs, amas)) /
        sum(amas))

    return cm_xyz


def reduced_mass(geo1, geo2):
    """ Calculate the reduced mass for two species.

        :param geo1: geometry of species 1 (Bohr)
        :type geo1: list(float)
        :param geo2: geometry of species 2 (Bohr)
        :type geo2: list(float)
        :return: reduced mass (amu)
        :rtype: float
    """

    mass1 = total_mass(geo1)
    mass2 = total_mass(geo2)

    return (mass1 * mass2) / (mass1 + mass2)


def inertia_tensor(geo, amu=True):
    """ Build the moment-of-inertia tensor for a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(tuple(float))
    """

    geo = mass_centered(geo)
    amas = masses(geo, amu=amu)
    xyzs = geom_base.coordinates(geo)
    ine = tuple(map(tuple, sum(
        ama * (numpy.vdot(xyz, xyz) * numpy.eye(3) - numpy.outer(xyz, xyz))
        for ama, xyz in zip(amas, xyzs))))

    return ine


def principal_axes(geo, amu=True):
    """ Determine the principal axes of rotation for a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(tuple(float))
    """

    ine = inertia_tensor(geo, amu=amu)
    _, paxs = numpy.linalg.eigh(ine)
    paxs = tuple(map(tuple, paxs))

    return paxs


def moments_of_inertia(geo, amu=True):
    """ Calculate the moments of inertia along the xyz axes
        (these not sorted in to A,B,C).

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(tuple(float))
    """

    ine = inertia_tensor(geo, amu=amu)
    moms, _ = numpy.linalg.eigh(ine)
    moms = tuple(moms)

    return moms


def rotational_constants(geo, amu=True):
    """ Calculate the rotational constants.
        (these not sorted in to A,B,C).

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param amu: parameter to control electron mass -> amu conversion
        :type amu: bool
        :rtype: tuple(float)
    """

    moms = moments_of_inertia(geo, amu=amu)
    cons = numpy.divide(1., moms) / 4. / numpy.pi / phycon.SOL
    cons = tuple(cons)
    return cons


def permutation(geo, ref_geo, thresh=1e-4):
    """ Determine the permutation of one geometry that reproduces another
        (if there isn't one -- the geometries are not aligned, return None).

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param ref_geo: molecular geometry
        :type ref_geo: automol molecular geometry data structure
        :param thresh: theshold for assessing if permutation exists
        :type thresh: float
        :rtype: tuple(int)
    """

    natms = geom_base.count(geo)
    symbs = geom_base.symbols(geo)
    xyzs = geom_base.coordinates(geo)

    perm_idxs = [None] * natms
    for idx, (symb, xyz) in enumerate(zip(symbs, xyzs)):
        # Loop over atoms in the reference geometry with the same symbol
        ref_idxs = geom_base.atom_indices(ref_geo, symb=symb)
        ref_xyzs = geom_base.coordinates(ref_geo, idxs=ref_idxs)
        perm_idx = next(
            (ref_idx for ref_idx, ref_xyz in zip(ref_idxs, ref_xyzs)
             if util.vec.distance(xyz, ref_xyz) < thresh), None)
        perm_idxs[idx] = perm_idx

    perm_idxs = tuple(perm_idxs)

    if any(perm_idx is None for perm_idx in perm_idxs):
        perm_idxs = None

    return perm_idxs
