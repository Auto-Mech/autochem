"""
  Functions used for handling and comparing multiple geometries
"""
import functools
import itertools
from typing import List

import numpy

from phydat import ptab

from ._0core import (
    coordinates,
    distance,
    from_string,
    symbols,
    xyz_string,
)


# # properties used for comparisons
def coulomb_spectrum(geo):
    """Calculate a Coulomb matrix eigenvalue spectrum where
    the eignevalues are sorted in ascending order.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :rtype: tuple(float)
    """

    mat = _coulomb_matrix(geo)
    vals = tuple(sorted(numpy.linalg.eigvalsh(mat)))
    return vals


def _coulomb_matrix(geo):
    """Calculate the Coulomb matrix wich describes the
    electrostatic interactions between nuclei:

    M[i,j] = 0.5Z_i^2.4 (i=j), Z_iZ_j/R_ij (i!=j

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :rtype: tuple(tuple(float))
    """

    nums = numpy.array(list(map(ptab.to_number, symbols(geo))))
    xyzs = numpy.array(coordinates(geo))

    _ = numpy.newaxis
    natms = len(nums)
    diag_idxs = numpy.diag_indices(natms)
    tril_idxs = numpy.tril_indices(natms, -1)
    triu_idxs = numpy.triu_indices(natms, 1)

    zxz = numpy.outer(nums, nums)
    rmr = numpy.linalg.norm(xyzs[:, _, :] - xyzs[_, :, :], axis=2)

    mat = numpy.zeros((natms, natms))
    mat[diag_idxs] = nums**2.4 / 2.0
    mat[tril_idxs] = zxz[tril_idxs] / rmr[tril_idxs]
    mat[triu_idxs] = zxz[triu_idxs] / rmr[triu_idxs]

    return mat


def distance_matrix(geo):
    """Form a Natom X Natom matrix containing the distance of all the
    atoms in a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: numpy.ndarray
    """

    mat = numpy.zeros((len(geo), len(geo)))
    for i in range(len(geo)):
        for j in range(len(geo)):
            mat[i][j] = distance(geo, i, j)

    return mat


# # comparisons
def almost_equal(geo1, geo2, rtol=2e-3, atol=2e-6):
    """Assess if the coordinates of two molecular geometries
    are numerically equal.

    :param geo1: molecular geometry 1
    :type geo1: automol molecular geometry data structure
    :param geo2: molecular geometry 2
    :type geo2: automol molecular geometry data structure
    :param rtol: Relative tolerance for the distances
    :type rtol: float
    :rtype: tuple(tuple(float))
    """

    ret = False
    if symbols(geo1) == symbols(geo2):
        ret = numpy.allclose(coordinates(geo1), coordinates(geo2), rtol=rtol, atol=atol)

    return ret


def almost_equal_coulomb_spectrum(geo1, geo2, rtol=1e-2):
    """Assess if two molecular geometries have similar coulomb spectrums.

    :param geo1: molecular geometry 1
    :type geo1: automol molecular geometry data structure
    :param geo2: molecular geometry 2
    :type geo2: automol molecular geometry data structure
    :param rtol: Relative tolerance for the distances
    :type rtol: float
    :rtype: bool
    """
    return numpy.allclose(coulomb_spectrum(geo1), coulomb_spectrum(geo2), rtol=rtol)


def argunique_coulomb_spectrum(geos, seen_geos=(), rtol=1e-2):
    """Get indices of unique geometries, by coulomb spectrum.

    :param geos: list of molecular geometries
    :type geos: tuple(automol molecular geometry data structure)
    :param seen_geos: geometries that have been assessed
    :type seen_geos: tuple(automol molecular geometry data structure)
    :param rtol: Relative tolerance for the distances
    :type rtol: float
    :rtype: tuple(int)
    """
    comp_ = functools.partial(almost_equal_coulomb_spectrum, rtol=rtol)
    idxs = _argunique(geos, comp_, seen_items=seen_geos)
    return idxs


def _argunique(items, comparison, seen_items=()):
    """Get the indices of unique items using some comparison function.

    :param items: items to assess for uniqueness
    :type items: tuple(obj)
    :param comparison: function used to compare items
    :type comparison: function object
    :param seen_items: items that have been assessed
    :type seen_items: tuple(obj)
    :rtype: tuple(int)
    """

    idxs = []
    seen_items = list(seen_items)
    for idx, item in enumerate(items):
        if not any(comparison(item, seen_item) for seen_item in seen_items):
            idxs.append(idx)
            seen_items.append(item)
    idxs = tuple(idxs)

    return idxs


def almost_equal_dist_matrix(geo1, geo2, thresh=0.1, idxs: List[int] = None):
    """form distance matrix for a set of xyz coordinates

    :param geo1: A geometry
    :type geo1: automol geom data structure
    :param geo2: Another geometry
    :type geo2: automol geom data structure
    :param idxs: Optionally, restrict this to a subset of indices
    :type idxs: Optional[List[int]]
    """
    idxs = list(range(len(geo1))) if idxs is None else idxs
    for i, j in itertools.combinations(idxs, 2):
        dist1 = distance(geo1, i, j)
        dist2 = distance(geo2, i, j)
        if abs(dist1 - dist2) > thresh:
            return False

    return True


def minimum_volume_geometry(geos):
    """Generate the geometry with smallest volume from a set
    of geometrties for a given species.

    :param geos: molecular geometries
    :type geos: tuple(geo obj)
    :rtype: geo obj
    """

    # Get lines for the output string
    geoms_string = "\n".join([xyz_string(geom) for geom in geos])
    lines = geoms_string.splitlines()

    # Get the number of atoms
    natom = int(lines[0])

    # loop over the lines to find the smallest geometry
    rrminmax = 1.0e10
    ngeom = 0
    small_geo_idx = 0
    while ngeom * (natom + 2) < len(lines):
        rrmax = 0.0
        for i in range(natom):
            for j in range(i + 1, natom):
                # Get the line
                xyz1 = lines[i + ngeom * (natom + 2) + 2].split()[1:]
                xyz2 = lines[j + ngeom * (natom + 2) + 2].split()[1:]
                # Get the coordinates
                atom1 = [float(val) for val in xyz1]
                atom2 = [float(val) for val in xyz2]
                # Calculate the interatomic distance
                rrtest = numpy.sqrt(
                    (atom1[0] - atom2[0]) ** 2
                    + (atom1[1] - atom2[1]) ** 2
                    + (atom1[2] - atom2[2]) ** 2
                )
                # Check and see if distance is more than max
                rrmax = max(rrmax, rrtest)
        # If max below moving threshold, set to smallest geom
        if rrmax < rrminmax:
            rrminmax = rrmax
            small_geo_idx = ngeom
        ngeom += 1

    # Set the output geometry
    geom_str = ""
    for i in range(natom):
        geom_str += lines[i + small_geo_idx * (natom + 2) + 2] + "\n"
    round_geom = from_string(geom_str)

    return round_geom
