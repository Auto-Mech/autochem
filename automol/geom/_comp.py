"""
  Functions used for handling and comparing multiple geometries
"""

import itertools
import functools
import numpy
from phydat import ptab
from phydat import bnd
from automol.zmat._zmat import value_dictionary
from automol.zmat._new import is_atom_closest_to_bond_atom
import automol.convert.geom
from automol import util
from automol.geom import _base as geom_base


# Assessments of transition states
def _ts_compare(ref_zma, zma, zrxn):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Initialize viable
    viable = True

    # Get the bond dists and calculate the distance of bond being formed
    ref_geo, _ = automol.convert.zmat.geometry(ref_zma)
    cnf_geo, _ = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)

    frm_bnd_keys = automol.reac.forming_bond_keys(grxn)
    brk_bnd_keys = automol.reac.breaking_bond_keys(grxn)

    cnf_dist_lst = []
    ref_dist_lst = []
    bnd_key_lst = []
    cnf_ang_lst = []
    ref_ang_lst = []
    for frm_bnd_key in frm_bnd_keys:
        frm_idx1, frm_idx2 = list(frm_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, frm_idx1, frm_idx2)
        ref_dist = automol.geom.distance(ref_geo, frm_idx1, frm_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(frm_bnd_key)

    for brk_bnd_key in brk_bnd_keys:
        brk_idx1, brk_idx2 = list(brk_bnd_key)
        cnf_dist = automol.geom.distance(cnf_geo, brk_idx1, brk_idx2)
        ref_dist = automol.geom.distance(ref_geo, brk_idx1, brk_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(brk_bnd_key)

    for frm_bnd_key in frm_bnd_keys:
        for brk_bnd_key in brk_bnd_keys:
            for frm_idx in frm_bnd_key:
                for brk_idx in brk_bnd_key:
                    if frm_idx == brk_idx:
                        idx2 = frm_idx
                        idx1 = list(frm_bnd_key - frozenset({idx2}))[0]
                        idx3 = list(brk_bnd_key - frozenset({idx2}))[0]
                        cnf_ang = automol.geom.central_angle(
                            cnf_geo, idx1, idx2, idx3)
                        ref_ang = automol.geom.central_angle(
                            ref_geo, idx1, idx2, idx3)
                        cnf_ang_lst.append(cnf_ang)
                        ref_ang_lst.append(ref_ang)

    # print('bnd_key_list', bnd_key_lst)
    # print('conf_dist', cnf_dist_lst)
    # print('ref_dist', ref_dist_lst)
    # print('conf_angle', cnf_ang_lst)
    # print('ref_angle', ref_ang_lst)

    # Set the maximum allowed displacement for a TS conformer
    max_disp = 0.6
    # better to check for bond-form length in bond scission with ring forming
    if 'addition' in grxn.class_:
        max_disp = 0.8
    if 'abstraction' in grxn.class_:
        # this was 1.4 - SJK reduced it to work for some OH abstractions
        max_disp = 1.0

    # Check forming bond angle similar to ini config
    if 'elimination' not in grxn.class_:
        for ref_angle, cnf_angle in zip(ref_ang_lst, cnf_ang_lst):
            if abs(cnf_angle - ref_angle) > .44:
                print('angle', ref_angle, cnf_angle)
                viable = False

    symbols = automol.geom.symbols(cnf_geo)
    lst_info = zip(ref_dist_lst, cnf_dist_lst, bnd_key_lst)
    for ref_dist, cnf_dist, bnd_key in lst_info:
        if 'add' in grxn.class_ or 'abst' in grxn.class_:
            bnd_key1, bnd_key2 = min(list(bnd_key)), max(list(bnd_key))
            symb1 = symbols[bnd_key1]
            symb2 = symbols[bnd_key2]

            if bnd_key in frm_bnd_keys:
                # Check if radical atom is closer to some atom
                # other than the bonding atom
                cls = is_atom_closest_to_bond_atom(
                    zma, bnd_key2, cnf_dist)
                if not cls:
                    print('distance', ref_dist, cnf_dist)
                    print(
                        ' - Radical atom now has a new nearest neighbor')
                    viable = False
                # check forming bond distance
                if abs(cnf_dist - ref_dist) > max_disp:
                    print('distance', ref_dist, cnf_dist)
                    viable = False

            # Check distance relative to equi. bond
            equi_bnd = automol.util.dict_.values_by_unordered_tuple(
                bnd.LEN_DCT, (symb1, symb2), fill_val=0.0)
            displace_from_equi = cnf_dist - equi_bnd
            dchk1 = abs(cnf_dist - ref_dist) > 0.1
            dchk2 = displace_from_equi < 0.2
            if dchk1 and dchk2:
                print(cnf_dist, equi_bnd)
                viable = False
        else:
            # check forming/breaking bond distance
            # if abs(cnf_dist - ref_dist) > 0.4:
            # max disp of 0.4 causes problems for bond scission w/ ring forming
            # not sure if setting it to 0.3 will cause problems for other cases
            if abs(cnf_dist - ref_dist) > 0.3:
                print('distance', ref_dist, cnf_dist)
                viable = False

    return viable


# Assessments of two geomeries
def minimum_distance(geo1, geo2):
    """ get the minimum distance between atoms in geo1 and those in geo2

        :param geo1: molecular geometry 1
        :type geo1: automol molecular geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol molecular geometry data structure
        :rtype: float
    """

    xyzs1 = geom_base.coordinates(geo1)
    xyzs2 = geom_base.coordinates(geo2)
    return min(util.vec.distance(xyz1, xyz2)
               for xyz1, xyz2 in itertools.product(xyzs1, xyzs2))


# Calculating quantities used for comparisons
def coulomb_spectrum(geo):
    """ Calculate a Coulomb matrix eigenvalue spectrum where
        the eignevalues are sorted in ascending order.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: tuple(float)
    """

    mat = _coulomb_matrix(geo)
    vals = tuple(sorted(numpy.linalg.eigvalsh(mat)))

    return vals


def _coulomb_matrix(geo):
    """ Calculate the Coulomb matrix wich describes the
        electrostatic interactions between nuclei:

        M[i,j] = 0.5Z_i^2.4 (i=j), Z_iZ_j/R_ij (i!=j

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: tuple(tuple(float))
    """

    nums = numpy.array(list(map(ptab.to_number, geom_base.symbols(geo))))
    xyzs = numpy.array(geom_base.coordinates(geo))

    _ = numpy.newaxis
    natms = len(nums)
    diag_idxs = numpy.diag_indices(natms)
    tril_idxs = numpy.tril_indices(natms, -1)
    triu_idxs = numpy.triu_indices(natms, 1)

    zxz = numpy.outer(nums, nums)
    rmr = numpy.linalg.norm(xyzs[:, _, :] - xyzs[_, :, :], axis=2)

    mat = numpy.zeros((natms, natms))
    mat[diag_idxs] = nums ** 2.4 / 2.
    mat[tril_idxs] = zxz[tril_idxs] / rmr[tril_idxs]
    mat[triu_idxs] = zxz[triu_idxs] / rmr[triu_idxs]

    return mat


def distance_matrix(geo):
    """ Form a Natom X Natom matrix containing the distance of all the
        atoms in a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: numpy.ndarray
    """

    mat = numpy.zeros((len(geo), len(geo)))
    for i in range(len(geo)):
        for j in range(len(geo)):
            mat[i][j] = geom_base.distance(geo, i, j)

    return mat


# Comparisons between two geometries
def almost_equal_coulomb_spectrum(geo1, geo2, rtol=1e-2):
    """ Assess if two molecular geometries have similar coulomb spectrums.

        :param geo1: molecular geometry 1
        :type geo1: automol molecular geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol molecular geometry data structure
        :param rtol: Relative tolerance for the distances
        :type rtol: float
        :rtype: bool
    """
    return numpy.allclose(
        coulomb_spectrum(geo1), coulomb_spectrum(geo2), rtol=rtol)


def argunique_coulomb_spectrum(geos, seen_geos=(), rtol=1e-2):
    """ Get indices of unique geometries, by coulomb spectrum.

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
    """ Get the indices of unique items using some comparison function.

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


def almost_equal(geo1, geo2, rtol=2e-3):
    """ Assess if the coordinates of two molecular geometries
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
    if geom_base.symbols(geo1) == geom_base.symbols(geo2):
        ret = numpy.allclose(
            geom_base.coordinates(geo1),
            geom_base.coordinates(geo2),
            rtol=rtol)

    return ret


def almost_equal_dist_matrix(geo1, geo2, thresh=0.1):
    """form distance matrix for a set of xyz coordinates
    """

    natoms = len(geo1)
    for i in range(natoms):
        for j in range(natoms):
            dist_mat1_ij = geom_base.distance(geo1, i, j)
            dist_mat2_ij = geom_base.distance(geo2, i, j)
            if abs(dist_mat1_ij - dist_mat2_ij) > thresh:
                return False

    return True


def are_torsions_same2(geo, geoi, idxs_lst):
    """ Are torsions the same with torsions identified
        by a list of 1x4 lists of atom indices
    """
    dtol = 0.09
    same_dihed = True
    for idxs in idxs_lst:
        val = automol.geom.dihedral_angle(geo, *idxs)
        vali = automol.geom.dihedral_angle(geoi, *idxs)
        valip = vali+2.*numpy.pi
        valim = vali-2.*numpy.pi
        vchk1 = abs(val - vali)
        vchk2 = abs(val - valip)
        vchk3 = abs(val - valim)
        if vchk1 > dtol and vchk2 > dtol and vchk3 > dtol:
            same_dihed = False
    return same_dihed


def are_torsions_same(geo, geoi, ts_bnds=()):
    """ compare all torsional angle values
    """

    dtol = 0.09
    same_dihed = True

    # Build the Z-Matrix torsion names
    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(
        geo, ts_bnds=ts_bnds)
    zmai = automol.geom.zmatrix(geoi)
    tors_namesi = automol.geom.zmatrix_torsion_coordinate_names(
        geoi, ts_bnds=ts_bnds)

    # Compare the torsions
    for idx, tors_name in enumerate(tors_names):
        val = value_dictionary(zma)[tors_name]
        vali = value_dictionary(zmai)[tors_namesi[idx]]
        valip = vali+2.*numpy.pi
        valim = vali-2.*numpy.pi
        vchk1 = abs(val - vali)
        vchk2 = abs(val - valip)
        vchk3 = abs(val - valim)
        if vchk1 > dtol and vchk2 > dtol and vchk3 > dtol:
            same_dihed = False

    return same_dihed


# Checks
def is_unique(geo, geo_lst, check_dct=None):
    """ Compare one of many structure features of a geometry to that of
        a list of geometries to see if it is unique.

        order of atoms also impacts the comparison as well
    """

    # Set default check values if none are provided
    if check_dct is None:
        check_dct = CHECK_DEFAULT_DCT

    unique = True
    like_idx = None
    for idx, geoi in enumerate(geo_lst):
        # Perform all of the desired comparison checks for similarity
        sim_chk_results = []
        for key, val in check_dct.items():
            kwargs = {'arg': val} if val is not None else {}
            sim_chk_results.append(CHECK_FXN_DCT[key](geo, geoi, **kwargs))
        # If all checks come back as True, than geoms are the same
        if all(sim_chk_results):
            unique = False
            like_idx = idx
            break

    return unique, like_idx


def _similar_dist(geo, geoi, arg=3e-1):
    """ Compare the distance matrices of two geometries.
    """
    return almost_equal_dist_matrix(geo, geoi, thresh=arg)


def _similar_tors(geo, geoi, arg=()):
    """ Compare the torsions of two geometries
    """
    return are_torsions_same(geo, geoi, ts_bnds=arg)


def _similar_stereo(geo, geoi, arg=None):
    """ Compare the stereochemistry of two geometries
    """
    _ = arg  # Added just to make wrapper function work
    ich = automol.convert.geom.inchi(geo)
    ichi = automol.convert.geom.inchi(geoi)
    return bool(ich == ichi)


def _similar_coulomb(geo, geoi, arg=1e-2):
    """ Compare the Coulomb spectrum of geometries.
    """
    return almost_equal_coulomb_spectrum(geo, geoi, rtol=arg)


CHECK_DEFAULT_DCT = {
    'dist': 3.5e-1,
    'coulomb': 1.5e-2,
    'stereo': None,
    'tors': None
}

CHECK_FXN_DCT = {
    'dist': _similar_dist,
    'tors': _similar_tors,
    'stereo': _similar_stereo,
    'coulomb': _similar_coulomb
}
