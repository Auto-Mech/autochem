""" extra high-level geometry library functions
"""
import numpy
import automol.graph
import automol.zmat.base
from automol.geom._conv import graph
from automol.geom._conv import inchi
from automol.geom._conv import x2z_zmatrix
from automol.geom._conv import x2z_torsion_coordinate_names
from automol.geom.base import swap_coordinates
from automol.geom.base import dihedral_angle
from automol.geom.base import central_angle
from automol.geom.base import almost_equal_dist_matrix
from automol.geom.base import almost_equal_coulomb_spectrum
from automol.geom.base import distance_matrix
from automol.geom.base import count


CHECK_DEFAULT_DCT = {
    'dist': 3.5e-1,
    'coulomb': 1.5e-2,
    'stereo': None,
    'tors': None
}


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
    ich = inchi(geo)
    ichi = inchi(geoi)
    return bool(ich == ichi)


def _similar_coulomb(geo, geoi, arg=1e-2):
    """ Compare the Coulomb spectrum of geometries.
    """
    return almost_equal_coulomb_spectrum(geo, geoi, rtol=arg)


CHECK_FXN_DCT = {
    'dist': _similar_dist,
    'tors': _similar_tors,
    'stereo': _similar_stereo,
    'coulomb': _similar_coulomb
}


def components_graph(geo, stereo=True):
    """ Generate a list of molecular graphs where each element is a graph that
        consists of fully connected (bonded) atoms. Stereochemistry is included
        if requested.
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph data structure
    """
    return automol.graph.connected_components(graph(geo, stereo=stereo))


def connected(geo, stereo=True):
    """ Determine if all atoms in geometry are completely connected.
        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: bool
    """
    return len(components_graph(geo, stereo=stereo)) == 1


def rot_permutated_geoms(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Convert an input geometry to a list of geometries
        corresponding to the rotational permuations of all the terminal groups.
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param frm_bnd_keys: keys denoting atoms forming bond in TS
        :type frm_bnd_keys: frozenset(int)
        :param brk_bnd_keys: keys denoting atoms breaking bond in TS
        :type brk_bnd_keys: frozenset(int)
        :rtype: tuple(automol geom data structure)
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = graph(geo, stereo=False)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atoms_neighbor_atom_keys(gra)
    ts_atms = []
    for bnd_ in frm_bnd_keys:
        ts_atms.extend(list(bnd_))
    for bnd_ in brk_bnd_keys:
        ts_atms.extend(list(bnd_))

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.radical_atom_keys(gra, sing_res=True)
        res_rad_atms = automol.graph.radical_atom_keys(gra)
        rad_atms = [atm for atm in rad_atms if atm not in res_rad_atms]
    else:
        rad_atms = []

    gra = gra[0]
    for atm in gra:
        if gra[atm][0] == 'H':
            all_hyds.append(atm)
    for atm in gra:
        if atm in unsat_atms and atm not in rad_atms:
            pass
        else:
            if atm not in ts_atms:
                nonh_neighs = []
                h_neighs = []
                neighs = neighbor_dct[atm]
                for nei in neighs:
                    if nei in all_hyds:
                        h_neighs.append(nei)
                    else:
                        nonh_neighs.append(nei)
                if len(nonh_neighs) < 2 and len(h_neighs) > 1:
                    term_atms[atm] = h_neighs
    geo_final_lst = [geo]
    for hyds in term_atms.values():
        geo_lst = []
        for geom in geo_final_lst:
            geo_lst.extend(_swap_for_one(geom, hyds))
        geo_final_lst = geo_lst

    return geo_final_lst


def _swap_for_one(geo, hyds):
    """ Rotational permuation for one rotational group.
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param hyd: list of hydrogen atom indices
        :type hyd: tuple(int)
        :rtype: tuple(automol molecular geometry data structure)
    """

    geo_lst = []
    if len(hyds) > 1:
        new_geo = geo
        if len(hyds) > 2:
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
        else:
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            geo_lst.append(new_geo)

    return geo_lst


def are_torsions_same2(geo, geoi, idxs_lst):
    """ Are torsions the same with torsions identified
        by a list of 1x4 lists of atom indices
    """
    dtol = 0.09
    same_dihed = True
    for idxs in idxs_lst:
        val = dihedral_angle(geo, *idxs)
        vali = dihedral_angle(geoi, *idxs)
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
    zma = x2z_zmatrix(geo, ts_bnds=ts_bnds)
    tors_names = x2z_torsion_coordinate_names(
        geo, ts_bnds=ts_bnds)
    zmai = x2z_zmatrix(geoi)
    tors_namesi = x2z_torsion_coordinate_names(
        geoi, ts_bnds=ts_bnds)

    # Compare the torsions
    for idx, tors_name in enumerate(tors_names):
        val = automol.zmat.base.value_dictionary(zma)[tors_name]
        vali = automol.zmat.base.value_dictionary(zmai)[tors_namesi[idx]]
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


def hydrogen_bonded_structure(
        geo, dist_thresh=4.82, angle_thresh=1.92,
        grxn=None):
    """ Compare bond lengths in structure to determine if there
        is a hydrogen bond

        :param geo: geometry object
        :type geo: geo object (tuple of tuples)
        :param grxn: reaction object (geo indexing)
        :type grxn: automol.reac.Reaction object
        :param dist_thresh: cutoff value for hbond length (Bohr)
        :type dist_thresh: float
        :param angle_thresh: cutoff value for hbond angle (Radian)
        :type angle_thresh: float
        :rtype: boolean
    """
    hydrogen_bond = hydrogen_bonded_idxs(
        geo, dist_thresh, angle_thresh, grxn)
    return hydrogen_bond is not None


def hydrogen_bonded_idxs(
        geo, dist_thresh=5.3, angle_thresh=1.92,
        grxn=None):
    """ Compare bond lengths in structure to determine if there
        is a hydrogen bond.

        :param geo: geometry object
        :type geo: geo object (tuple of tuples)
        :param grxn: reaction object (geo indexing)
        :type grxn: automol.reac.Reaction object
        :param dist_thresh: cutoff value for hbond length (Bohr)
        :type dist_thresh: float
        :param angle_thresh: cutoff value for hbond angle (Radian)
        :type angle_thresh: float
        :rtype: tuple
    """
    # Initialize the hydrogen bond list to None
    hydrogen_bond = None
    if count(geo) > 1:
        # Get the forming/breaking bond idxs if possible
        if grxn is not None:
            frm_bnd_keys = automol.graph.ts.ts_forming_bond_keys(
                grxn.ts_graph)
            brk_bnd_keys = automol.graph.ts.ts_breaking_bond_keys(
                grxn.ts_graph)
            rxn_keys = set()
            for key in frm_bnd_keys:
                rxn_keys = rxn_keys | key
            for key in brk_bnd_keys:
                rxn_keys = rxn_keys | key
            rxn_h_idxs = tuple(rxn_keys)
        else:
            rxn_h_idxs = ()

        # Get all potential indices for HB interactions
        gra = graph(geo)
        dist_mat = distance_matrix(geo)
        adj_atm_dct = automol.graph.atoms_neighbor_atom_keys(gra)
        h_idxs = automol.graph.atom_keys(gra, symb='H')
        acceptor_idxs = list(automol.graph.radical_atom_keys(gra))
        acceptor_idxs.extend(list(automol.graph.atom_keys(gra, symb='O')))
        # Loop over indices, ignoring H-idxs in reacting bonds
        hb_idxs = tuple(idx for idx in h_idxs
                        if idx not in rxn_h_idxs)
        for h_idx in hb_idxs:
            for acceptor_idx in acceptor_idxs:
                donor_idx = list(adj_atm_dct[h_idx])[0]
                if acceptor_idx in adj_atm_dct[donor_idx]:
                    continue
                if dist_mat[h_idx][acceptor_idx] < dist_thresh:
                    ang = central_angle(geo, donor_idx, h_idx, acceptor_idx)
                    if ang > angle_thresh:
                        hydrogen_bond = (donor_idx, h_idx, acceptor_idx,)
                        dist_thresh = dist_mat[h_idx][acceptor_idx]
    return hydrogen_bond
