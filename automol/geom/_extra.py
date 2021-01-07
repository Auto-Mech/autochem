""" extra geometry functions
"""
import numpy
import automol.convert.geom
import automol.graph
from automol.geom._geom import remove
from automol.geom._geom import swap_coordinates
from automol.geom._geom import distance


def end_group_sym_factor(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Determine sym factor for terminal groups in a geometry
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = automol.convert.geom.graph(geo, remove_stereo=True)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atom_neighbor_keys(gra)

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.sing_res_dom_radical_atom_keys(gra)
        res_rad_atms = automol.graph.resonance_dominant_radical_atom_keys(gra)
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
            if atm not in frm_bnd_keys and atm not in brk_bnd_keys:
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
    factor = 1.
    remove_atms = []
    for atm in term_atms:
        hyds = term_atms[atm]
        if len(hyds) > 1:
            factor *= len(hyds)
            remove_atms.extend(hyds)
    geo = remove(geo, remove_atms)
    return geo, factor


def rot_permutated_geoms(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ convert an input geometry to a list of geometries
        corresponding to the rotational permuations of all the terminal groups
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = automol.convert.geom.graph(geo, remove_stereo=True)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atom_neighbor_keys(gra)

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.sing_res_dom_radical_atom_keys(gra)
        res_rad_atms = automol.graph.resonance_dominant_radical_atom_keys(gra)
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
            if atm not in frm_bnd_keys and atm not in brk_bnd_keys:
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
    for atm in term_atms:
        hyds = term_atms[atm]
        geo_lst = []
        for geom in geo_final_lst:
            geo_lst.extend(_swap_for_one(geom, hyds))
        geo_final_lst = geo_lst
    return geo_final_lst


def _swap_for_one(geo, hyds):
    """ rotational permuation for one rotational group
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


def almost_equal_dist_matrix(geo1, geo2, thresh=0.1):
    """form distance matrix for a set of xyz coordinates
    """
    for i in range(len(geo1)):
        for j in range(len(geo1)):
            dist_mat1_ij = distance(geo1, i, j)
            dist_mat2_ij = distance(geo2, i, j)
            if abs(dist_mat1_ij - dist_mat2_ij) > thresh:
                return False
    return True


def find_xyzp_using_internals(xyz1, xyz2, xyz3, pdist, pangle, pdihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates

    TODO: This already exists: util.vec.from_internals does exactly the same
    thing; find where this is used and replace it with that
    """

    # Set to numpy arrays
    xyz1 = numpy.array(xyz1)
    xyz2 = numpy.array(xyz2)
    xyz3 = numpy.array(xyz3)

    # Set the coordinates of Point P in the RT system
    xyzp_rt = numpy.array([pdist * numpy.sin(pangle) * numpy.cos(pdihed),
                           pdist * numpy.cos(pangle),
                           -(pdist * numpy.sin(pangle) * numpy.sin(pdihed))
                           ])

    # Set the coordinates of the Point 2 and 3 in the RT system
    dist12 = numpy.linalg.norm(xyz1 - xyz2)
    dist13 = numpy.linalg.norm(xyz1 - xyz3)
    dist23 = numpy.linalg.norm(xyz2 - xyz3)
    xyz2_rt = numpy.array([0.0, dist12, 0.0])

    val = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    valx3 = numpy.sqrt(dist13**2 - val**2)
    valy3 = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    xyz3_rt = numpy.array([valx3, valy3, 0.0])

    # Translate original frame of ref coors so that xyz1 is at (0, 0, 0)
    xyz2_t = xyz2 - xyz1
    xyz3_t = xyz3 - xyz1

    # Rotation matrix to rotate back to the original ref system
    r12 = (xyz2[0] - xyz1[0]) / xyz2_rt[1]
    r22 = (xyz2[1] - xyz1[1]) / xyz2_rt[1]
    r32 = (xyz2[2] - xyz1[2]) / xyz2_rt[1]

    r11 = (xyz3[0] - xyz1[0] - xyz3_rt[1]*r12) / xyz3_rt[0]
    r21 = (xyz3[1] - xyz1[1] - xyz3_rt[1]*r22) / xyz3_rt[0]
    r31 = (xyz3[2] - xyz1[2] - xyz3_rt[1]*r32) / xyz3_rt[0]

    anum_aconst = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
    den_aconst = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]

    if abs(anum_aconst) < 1.0e-6 and abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    elif abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    else:
        anum = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
        aden = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]
        aconst = anum / aden

    den1 = (xyz3_t[1] / xyz3_t[0]) - aconst * (xyz3_t[2] / xyz3_t[0])
    if den1 == 0.0:
        den1 = 1.0e-20
    bconst = 1.0 / den1

    # Set vals for another point
    valx = -(1.0 / numpy.sqrt(1.0 + (bconst**2) * (1.0 + aconst**2)))
    valy = -(valx * bconst)
    xyz4_t = numpy.array([valx, valy, -(valy * aconst)])

    r13 = xyz4_t[0]
    r23 = xyz4_t[1]
    r33 = xyz4_t[2]
    r13n = -r13
    r23n = -r23
    r33n = -r33

    # Now rotate and translate back
    # Here I check  the (001) vector direction to decide whether
    # To take the positive of negative results of square root taken above
    xap = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13 * xyzp_rt[2]))
    yap = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))
    zap = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))

    xan = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13n * xyzp_rt[2]))
    yan = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r23n * xyzp_rt[2]))
    zan = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33n * xyzp_rt[2]))

    bvec = xyz1 - xyz2
    cvec = xyz2 - xyz3
    vec1 = (bvec[1] * cvec[2]) - (bvec[2] * cvec[1])
    vec2 = (bvec[2] * cvec[0]) - (bvec[0] * cvec[2])
    vec3 = (bvec[0] * cvec[1]) - (bvec[1] * cvec[0])

    if abs(xyz4_t[0]) > 1.0e-5:
        checkv = vec1 / xyz4_t[0]
    elif abs(xyz4_t[1]) > 1.0e-5:
        checkv = vec2 / xyz4_t[1]
    else:
        checkv = vec3 / xyz4_t[2]

    if checkv >= 0.0:
        xyzp = numpy.array([xap, yap, zap])
    else:
        xyzp = numpy.array([xan, yan, zan])

    return xyzp[0], xyzp[1], xyzp[2]
