""" extra geometry functions
"""

import automol.convert.geom
import automol.graph
from automol.geom import _trans as trans


def end_group_symmetry_factor(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Determine sym factor for terminal groups in a geometry

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param frm_bnd_keys: keys denoting atoms forming bond in TS
        :type frm_bnd_keys: frozenset(int)
        :param brk_bnd_keys: keys denoting atoms breaking bond in TS
        :type brk_bnd_keys: frozenset(int)
        :rtype: (automol geom data structure, float)
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = automol.convert.geom.graph(geo, stereo=False)
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
    geo = trans.remove(geo, remove_atms)

    return geo, factor


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

    gra = automol.convert.geom.graph(geo, stereo=False)
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
            new_geo = trans.swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = trans.swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
            new_geo = trans.swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = trans.swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
        else:
            geo_lst.append(new_geo)
            new_geo = trans.swap_coordinates(new_geo, hyds[0], hyds[1])
            geo_lst.append(new_geo)

    return geo_lst
