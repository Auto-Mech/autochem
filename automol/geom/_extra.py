""" extra geometry functions
"""

import automol.convert.geom
import automol.graph
from automol.geom import _trans as trans


def int_sym_num_from_sampling(sym_geos, grxn=None):
    """ Determine the symmetry number for a given conformer geometry.
    (1) Explore the saved conformers to find the list of similar conformers -
        i.e. those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.
    """

    if grxn is not None:
        frm_bnd_keys = automol.reac.forming_bond_keys(grxn)
        brk_bnd_keys = automol.reac.breaking_bond_keys(grxn)
    else:
        frm_bnd_keys, brk_bnd_keys = frozenset({}), frozenset({})

    int_sym_num = 0
    # modify geometries to remove H's from rotatable XHn end group
    # this will be accounted for separately as multiplicative factor
    mod_sym_geos = []
    for geo_sym_i in sym_geos:
        ret = automol.geom.end_group_symmetry_factor(
            geo_sym_i, frm_bnd_keys, brk_bnd_keys)
        mod_geo_sym_i, end_group_factor = ret
        # ioprinter.info_message('end_group_factor test:', end_group_factor)

        new_geom = True
        for mod_geo_sym_j in mod_sym_geos:
            if automol.geom.almost_equal_dist_matrix(
                    mod_geo_sym_i, mod_geo_sym_j, thresh=3e-1):
                if grxn is not None:
                    new_geom = False
                    break
                tors_same = automol.geom.are_torsions_same(
                    mod_geo_sym_i, mod_geo_sym_j, ts_bnds=())
                if tors_same:
                    new_geom = False
                    break
        if new_geom:
            mod_sym_geos.append(mod_geo_sym_i)
            int_sym_num += 1

    int_sym_num *= end_group_factor

    return int_sym_num


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
    neighbor_dct = automol.graph.atoms_neighbor_atom_keys(gra)

    ts_atms = []
    for bnd in frm_bnd_keys:
        ts_atms.extend(list(bnd))
    for bnd in brk_bnd_keys:
        ts_atms.extend(list(bnd))
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
            if atm not in ts_atms:
                nonh_neighs = []
                h_neighs = []
                neighs = neighbor_dct[atm]
                for nei in neighs:
                    if nei in all_hyds:
                        h_neighs.append(nei)
                    else:
                        nonh_neighs.append(nei)
                if len(nonh_neighs) == 1 and len(h_neighs) > 1:
                    term_atms[atm] = h_neighs
                    print('terminal atom accepted:', atm, h_neighs)
    factor = 1.
    remove_atms = []
    for atm in term_atms:
        hyds = term_atms[atm]
        if len(hyds) > 1:
            factor *= len(hyds)
            remove_atms.extend(hyds)
    geo = trans.remove_coordinates(geo, remove_atms)

    return geo, factor, remove_atms


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
    neighbor_dct = automol.graph.atoms_neighbor_atom_keys(gra)
    ts_atms = []
    for bnd in frm_bnd_keys:
        ts_atms.extend(list(bnd))
    for bnd in brk_bnd_keys:
        ts_atms.extend(list(bnd))

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
            print('atm test:', atm, ts_atms)
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
