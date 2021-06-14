""" Handle symmetry factor stuff
"""

import automol


# internal symmetry number
def internal_symm_from_sampling(symm_geos, rotors, grxn=None, zma=None):
    """ Determines the internal symmetry number for a given conformer geometry
        by assessing a set of symmetrically similar structures that have been
        obtained by previous conformer sampling processes.

        (1) Explore saved conformers to find the list of similar conformers,
            i.e., those with a coulomb matrix and energy that are equivalent
            to those for the reference geometry.
        (2) Expand each of those similar conformers by applying
            rotational permutations to each of the terminal groups.
        (3) Count how many distinct distance matrices there are in
            the fully expanded conformer list.

        :param symm_geos: geometries that are symmetrically similar to one another
    """

    if grxn is not None:
        frm_bnd_keys = automol.reac.forming_bond_keys(grxn)
        brk_bnd_keys = automol.reac.breaking_bond_keys(grxn)
        tors_names = automol.rotor.names(rotors, flat=True)
        tors_idxs = [automol.zmat.coord_idxs(zma, name) for name in tors_names]
    else:
        frm_bnd_keys, brk_bnd_keys = frozenset({}), frozenset({})

    # Modify geometries to remove H's from rotatable XHn end group;
    # this will be accounted for separately as multiplicative factor
    int_sym_num = 0
    mod_symm_geos = []
    for geo_sym_i in symm_geos:
        ret = automol.geom.end_group_symmetry_factor(
            geo_sym_i, frm_bnd_keys, brk_bnd_keys)
        mod_geo_sym_i, end_group_factor, removed_atms = ret
        if grxn is not None:
            mod_tors_idxs = _modify_idxs(
                tors_idxs, removed_atms, automol.zmat.dummy_keys(zma))

        new_geom = True
        for mod_geo_sym_j in mod_symm_geos:
            if automol.geom.almost_equal_dist_matrix(
                    mod_geo_sym_i, mod_geo_sym_j, thresh=3e-1):
                if grxn is None:
                    tors_same = automol.geom.are_torsions_same(
                        mod_geo_sym_i, mod_geo_sym_j, ts_bnds=())
                else:
                    tors_same = automol.geom.are_torsions_same2(
                        mod_geo_sym_i, mod_geo_sym_j, mod_tors_idxs)
                if tors_same:
                    new_geom = False
                    break
        if new_geom:
            mod_symm_geos.append(mod_geo_sym_i)
            int_sym_num += 1

    int_sym_num *= end_group_factor

    return int_sym_num, end_group_factor


def reduce_internal_symm(int_symm, ext_symm, end_group_factor, geo):
    """ Reduce symm if external sym is 3??
    """
    if ext_symm % 3 == 0 and end_group_factor > 1:
        if not automol.graph.is_branched(automol.geom.graph(geo)):
            int_symm = int_symm / 3

    return int_symm


def rotor_reduced_symm_factor(sym_factor, rotor_symms):
    """ Decrease the overall molecular symmetry factor by the
        torsional mode symmetry numbers
    """
    for symm in rotor_symms:
        sym_factor /= symm

    return sym_factor


# Helper functions
def _modify_idxs(idxs_lst, removed_atms, dummy_atms):
    mod_idxs_lst = []
    no_dummy_idxs_lst = []
    for idxs in idxs_lst:
        mod_idxs = []
        for idx in idxs:
            mod_idx = idx
            for atm in dummy_atms:
                if atm < idx:
                    mod_idx -= 1
            mod_idxs.append(mod_idx)
        no_dummy_idxs_lst.append(mod_idxs)

    for idxs in no_dummy_idxs_lst:
        in_lst = True
        for idx in idxs:
            if idx in removed_atms:
                in_lst = False
        if in_lst:
            mod_idxs = []
            for idx in idxs:
                mod_idx = idx
                for atm in removed_atms:
                    if atm < idx:
                        mod_idx -= 1
                mod_idxs.append(mod_idx)
            mod_idxs_lst.append(mod_idxs)
    return mod_idxs_lst
