""" Handle symmetry factor stuff
"""


import random
import string

from automol import geom, graph, reac, rotor, zmat


# external
def external_symm(geo):
    """Determine the external symmetry factor which relates to the external
    point group symmetry
    """
    return geom.external_symmetry_factor(geo)


# internal symmetry number
def internal_symm_from_sampling(symm_geos, rotors, grxn=None, zma=None):
    """Determines the internal symmetry number for a given conformer geometry
    by assessing a set of symmetrically similar structures that have been
    obtained by previous conformer sampling processes.

    (1) Explore saved conformers to find the list of similar conformers,
        i.e., those with a coulomb matrix and energy that are equivalent
        to those for the reference geometry.
    (2) Expand each of those similar conformers by applying
        rotational permutations to each of the terminal groups.
    (3) Count how many distinct distance matrices there are in
        the fully expanded conformer list.

    :param symm_geos: geometries symmetrically similar to one another
    """

    if grxn is not None:
        frm_bnd_keys = graph.ts.forming_bond_keys(reac.ts_graph(grxn))
        brk_bnd_keys = graph.ts.breaking_bond_keys(reac.ts_graph(grxn))
        tors_names = rotor.names(rotors, flat=True)
        tors_idxs = [zmat.coord_idxs(zma, name) for name in tors_names]
    else:
        frm_bnd_keys, brk_bnd_keys = frozenset({}), frozenset({})

    # Modify geometries to remove H's from rotatable XHn end group;
    # this will be accounted for separately as multiplicative factor
    int_sym_num = 0
    mod_symm_geos = []
    for geo_sym_i in symm_geos:
        ret = end_group_symmetry_factor(geo_sym_i, frm_bnd_keys, brk_bnd_keys)
        mod_geo_sym_i, end_group_factor, removed_atms = ret
        if grxn is not None:
            mod_tors_idxs = _modify_idxs(tors_idxs, removed_atms, zmat.dummy_keys(zma))

        new_geom = True
        for mod_geo_sym_j in mod_symm_geos:
            if geom.almost_equal_dist_matrix(mod_geo_sym_i, mod_geo_sym_j, thresh=3e-1):
                if grxn is None:
                    tors_same = geom.are_torsions_same(
                        mod_geo_sym_i, mod_geo_sym_j, ts_bnds=()
                    )
                else:
                    tors_same = geom.are_torsions_same2(
                        mod_geo_sym_i, mod_geo_sym_j, mod_tors_idxs
                    )
                if tors_same:
                    new_geom = False
                    break
        if new_geom:
            mod_symm_geos.append(mod_geo_sym_i)
            int_sym_num += 1

    int_sym_num *= end_group_factor

    return int_sym_num, end_group_factor


def reduce_internal_symm(geo, int_symm, ext_symm, end_group_factor):
    """Reduce symm if external sym is 3??"""
    if ext_symm % 3 == 0 and end_group_factor > 1:
        if not graph.is_branched(geom.graph(geo)):
            int_symm = int_symm / 3

    return int_symm


def rotor_reduced_symm_factor(sym_factor, rotor_symms):
    """Decrease the overall molecular symmetry factor by the
    torsional mode symmetry numbers
    """
    for symm in rotor_symms:
        sym_factor /= symm

    return sym_factor


# End group symmetry calculators
def end_group_symmetry_factor(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """Determine sym factor for terminal groups in a geometry
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

    gra = geom.graph(geo, stereo=False)
    term_atms = {}
    all_hyds = []
    neighbor_dct = graph.atoms_neighbor_atom_keys(gra)

    ts_atms = []
    for bnd_ in frm_bnd_keys:
        ts_atms.extend(list(bnd_))
    for bnd_ in brk_bnd_keys:
        ts_atms.extend(list(bnd_))
    # determine if atom is a part of a double bond
    unsat_atms = graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = graph.radical_atom_keys(gra, sing_res=True)
        res_rad_atms = graph.radical_atom_keys(gra)
        rad_atms = [atm for atm in rad_atms if atm not in res_rad_atms]
    else:
        rad_atms = []

    gra = gra[0]
    for atm in gra:
        if gra[atm][0] == "H":
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
    factor = 1.0
    remove_atms = []
    for hyds in term_atms.values():
        if len(hyds) > 1:
            factor *= len(hyds)
            remove_atms.extend(hyds)
    geo = geom.remove(geo, remove_atms)

    return geo, factor, remove_atms


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


ICH_DCT = {"C": "InChI=1S/C", "O": "InChI=1S/O"}


def oxygenated_hydrocarbon_symm_num(geo, zrxn=None, racemic=True):
    """determine the symmetry number of a CHO molecule"""
    int_symm = 1.0
    chiral_center = 0
    if zrxn is not None:
        gra = reac.ts_graph(zrxn)
    else:
        gra = geom.graph(geo)
    ethane_gra = (
        {0: ("C", 3, None), 1: ("C", 3, None)},
        {frozenset({0, 1}): (1, None)},
    )
    if graph.implicit(gra) == ethane_gra:
        int_symm = 3.0
        ext_symm = geom.external_symmetry_factor(geo)
    else:
        gra = graph.explicit(gra)
        atms = graph.atom_keys(gra)
        atm_vals = graph.unsaturated_atom_keys(gra)
        ring_atms = graph.rings_atom_keys(gra)
        ring_atms = [x for ring in ring_atms for x in ring]
        atm_rads = graph.nonresonant_radical_atom_keys(gra)
        atm_syms = graph.atom_symbols(gra)
        atms = [x for x in atms if atm_syms[x] != "H"]
        for atm in atms:
            if atm in atm_vals and atm not in atm_rads:
                continue
            if atm in ring_atms:
                atm_groups = graph.ring_atom_chirality(gra, atm, ring_atms)
            else:
                atm_groups = graph.branches(gra, atm)
            group_dct = {}
            for group in atm_groups:
                try:
                    group_ich = graph.inchi(group)
                except Exception as err:
                    # Excepts rdkit errors, assumes group is complicated enough
                    # that is is unique
                    if not str(err).startswith("Python argument types in"):
                        print("Error evaluating atom group in symm number routine")
                        print(
                            "Symmetry number may be incorrect as a result, group is",
                            group,
                        )
                    group_ich = "".join(
                        random.choice(string.ascii_letters) for i in range(10)
                    )
                if group_ich in group_dct:
                    group_dct[group_ich] += 1
                else:
                    group_dct[group_ich] = 1
            # remove atom inchi from dct
            group_dct = {x: y for x, y in group_dct.items() if y != 0}
            if len(group_dct) == 4:
                chiral_center += 1
            elif len(group_dct) == 3 and atm in atm_rads:
                chiral_center += 1
            if atm in ring_atms:
                continue
            if len(group_dct) == 2:
                chain_group = None
                symm_groups = None
                for group in group_dct.keys():
                    if group_dct[group] == 1:
                        if group not in ("InChI=1S/H", "InChI=1S/O"):
                            chain_group = group
                    else:
                        symm_groups = group
                if chain_group and symm_groups:
                    atm_symm = group_dct[symm_groups]
                    # if atm_symm == 2:
                    #     if chain_group != 'InChI=1S/H':
                    #         atm_symm = 1
                    int_symm *= atm_symm
        ext_symm = geom.external_symmetry_factor(
            geo, chiral_center=chiral_center > 0.0 and racemic
        )

    return int_symm, ext_symm
