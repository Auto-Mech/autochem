""" Handle symmetry factor stuff
"""
import random
import string

from automol import geom, graph, reac, zmat
from automol.data import rotor
from automol.util import zmat_conv


# internal symmetry number
def symmetry_factors_from_sampling(geos, rotors, grxn=None):
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
    # First, get the end-group symmetry factor
    geo0 = geos[0]
    gra = geom.graph(geo0, stereo=False) if grxn is None else reac.ts_graph(grxn)
    end_sym_fac = end_group_symmetry_factor(geos[0], gra=gra)

    # Get indices for the torsions
    zma = rotor.rotors_zmatrix(rotors)
    zc_ = zmat.conversion_info(zma)
    tors_names = rotor.rotors_torsion_names(rotors, flat=True)
    tors_zidxs = [zmat.coordinate(zma, n) for n in tors_names]
    tors_gidxs = zmat_conv.relabel_zmatrix_key_sequence(zc_, tors_zidxs)

    # Ignore terminal hydrogens
    gra = graph.implicit(gra, atm_keys=graph.terminal_atom_keys(gra))
    idx_pool = graph.atom_keys(gra)
    tors_gidxs = [ixs for ixs in tors_gidxs if set(ixs) <= idx_pool]

    # For filtering out identical (not just symmetrically equivalent) structures
    def _are_equivalent(geo1, geo2):
        same_dist = geom.almost_equal_dist_matrix(
            geo1, geo2, thresh=3e-1, idxs=idx_pool
        )
        same_tors = geom.are_torsions_same(geo1, geo2, tors_gidxs)
        return same_dist and same_tors

    # Identify unique geometries
    uniq_geos = []
    for geo in geos:
        if not any(_are_equivalent(geo, g) for g in uniq_geos):
            uniq_geos.append(geo)

    int_sym_fac = len(uniq_geos) * end_sym_fac

    return int_sym_fac, end_sym_fac


def end_group_symmetry_factor(geo, gra=None) -> float:
    """Determine the symmetry factor for terminal groups in a geometry

    :param geo: A geometry
    :type geo: automol geom data structure
    :param gra: A graph describing connectivity, defaults to None
    :type gra: automol graph data structure, optional
    :returns: The symmetry factor for the terminal atoms
    :rtype: float
    """
    gra = geom.graph(geo, stereo=False) if gra is None else gra
    gra = graph.implicit(gra)

    term_nkey_dct = graph.terminal_atom_neighbors(gra, atom=False)

    bords_dct = graph.kekules_bond_orders_collated(gra)
    nhyd_dct = graph.atom_implicit_hydrogens(gra)

    factor = 1.0
    for term_key, term_nkey in term_nkey_dct.items():
        # If this can only be single-bonded, it is a rotational end-group For TS graphs,
        # this will exclude reacting bonds
        bkey = frozenset({term_key, term_nkey})
        is_single = set(bords_dct[bkey]) == {1}

        nhyd = nhyd_dct[term_key]
        if nhyd and is_single:
            factor *= nhyd

    return factor


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
