""" Deal with comparisons
"""

from phydat import bnd

from .. import geom, graph, zmat
from ..util import dict_
from ._0core import class_, ts_graph


def similar_saddle_point_structure(zma, ref_zma, zrxn, sens=1.0):
    """Perform a series of checks to assess the viability
    of a transition state geometry prior to saving
    """

    # Initialize viable
    viable = True

    # Get the bond dists and calculate the distance of bond being formed
    ref_geo = zmat.geometry(ref_zma, dummy=True)
    cnf_geo = zmat.geometry(zma, dummy=True)
    ts_gra = ts_graph(zrxn)

    frm_bnd_keys = graph.ts.forming_bond_keys(ts_gra)
    brk_bnd_keys = graph.ts.breaking_bond_keys(ts_gra)
    tra_dct = graph.ts.transferring_atoms(ts_gra)

    cnf_dist_lst = []
    ref_dist_lst = []
    bnd_key_lst = []
    cnf_ang_lst = []
    ref_ang_lst = []
    for frm_bnd_key in frm_bnd_keys:
        frm_idx1, frm_idx2 = frm_bnd_key
        cnf_dist = geom.distance(cnf_geo, frm_idx1, frm_idx2)
        ref_dist = geom.distance(ref_geo, frm_idx1, frm_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(frm_bnd_key)

    for brk_bnd_key in brk_bnd_keys:
        brk_idx1, brk_idx2 = brk_bnd_key
        cnf_dist = geom.distance(cnf_geo, brk_idx1, brk_idx2)
        ref_dist = geom.distance(ref_geo, brk_idx1, brk_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(brk_bnd_key)

    for tra_key, (don_key, acc_key) in tra_dct.items():
        cnf_ang = geom.central_angle(cnf_geo, don_key, tra_key, acc_key)
        ref_ang = geom.central_angle(ref_geo, don_key, tra_key, acc_key)
        cnf_ang_lst.append(cnf_ang)
        ref_ang_lst.append(ref_ang)

    # Set the maximum allowed displacement for a TS conformer
    max_disp = 0.6 * sens
    # better to check for bond-form length in bond scission with ring forming
    if "addition" in class_(zrxn):
        max_disp = 0.8 * sens
    if "abstraction" in class_(zrxn):
        # this was 1.4 - SJK reduced it to work for some OH abstractions
        max_disp = 1.0 * sens

    # Check forming bond angle similar to ini config
    if "elimination" not in class_(zrxn):
        for ref_angle, cnf_angle in zip(ref_ang_lst, cnf_ang_lst):
            if abs(cnf_angle - ref_angle) > 0.44 * sens:
                print("transitioning bond angle has diverged", ref_angle, cnf_angle)
                viable = False

    symbs = geom.symbols(cnf_geo)
    lst_info = zip(ref_dist_lst, cnf_dist_lst, bnd_key_lst)
    for ref_dist, cnf_dist, bnd_key in lst_info:
        if "add" in class_(zrxn) or "abst" in class_(zrxn):
            bnd_key1, bnd_key2 = min(list(bnd_key)), max(list(bnd_key))
            symb1 = symbs[bnd_key1]
            symb2 = symbs[bnd_key2]

            if bnd_key in frm_bnd_keys:
                # Check if radical atom is closer to some atom
                # other than the bonding atom
                is_ok = geom.could_be_forming_bond(cnf_geo, *bnd_key, gra=ts_gra)
                if not is_ok:
                    # ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                    print(" - Radical atom now has a new nearest neighbor")
                    viable = False
                # check forming bond distance
                if abs(cnf_dist - ref_dist) > max_disp:
                    print(
                        "distance of transitioning bonds has diverged",
                        ref_dist,
                        cnf_dist,
                    )
                    viable = False

            # Check distance relative to equi. bond
            equi_bnd = dict_.value_by_unordered_key(
                bnd.LEN_DCT, (symb1, symb2), fill_val=0.0
            )
            displace_from_equi = cnf_dist - equi_bnd
            dchk1 = abs(cnf_dist - ref_dist) > 0.1 * sens
            dchk2 = displace_from_equi < 0.2 / sens
            if dchk1 and dchk2:
                print("transitioning bond differs from equi", cnf_dist, equi_bnd)
                viable = False
        else:
            # check forming/breaking bond distance
            # if abs(cnf_dist - ref_dist) > 0.4:
            # max disp of 0.4 causes problems for bond scission w/ ring forming
            # not sure if setting it to 0.3 will cause problems for other cases
            if abs(cnf_dist - ref_dist) > 0.3 * sens:
                print(
                    "transitioning bond distance has diverged from starting guess",
                    ref_dist,
                    cnf_dist,
                )
                viable = False

    if not _check_stereo_parities(zma, ref_zma, zrxn):
        viable = False

    return viable


def _check_stereo_parities(zma, ref_zma, zrxn):
    """make sure stereo is consistent with ref_zma"""
    gra = graph.without_stereo(ts_graph(zrxn))
    ste_keys = graph.unassigned_stereocenter_keys(gra, atom=False)
    geo = zmat.geometry(zma, dummy=True)
    ref_geo = zmat.geometry(ref_zma, dummy=True)

    lin_keys = set(graph.linear_atom_keys(gra))
    bad_keys = graph.geometries_parity_mismatches(gra, geo, ref_geo, ste_keys)

    viable = True
    for key in bad_keys:
        if isinstance(key, int):
            viable = False
            print("Invalid stereo at atom ", key)
        elif isinstance(key, frozenset) and not key <= lin_keys:
            viable = False
            print("Invalid stereo at bond ", key)

    return viable
