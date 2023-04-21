""" Deal with comparisons
"""

from phydat import bnd
import automol.zmat
from automol.geom import symbols
from automol.geom import distance
from automol.geom import central_angle
from automol.util import dict_
from automol.reac._reac import forming_bond_keys
from automol.reac._reac import breaking_bond_keys
from automol.reac._reac import relabel_for_geometry
from automol.graph.base._core import atom_keys
from automol.graph.base._core import without_stereo_parities
from automol.graph.base._core import from_ts_graph
from automol.graph.base._canon import stereogenic_bond_keys
from automol.graph.base._canon import parity_evaluator_from_geometry_


def similar_saddle_point_structure(zma, ref_zma, zrxn, sens=1.0):
    """ Perform a series of checks to assess the viability
        of a transition state geometry prior to saving
    """

    # Initialize viable
    viable = True

    # Get the bond dists and calculate the distance of bond being formed
    ref_geo = automol.zmat.geometry(ref_zma)
    cnf_geo = automol.zmat.geometry(zma)
    grxn = relabel_for_geometry(zrxn)

    frm_bnd_keys = forming_bond_keys(grxn)
    brk_bnd_keys = breaking_bond_keys(grxn)

    cnf_dist_lst = []
    ref_dist_lst = []
    bnd_key_lst = []
    cnf_ang_lst = []
    ref_ang_lst = []
    for frm_bnd_key in frm_bnd_keys:
        frm_idx1, frm_idx2 = list(frm_bnd_key)
        cnf_dist = distance(cnf_geo, frm_idx1, frm_idx2)
        ref_dist = distance(ref_geo, frm_idx1, frm_idx2)
        cnf_dist_lst.append(cnf_dist)
        ref_dist_lst.append(ref_dist)
        bnd_key_lst.append(frm_bnd_key)

    for brk_bnd_key in brk_bnd_keys:
        brk_idx1, brk_idx2 = list(brk_bnd_key)
        cnf_dist = distance(cnf_geo, brk_idx1, brk_idx2)
        ref_dist = distance(ref_geo, brk_idx1, brk_idx2)
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
                        cnf_ang = central_angle(cnf_geo, idx1, idx2, idx3)
                        ref_ang = central_angle(ref_geo, idx1, idx2, idx3)
                        cnf_ang_lst.append(cnf_ang)
                        ref_ang_lst.append(ref_ang)

    # Set the maximum allowed displacement for a TS conformer
    max_disp = 0.6 * sens
    # better to check for bond-form length in bond scission with ring forming
    if 'addition' in grxn.class_:
        max_disp = 0.8 * sens
    if 'abstraction' in grxn.class_:
        # this was 1.4 - SJK reduced it to work for some OH abstractions
        max_disp = 1.0 * sens

    # Check forming bond angle similar to ini config
    if 'elimination' not in grxn.class_:
        for ref_angle, cnf_angle in zip(ref_ang_lst, cnf_ang_lst):
            if abs(cnf_angle - ref_angle) > .44 * sens:
                # ioprinter.diverged_ts('angle', ref_angle, cnf_angle)
                viable = False

    symbs = symbols(cnf_geo)
    lst_info = zip(ref_dist_lst, cnf_dist_lst, bnd_key_lst)
    for ref_dist, cnf_dist, bnd_key in lst_info:
        if 'add' in grxn.class_ or 'abst' in grxn.class_:
            bnd_key1, bnd_key2 = min(list(bnd_key)), max(list(bnd_key))
            symb1 = symbs[bnd_key1]
            symb2 = symbs[bnd_key2]

            if bnd_key in frm_bnd_keys:
                # Check if radical atom is closer to some atom
                # other than the bonding atom
                cls = automol.zmat.is_atom_closest_to_bond_atom(
                    zma, bnd_key2, cnf_dist)
                if not cls:
                    # ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                    print(
                        ' - Radical atom now has a new nearest neighbor')
                    viable = False
                # check forming bond distance
                if abs(cnf_dist - ref_dist) > max_disp:
                    # ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                    viable = False

            # Check distance relative to equi. bond
            equi_bnd = dict_.values_by_unordered_tuple(
                bnd.LEN_DCT, (symb1, symb2), fill_val=0.0)
            displace_from_equi = cnf_dist - equi_bnd
            dchk1 = abs(cnf_dist - ref_dist) > 0.1 * sens
            dchk2 = displace_from_equi < 0.2 * sens
            if dchk1 and dchk2:
                # ioprinter.bad_equil_ts(cnf_dist, equi_bnd)
                viable = False
        else:
            # check forming/breaking bond distance
            # if abs(cnf_dist - ref_dist) > 0.4:
            # max disp of 0.4 causes problems for bond scission w/ ring forming
            # not sure if setting it to 0.3 will cause problems for other cases
            if abs(cnf_dist - ref_dist) > 0.3 * sens:
                # ioprinter.diverged_ts('distance', ref_dist, cnf_dist)
                viable = False

    if not _check_stereo_parities(zma, ref_zma, zrxn):
        viable = False

    return viable


def _check_stereo_parities(zma, ref_zma, zrxn):
    """make sure stereo is consistent with ref_zma
    """
    fgra = without_stereo_parities(from_ts_graph(zrxn.forward_ts_graph))
    bgra = without_stereo_parities(from_ts_graph(zrxn.backward_ts_graph))
    forw_ste_keys = stereogenic_bond_keys(fgra)
    back_ste_keys = stereogenic_bond_keys(bgra)
    forw_idxs = zrxn.key_map(rev=False, stereo=False)
    back_idxs = zrxn.key_map(rev=True, stereo=False)
    fatm_keys = atom_keys(fgra)
    batm_keys = atom_keys(bgra)
    fpri_dct = dict(zip(fatm_keys, fatm_keys))
    bpri_dct = dict(zip(batm_keys, batm_keys))
    geo = automol.zmat.geometry(zma, dummy=True)
    ref_geo = automol.zmat.geometry(ref_zma, dummy=True)

    forw_par_eval_ = parity_evaluator_from_geometry_(
        fgra, geo, geo_idx_dct=forw_idxs)
    ref_forw_par_eval_ = parity_evaluator_from_geometry_(
        fgra, ref_geo, geo_idx_dct=forw_idxs)
    back_par_eval_ = parity_evaluator_from_geometry_(
        bgra, geo, geo_idx_dct=back_idxs)
    ref_back_par_eval_ = parity_evaluator_from_geometry_(
        bgra, ref_geo, geo_idx_dct=back_idxs)

    forw_par_ = forw_par_eval_(fpri_dct)
    ref_forw_par_ = ref_forw_par_eval_(fpri_dct)
    back_par_ = back_par_eval_(bpri_dct)
    ref_back_par_ = ref_back_par_eval_(bpri_dct)

    viable = True
    for bnd_key in forw_ste_keys:
        if forw_par_(bnd_key) != ref_forw_par_(bnd_key):
            viable = False
            print('Invalid stereo at bond ', bnd_key)
    for bnd_key in back_ste_keys:
        if back_par_(bnd_key) != ref_back_par_(bnd_key):
            viable = False
            print('Invalid stereo at backward bond ', bnd_key)

    return viable
