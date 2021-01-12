""" construct transition state z-matrices
"""

import functools
import numpy
import automol
from automol.graph._graph import atom_neighbor_keys as _atom_neighbor_keys
from automol.zmatrix._zmatrix import shift_row_to_end
from automol.zmatrix._util import shifted_standard_zmas_graphs
from automol.zmatrix._util import reorder_zmatrix_for_redef
from automol.zmatrix._util import shift_vals_from_dummy
from automol.zmatrix import _const as constbuild
import automol.graph.trans_old as trans_old


def hydrogen_migration(rct_zmas_lst, prd_zmas_lst, tras):
    """ z-matrix for a hydrogen migration reaction

        tras in the ZMATRIX coordinate system
    """

    # Find the zma and traj pair that has min dist among rxn coord
    init_zma, min_tra = min_dist(rct_zmas_lst, tras)
    _ = prd_zmas_lst  # unneeded for TS build
    frm_bnd_keys, = trans_old.formed_bond_keys(min_tra)
    brk_bnd_keys, = trans_old.broken_bond_keys(min_tra)

    # Build a properly ordered zma
    count = 0
    while True:

        # figure out which idx in frm_bnd_keys corresponds to the hydrogen
        symbols = automol.vmatrix.symbols(automol.zmatrix.var_(init_zma))
        dist_coo_key = tuple(reversed(sorted(frm_bnd_keys)))
        for idx in dist_coo_key:
            if symbols[idx] == 'H':
                h_idx = idx
            else:
                a1_idx = idx

        brk_dist_coo_key = tuple(reversed(sorted(brk_bnd_keys)))
        for idx in brk_dist_coo_key:
            if symbols[idx] != 'H':
                a2_idx = idx

        # determine if the zmatrix needs to be rebuilt by x2z
        # determines if the hydrogen atom is used to define other atoms
        init_keys = automol.zmatrix.key_matrix(init_zma)
        rebuild = False
        for i, _ in enumerate(init_keys):
            if i > h_idx and any(idx == h_idx for idx in init_keys[i]):
                rebuild = True
                break

        # rebuild zmat and go through while loop again if needed
        # shifts order of cartesian coords and rerun x2z to get a new zmat
        # else go to next stage
        if rebuild:
            re_zma, frm_bnd_keys, brk_bnd_keys = reorder_zmatrix_for_redef(
                init_zma, a1_idx, h_idx, frm_bnd_keys, brk_bnd_keys)
            rct_zma = [re_zma]
            count += 1
            if count == 6:
                break
        else:
            rct_zma = init_zma
            break

    if not rct_zma:
        return None

    # Determine the backbone atoms to redefine the z-matrix entry
    _, gras = shifted_standard_zmas_graphs([rct_zma], remove_stereo=True)
    gra = functools.reduce(automol.graph.union, gras)
    xgr1, = automol.graph.connected_components(gra)
    chains_dct = automol.graph.atom_longest_chains(xgr1)

    idx_found = True
    a3_idx = chains_dct[a2_idx][1]
    if a3_idx in (h_idx, a1_idx):
        idx_found = False
        a2_neighbors = _atom_neighbor_keys(xgr1)[a2_idx]
        for idx in a2_neighbors:
            if idx not in (h_idx, a1_idx):
                a3_idx = idx
                idx_found = True

    if not idx_found:
        a3_idx = chains_dct[a1_idx][1]
        if a3_idx in (h_idx, a2_idx):
            a1_neighbors = _atom_neighbor_keys(xgr1)[a1_idx]
            for idx in a1_neighbors:
                if idx not in (h_idx, a2_idx):
                    a3_idx = idx

    # Set the values of the coordinates of the migrating H atom
    rct_geo = automol.zmatrix.geometry(rct_zma)
    distance = automol.geom.distance(
        rct_geo, h_idx, a1_idx)
    angle = automol.geom.central_angle(
        rct_geo, h_idx, a1_idx, a2_idx)
    dihedral = automol.geom.dihedral_angle(
        rct_geo, h_idx, a1_idx, a2_idx, a3_idx)

    dihed_is_180 = numpy.isclose(
        abs(dihedral), (numpy.pi), rtol=(5.0 * numpy.pi / 180.))
    dihed_is_0 = numpy.isclose(
        dihedral, (0.0), rtol=(5.0 * numpy.pi / 180.))
    dihed_is_360 = numpy.isclose(
        abs(dihedral), (2*numpy.pi), rtol=(5.0 * numpy.pi / 180.))

    if dihed_is_180 or dihed_is_0 or dihed_is_360:
        dihedral -= (15.0 * numpy.pi / 180.)

    # Reset the keys for the migrating H atom
    new_idxs = (a1_idx, a2_idx, a3_idx)
    key_dct = {h_idx: new_idxs}
    ts_zma = automol.zmatrix.set_keys(rct_zma, key_dct)
    h_names = automol.zmatrix.name_matrix(ts_zma)[h_idx]
    ts_zma = automol.zmatrix.set_values(
        ts_zma, {h_names[0]: distance,
                 h_names[1]: angle,
                 h_names[2]: dihedral}
    )

    # standardize the ts zmat and get tors and dist coords
    coo_dct = automol.zmatrix.coordinates(ts_zma)
    for coo_name, coo_key in coo_dct.items():
        dist_coo_key_rev = dist_coo_key[::-1]
        if coo_key[0] in (dist_coo_key, dist_coo_key_rev):
            dist_name = coo_name
            break

    # if migrating H atom is not the final zmat entry, shift it to the end (if
    # needed)
    if h_idx != automol.zmatrix.count(ts_zma) - 1:
        ts_zma, h_idx, frm_bnd_keys, brk_bnd_keys = shift_row_to_end(
            ts_zma, h_idx, frm_bnd_keys, brk_bnd_keys)

    ts_name_dct = automol.zmatrix.standard_names(ts_zma)
    dist_name = ts_name_dct[dist_name]
    ts_zma = automol.zmatrix.standard_form(ts_zma)

    # Build the torsional coordinates
    const_bnd_keys, _, _ = constbuild.hydrogen_migration(
        ts_zma, rct_zma, h_idx,
        frm_bnd_keys, brk_bnd_keys)
    # tors_bnd_keys = [coord_idxs(ts_zma, name) for name in tors_names]

    # Build the bond coordinates
    frm_bnd_keys = shift_vals_from_dummy(frm_bnd_keys, ts_zma)
    brk_bnd_keys = shift_vals_from_dummy(brk_bnd_keys, ts_zma)
    const_bnd_keys = shift_vals_from_dummy(const_bnd_keys, ts_zma)

    # Set the return
    ret = {
        'ts_zma': ts_zma,
        'bnd_keys': (frm_bnd_keys, brk_bnd_keys),
        'const_keys': (const_bnd_keys, frozenset({}), frozenset({}))
    }

    return ret


def concerted_unimol_elimination(rct_zmas, prd_zmas, tras):
    """ z-matrix for a concerted unimolecular elimination reaction
    """

    # Find the zma and traj pair that has min dist among rxn coord
    _ = prd_zmas  # Unneeded
    rct_zma, min_tra = min_dist(rct_zmas, tras)
    frm_bnd_keys, = trans_old.formed_bond_keys(min_tra)
    brk_bnd_keys = trans_old.broken_bond_keys(min_tra)
    brk_bnd_key1, brk_bnd_key2 = brk_bnd_keys

    count = 1
    while True:

        init_zma = rct_zma
        print('init_zma\n', automol.zmatrix.string(init_zma))
        # print('init_geo\n', automol.geom.string(automol.zmatrix.geometry(init_zma)))

        # Get index for migrating atom (or bond-form atom in group)
        brk_bnd_key1, brk_bnd_key2 = brk_bnd_keys
        for bnd_key in (brk_bnd_key1, brk_bnd_key2):
            if bnd_key & frm_bnd_keys:
                mig_key = next(iter(bnd_key & frm_bnd_keys))
        for key in frm_bnd_keys:
            if key != mig_key:
                a1_idx = key

        # Get chain for redefining the rc1_atm1_key z-matrix entries
        _, gras = shifted_standard_zmas_graphs(
            [init_zma], remove_stereo=True)
        gra = functools.reduce(automol.graph.union, gras)
        xgr1, = automol.graph.connected_components(gra)
        atm1_neighbors = _atom_neighbor_keys(xgr1)[a1_idx]
        for idx in atm1_neighbors:
            num_keys = len(_atom_neighbor_keys(xgr1)[idx])
            if idx != mig_key and num_keys > 1:
                a2_idx = idx
        atm2_neighbors = _atom_neighbor_keys(xgr1)[a2_idx]
        for idx in atm2_neighbors:
            if idx not in (mig_key, a1_idx):
                a3_idx = idx

        mig_redef_keys = (a1_idx, a2_idx, a3_idx)
        print('mig key', mig_key)
        print('mig redef keys', mig_redef_keys)
        print('frm keys', frm_bnd_keys)
        print('brk keys', brk_bnd_keys)
        # determine if the zmatrix needs to be rebuilt by x2z
        # determines if the hydrogen atom is used to define other atoms
        rebuild = False
        if any(idx > mig_key for idx in mig_redef_keys):
            rebuild = True

        # rebuild zmat and go through while loop again if needed
        # shift order of cartesian coords & rerun x2z to get a new zmat
        # else go to next stage
        if rebuild:
            # multiple bond keys need to be passed
            reord_zma, frm_bnd_keys, brk_bnd_keys = reorder_zmatrix_for_redef(
                init_zma, a1_idx, mig_key, frm_bnd_keys, brk_bnd_keys)
            rct_zma = reord_zma
            count += 1
            if count == 3:
                finish_build = False
                break
            print('rebuild')
        else:
            rct_zma = init_zma
            finish_build = True
            print('no rebuild')
            break

    # If z-mat with good order not found, exit function
    if not finish_build:
        return None

    # determine the new coordinates
    rct_geo = automol.zmatrix.geometry(rct_zma)
    distance = automol.geom.distance(
        rct_geo, mig_key, a1_idx)
    angle = automol.geom.central_angle(
        rct_geo, mig_key, a1_idx, a2_idx)
    dihedral = automol.geom.dihedral_angle(
        rct_geo, mig_key, a1_idx, a2_idx, a3_idx)

    # Reset the keys for the migrating H atom
    new_idxs = (a1_idx, a2_idx, a3_idx)
    key_dct = {mig_key: new_idxs}
    ts_zma = automol.zmatrix.set_keys(rct_zma, key_dct)

    # Reset the values in the value dict
    mig_names = automol.zmatrix.name_matrix(ts_zma)[mig_key]
    ts_zma = automol.zmatrix.set_values(
        ts_zma, {mig_names[0]: distance,
                 mig_names[1]: angle,
                 mig_names[2]: dihedral}
    )

    # standardize the ts zmat and get tors and dist coords
    coo_dct = automol.zmatrix.coordinates(ts_zma)
    dist_coo_key = tuple(reversed(sorted(frm_bnd_keys)))
    dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
                     if dist_coo_key in coo_keys)
    ts_name_dct = automol.zmatrix.standard_names(ts_zma)
    dist_name = ts_name_dct[dist_name]
    ts_zma = automol.zmatrix.standard_form(ts_zma)

    brk_bnd_keys = frozenset({
        shift_vals_from_dummy(brk_bnd_key1, ts_zma),
        shift_vals_from_dummy(brk_bnd_key2, ts_zma)})
    frm_bnd_keys = shift_vals_from_dummy(frm_bnd_keys, ts_zma)

    # Set the return
    ret = {
        'ts_zma': ts_zma,
        'bnd_keys': (frm_bnd_keys, brk_bnd_keys),
        'const_keys': (frozenset({}), frozenset({}), frozenset({}))
    }

    return ret


def ring_forming_scission(rct_zmas, prd_zmas, tras):
    """zmatrix for a ring forming bond scission
    """

    _ = prd_zmas  # Unneeded
    tra = tras[0]
    brk_bnd_keys, = trans_old.broken_bond_keys(tra)
    frm_bnd_keys, = trans_old.formed_bond_keys(tra)
    ts_zma = rct_zmas[0]

    # set up radical atom, leaving atom, newly formed radical atom
    # also set up chain between radical atom and newly formed radical atom
    ts_gra = automol.zmatrix.graph(ts_zma)
    rad_atm = list(automol.graph.sing_res_dom_radical_atom_keys(ts_gra))[0]
    for atm in brk_bnd_keys:
        if atm not in frm_bnd_keys:
            leave_atm = atm
        else:
            new_rad_atm = atm

    chain_between = automol.zmatrix.chain_between(ts_zma, new_rad_atm, rad_atm)

    tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)
    coo_dct = automol.zmatrix.coordinates(ts_zma)
    ang_90 = numpy.pi/2.
    ts_tors_names = []
    const_tors_names = []

    # set torsion from rad atom towards chain to 90
    for tors_name in tors_names:
        axis = coo_dct[tors_name][0][1:3]
        # (ii) remove torsions in chain_between from final torsion sampling lst
        if ((axis[0] not in chain_between) or (axis[1] not in chain_between)):
            ts_tors_names.append(tors_name)
        if ((rad_atm == axis[0] and axis[1] in chain_between) or
                (rad_atm == axis[1] and axis[0] in chain_between)):
            ts_zma_p = automol.zmatrix.set_values(ts_zma, {tors_name: ang_90})

    # vary torsion in chain_between to minimize dist b/w rad_atm & new_rad_atm
    preopt_tors_names = []
    for tors_name in tors_names:
        axis = coo_dct[tors_name][0][1:3]
        if ((axis[0] in chain_between) and (axis[1] in chain_between) and
                (rad_atm not in axis) and (new_rad_atm not in axis)):
            preopt_tors_names.append(tors_name)
            # add ring-form torsions to constraints to ensure 0 diheds for ring
            const_tors_names.append(tors_name)

    angles = [0., 2.*numpy.pi/3, 4.*numpy.pi/3]
    trial_zmas = [ts_zma_p]
    for preopt_tors_name in preopt_tors_names:
        new_trial_zmas = []
        for zma_i in trial_zmas:
            for ang in angles:
                new_trial_zmas.append(
                    automol.zmatrix.set_values(
                        zma_i, {preopt_tors_name: ang}))
        trial_zmas = new_trial_zmas

    dist_min = 1.0e30
    for trial_zma in trial_zmas:
        geo_i = automol.zmatrix.geometry(trial_zma)
        dist = automol.geom.distance(geo_i, rad_atm, new_rad_atm)
        if dist < dist_min:
            dist_min = dist
            ts_zma = trial_zma

    ang_stp = 2.*numpy.pi/6.
    # vary torsion from new_rad_atm to leaving atom so that leave_atm is far from rad_atm
    for tors_name in tors_names:
        ang = -ang_stp
        axis = coo_dct[tors_name][0][1:3]
        if ((new_rad_atm == axis[0] and axis[1] in chain_between) or
                (new_rad_atm == axis[1] and axis[0] in chain_between)):
            dist_max = 0.0
            for _ in range(6):
                ang += ang_stp
                ts_zma_i = automol.zmatrix.set_values(ts_zma, {tors_name: ang})
                geo_i = automol.zmatrix.geometry(ts_zma_i)
                dist = automol.geom.distance(geo_i, rad_atm, leave_atm)
                if dist > dist_max:
                    dist_max = dist
                    ts_zma_max = ts_zma_i
            const_tors_names.append(tors_name)
        # Set up ts torsions: remove ones with axis in chain b/w new & old rad atoms
        if ((axis[0] not in chain_between) or (axis[1] not in chain_between)):
            ts_tors_names.append(tors_name)

    # Set TS ZMA to one with max dist
    ts_zma = ts_zma_max

    # vary angles to decrease rad_atm to new_rad_atm to < 2.25 Ang
    dist_thresh = 4.25
    ang_names = automol.zmatrix.central_angle_names(ts_zma)
    ring_angs = []
    const_angs_names = []
    for ang_name in ang_names:
        ang_atms = coo_dct[ang_name][0]
        if ((ang_atms[0] in chain_between) and
               (ang_atms[1] in chain_between) and
               (ang_atms[2] in chain_between)):
            ring_angs.append(ang_name)
            const_angs_names.append(ang_name)
    dist = 1.e30
    ang_stp = numpy.pi/360.
    counter = 0
    while ((dist > dist_thresh) and (counter < 30)):
        counter += 1
        values = automol.zmatrix.values(ts_zma)
        for ang_name in ring_angs:
            ang = values[ang_name] - ang_stp
            ts_zma = automol.zmatrix.set_values(ts_zma, {ang_name: ang})
        geo_i = automol.zmatrix.geometry(ts_zma)
        dist = automol.geom.distance(geo_i, rad_atm, new_rad_atm)

    ret = {
        'ts_zma': ts_zma,
        'bnd_keys': (frozenset({}), brk_bnd_keys),
        'const_keys': (frozenset({}), frozenset({}), frozenset({}))
    }

    return ret


def beta_scission(rct_zmas, prd_zmas, tras):
    """ z-matrix for a beta-scission reaction
    """

    # Obtain frm and brk keys from arbtrary tras
    tra = tras[0]
    brk_bnd_keys, = trans_old.broken_bond_keys(tra)
    _ = prd_zmas  # Unneeded

    # Obtain the TS ZMA
    ts_zma, = rct_zmas  # wrong since rct_zmas is a lst
    ts_zma = automol.zmatrix.standard_form(ts_zma)
    brk_bnd_keys = shift_vals_from_dummy(brk_bnd_keys, ts_zma)

    ret = {
        'ts_zma': ts_zma,
        'bnd_keys': (frozenset({}), brk_bnd_keys),
        'const_keys': (frozenset({}), frozenset({}), frozenset({}))
    }

    return ret


#  Funtions to select a certain ZMA for the build
def min_dist(rct_zmas_lst, tras, dthresh=1000.0):
    """ determine a structure that has minimum distance for a forming bond
    """

    min_tra = None
    for tra in tras:
        for rct_zmas in rct_zmas_lst:
            # Calculate distance for rct1 (assuming unimol) along forming bond
            frm_bnd_keys, = trans_old.formed_bond_keys(tra)
            zma = rct_zmas[0]  # assumes unimolecular list, may need change?
            geo = automol.zmatrix.geometry(zma)
            dist = automol.geom.distance(geo, *list(frm_bnd_keys))
            # Determine if distance is a new minimum
            if dist < dthresh:
                dthresh = dist
                min_tra = tra
                min_zma = zma

    return min_zma, min_tra


####
# geo_ret = {
#     'rcts_gra': rcts_gra,
#     'bnd_keys': (frozenset({}), brk_bnd_keys),
# }
