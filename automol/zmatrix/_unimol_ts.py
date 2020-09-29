""" construct transition state z-matrices
"""

import functools
import numpy
from qcelemental import constants as qcc
import automol
from automol.graph._graph import atom_neighbor_keys as _atom_neighbor_keys
from automol.zmatrix._util import shifted_standard_zmas_graphs
from automol.zmatrix._util import reorder_zmatrix_for_migration


ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')


def min_hyd_mig_dist(rct_zmas, prd_zmas):
    """ determines distance coordinate to minimize for hydrogen migration reaction
    """
    min_frm_bnd_key = None
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    if len(rct_zmas) == 1:
        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)
        tras, _, _ = automol.graph.reac.hydrogen_migration(rct_gras, prd_gras)
        # If reaction found, then proceed
        if tras:
            min_key = None
            min_dist = 100
            for tra in tras:
                frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
                geo = automol.zmatrix.geometry(rct_zmas[0])
                dist = automol.geom.distance(geo, *list(frm_bnd_key))
                if dist < min_dist:
                    min_dist = dist
                    min_key = frm_bnd_key
                if min_key:
                    min_frm_bnd_key = min_key

    return min_frm_bnd_key


def hydrogen_migration(rct_zmas, prd_zmas):
    """ z-matrix for a hydrogen migration reaction
    """
    ret = None, None, None, None, None
    rct_zma = None

    # set products which will be unchanged to ts algorithm
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    count = 1
    while True:
        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)

        # try an atom migration then try a proton migration
        tras, _, _ = automol.graph.reac.hydrogen_migration(rct_gras, prd_gras)

        # If reaction found, then proceed
        if tras:
            # Get the bond formation keys and the reactant zmatrix
            min_dist = 100.
            frm_bnd_key = None
            brk_bnd_key = None
            for tra_i in tras:
                # Get the bond formation and breaking keys
                bnd_key, = automol.graph.trans.formed_bond_keys(tra_i)
                geo = automol.zmatrix.geometry(rct_zmas[0])
                dist = automol.geom.distance(geo, *list(bnd_key))
                if dist < min_dist:
                    min_dist = dist
                    frm_bnd_key = bnd_key
                    brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra_i)
            init_zma, = rct_zmas
            # print('init_zma test:', init_zma)

            # figure out which idx in frm_bnd_keys corresponds to the hydrogen
            symbols = automol.vmatrix.symbols(automol.zmatrix.var_(init_zma))
            dist_coo_key = tuple(reversed(sorted(frm_bnd_key)))
            for idx in dist_coo_key:
                if symbols[idx] == 'H':
                    h_idx = idx
                else:
                    a1_idx = idx

            # print('h_idx test:', h_idx, a1_idx)

            brk_dist_coo_key = tuple(reversed(sorted(brk_bnd_key)))
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
                reord_zma = reorder_zmatrix_for_migration(
                    init_zma, a1_idx, h_idx)
                rct_zmas = [reord_zma]
                count += 1
                if count == 6:
                    break
            else:
                rct_zma = init_zma
                break
        else:
            return None

    if not rct_zma:
        return None

    # determine the backbone atoms to redefine the z-matrix entry
    _, gras = shifted_standard_zmas_graphs([rct_zma], remove_stereo=True)
    gra = functools.reduce(automol.graph.union, gras)
    xgr1, = automol.graph.connected_components(gra)
    chains_dct = automol.graph.atom_longest_chains(xgr1)

    tors_names = automol.zmatrix.torsion_coordinate_names(rct_zma)

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

    # determine the new coordinates
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
    # Reset the values in the value dict
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
        ts_zma, new_h_idx, frm_bnd_key, brk_bnd_key = automol.zmatrix.shift_row_to_end(
                ts_zma, h_idx, frm_bnd_key, brk_bnd_key)
    else:
        new_h_idx = h_idx

    ts_name_dct = automol.zmatrix.standard_names(ts_zma)
    dist_name = ts_name_dct[dist_name]
    ts_zma = automol.zmatrix.standard_form(ts_zma)

    pot_tors_names = []
    for tors_name in tors_names:
        pot_tors_names.append(ts_name_dct[tors_name])

    # get full set of potential torsional coordinates

    # remove the torsional coordinates that are in the ring
    chn_end = new_h_idx
    for key in frm_bnd_key:
        if key != new_h_idx:
            chn_end = key
    ring_atoms1 = automol.zmatrix.chain_between(ts_zma, chn_end, new_h_idx)

    # Use this ring information to constrain the beta bond from the attacking atom
    const_bnd_key = None
    if len(ring_atoms1) > 3:
        const_bnd_key = frozenset({ring_atoms1[1], ring_atoms1[2]})

    chn_end = new_h_idx
    for key in brk_bnd_key:
        if key != new_h_idx:
            chn_end = key
    ring_atoms2 = automol.zmatrix.chain_between(ts_zma, chn_end, new_h_idx)
    ring_atoms = ring_atoms1
    if len(ring_atoms2) > len(ring_atoms1):
        ring_atoms = ring_atoms2

    tors_names = []
    coo_dct = automol.zmatrix.coordinates(ts_zma)
    for tors_name in pot_tors_names:
        axis = coo_dct[tors_name][0][1:3]
        if axis[0] not in ring_atoms or axis[1] not in ring_atoms:
            tors_names.append(tors_name)

    # Build the reacts graph
    if h_idx != new_h_idx:
        atom_idxs = (x for x in automol.graph.atom_keys(rct_gras[0]))
        atm_key_dct = {}
        for idx in atom_idxs:
            if idx != h_idx:
                if idx > h_idx:
                    new_idx = idx - 1
                else:
                    new_idx = idx
                atm_key_dct[idx] = new_idx
        atm_key_dct[h_idx] = new_h_idx

        rcts_gra = automol.graph.relabel(rct_gras[0], atm_key_dct)
    else:
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

    frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
    brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)
    if const_bnd_key:
        const_bnd_key = shift_vals_from_dummy(const_bnd_key, ts_zma)

    ret = ts_zma, dist_name, frm_bnd_key, brk_bnd_key, const_bnd_key, tors_names, rcts_gra

    return ret


def min_unimolecular_elimination_dist(rct_zmas, prd_zmas):
    """ determines distance coordinate to minimize for a
        concerted unimolecular elimination reaction
    """
    min_frm_bnd_key = None
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    if len(rct_zmas) == 1:
        rct_zmas, rct_gras = shifted_standard_zmas_graphs(
            rct_zmas, remove_stereo=True)
        tras, _, _ = automol.graph.reac.elimination(rct_gras, prd_gras)
        if tras:
            if len(tras[0]) == 1:
                tras = [tras]
            min_frm_bnd_key = None
            min_dist = 10
            for tra in tras:
                frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
                geo = automol.zmatrix.geometry(rct_zmas[0])
                dist = automol.geom.distance(geo, *list(frm_bnd_key))
                if dist < min_dist:
                    min_dist = dist
                    min_frm_bnd_key = frm_bnd_key
    # print('AT ELIM1')
    # import sys
    # sys.exit()

    return min_frm_bnd_key


def concerted_unimolecular_elimination(rct_zmas, prd_zmas):
    """ z-matrix for a concerted unimolecular elimination reaction
    """

    # Initialize info for the returns
    ret = None, None, None, None, None
    finish_build = True

    # Attempt to build appropriate z-matrix
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    if len(rct_zmas) == 1:
        count = 1
        while True:
            rct_zmas, rct_gras = shifted_standard_zmas_graphs(
                rct_zmas, remove_stereo=True)
            init_zma, = rct_zmas

            tras, _, _ = automol.graph.reac.elimination(rct_gras, prd_gras)
            if tras is not None:
                if len(tras[0]) == 1:
                    tras = [tras]
                min_dist = 100.
                frm_bnd_key = None
                for tra_i in tras:
                    # Get the bond formation and breaking keys
                    bnd_key, = automol.graph.trans.formed_bond_keys(tra_i)
                    geo = automol.zmatrix.geometry(rct_zmas[0])
                    dist = automol.geom.distance(geo, *list(bnd_key))
                    if dist < min_dist:
                        min_dist = dist
                        frm_bnd_key = bnd_key
                        tra = tra_i
                brk_keys = automol.graph.trans.broken_bond_keys(tra)
                brk_bnd_key1, brk_bnd_key2 = brk_keys
                init_zma, = rct_zmas


                # Get index for migrating atom (or bond-form atom in group)
                for bnd_key in (brk_bnd_key1, brk_bnd_key2):
                    if bnd_key & frm_bnd_key:
                        mig_key = next(iter(bnd_key & frm_bnd_key))
                for key in frm_bnd_key:
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

                # determine if the zmatrix needs to be rebuilt by x2z
                # determines if the hydrogen atom is used to define other atoms
                rebuild = False
                if any(idx > mig_key for idx in mig_redef_keys):
                    rebuild = True

                # rebuild zmat and go through while loop again if needed
                # shift order of cartesian coords & rerun x2z to get a new zmat
                # else go to next stage
                if rebuild:
                    reord_zma = reorder_zmatrix_for_migration(
                        init_zma, a1_idx, mig_key)
                    rct_zmas = [reord_zma]
                    count += 1
                    if count == 3:
                        finish_build = False
                        break
                else:
                    rct_zma = init_zma
                    finish_build = True
                    break
            else:
                finish_build = False

    # If z-mat with good order found, finish building it
    if finish_build:

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
        dist_coo_key = tuple(reversed(sorted(frm_bnd_key)))
        dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
                         if dist_coo_key in coo_keys)
        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        dist_name = ts_name_dct[dist_name]
        ts_zma = automol.zmatrix.standard_form(ts_zma)

        # Get the name of the coordinate of the other bond that is breaking
        brk_dist_name = None
        for brk_key in (brk_bnd_key1, brk_bnd_key2):
            if not brk_key.intersection(frm_bnd_key):
                brk_dist_name = automol.zmatrix.bond_key_from_idxs(
                    ts_zma, brk_key)

        # Add second attempt to get brk_dist_name
        if brk_dist_name is None:
            brk_dist_names = [
                automol.zmatrix.bond_key_from_idxs(ts_zma, brk_bnd_key1),
                automol.zmatrix.bond_key_from_idxs(ts_zma, brk_bnd_key2)
            ]
            # Grab the name that is not None
            for name in brk_dist_names:
                if name is not None:
                    brk_dist_name = name

        # get full set of potential torsional coordinates
        pot_tors_names = automol.zmatrix.torsion_coordinate_names(rct_zma)

        # remove the torsional coordinates that would break reaction coordinate
        gra = automol.zmatrix.graph(ts_zma, remove_stereo=True)
        coo_dct = automol.zmatrix.coordinates(ts_zma)
        tors_names = []
        for tors_name in pot_tors_names:
            axis = coo_dct[tors_name][0][1:3]
            grp1 = [axis[1]] + (
                list(automol.graph.branch_atom_keys(gra, axis[0], axis) -
                     set(axis)))
            grp2 = [axis[0]] + (
                list(automol.graph.branch_atom_keys(gra, axis[1], axis) -
                     set(axis)))
            if not ((mig_key in grp1 and a1_idx in grp2) or
                    (mig_key in grp2 and a1_idx in grp1)):
                tors_names.append(tors_name)

        # Get reactants graph
        _, rct_gras = shifted_standard_zmas_graphs(
            [rct_zma], remove_stereo=True)
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

        brk_bnd_key1 = shift_vals_from_dummy(brk_bnd_key1, ts_zma)
        brk_bnd_key2 = shift_vals_from_dummy(brk_bnd_key2, ts_zma)
        brk_bnd_keys = frozenset({brk_bnd_key1, brk_bnd_key2})
        frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)

        ret = ts_zma, dist_name, brk_dist_name, brk_bnd_keys, frm_bnd_key, tors_names, rcts_gra

    return ret


def ring_forming_scission(rct_zmas, prd_zmas):
    """zmatrix for a ring forming bond scission
    """
    ret = None
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    tras, _, _ = automol.graph.reac.ring_forming_scission(rct_gras, prd_gras)
    if tras:
        tra = tras[0]
        brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)
        frm_bnd_key, = automol.graph.trans.formed_bond_keys(tra)
        ts_zma = rct_zmas[0]

        # set up radical atom, leaving atom, newly formed radical atom
        # also set up chain between radical atom and newly formed radical atom
        ts_gra = automol.zmatrix.graph(ts_zma)
        rad_atm = list(automol.graph.sing_res_dom_radical_atom_keys(ts_gra))[0]
        for atm in brk_bnd_key:
            if atm not in frm_bnd_key:
                leave_atm = atm
            else:
                new_rad_atm = atm

        chain_between = automol.zmatrix.chain_between(ts_zma, new_rad_atm, rad_atm)

        tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)
        coo_dct = automol.zmatrix.coordinates(ts_zma)
        ang_90 = numpy.pi/2.
        ts_tors_names = []
        const_tors_names = []
        # (i) set torsion from rad atom towards chain to 90
        for tors_name in tors_names:
            axis = coo_dct[tors_name][0][1:3]
            # (ii) remove torsions in chain_between from final torsion sampling list
            if ((axis[0] not in chain_between) or (axis[1] not in chain_between)):
                ts_tors_names.append(tors_name)
            if ((rad_atm == axis[0] and axis[1] in chain_between) or
                    (rad_atm == axis[1] and axis[0] in chain_between)):
                ts_zma_p = automol.zmatrix.set_values(ts_zma, {tors_name: ang_90})
                # const_tors_names.append(tors_name)

        # (iii) vary torsions in chain_between to minimize distance from rad_atm to new_rad_atm
        preopt_tors_names = []
        for tors_name in tors_names:
            axis = coo_dct[tors_name][0][1:3]
            if ((axis[0] in chain_between) and (axis[1] in chain_between) and
                    (rad_atm not in axis) and (new_rad_atm not in axis)):
                preopt_tors_names.append(tors_name)
                # add any ring forming torsions to constraints to ensure 0 dihedrals for the ring
                const_tors_names.append(tors_name)

        print('preopt_tors_names:', preopt_tors_names)

        angles = [0., 2.*numpy.pi/3, 4.*numpy.pi/3]
        # angles = [0., numpy.pi/3., 2.*numpy.pi/3, 3.*numpy.pi/3., 4.*numpy.pi/3, 5*numpy.pi/3.]
        trial_zmas = [ts_zma_p]
        for preopt_tors_name in preopt_tors_names:
            new_trial_zmas = []
            for zma_i in trial_zmas:
                for ang in angles:
                    new_trial_zmas.append(automol.zmatrix.set_values(zma_i, {preopt_tors_name: ang}))
            trial_zmas = new_trial_zmas

        dist_min = 1.0e30
        for trial_zma in trial_zmas:
            geo_i = automol.zmatrix.geometry(trial_zma)
            dist = automol.geom.distance(geo_i, rad_atm, new_rad_atm)
            if dist < dist_min:
                dist_min = dist
                ts_zma = trial_zma

        ang_stp = 2.*numpy.pi/6.
        # (iv) vary torsion from new_rad_atm to leaving atom so that leave_atm is far from rad_atm
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
            # set up ts torsions - remove ones with axis in the chain between new and old rad atoms
            if ((axis[0] not in chain_between) or (axis[1] not in chain_between)):
                ts_tors_names.append(tors_name)
            # elif (axis[0] in chain_between) and (axis[1] in chain_between):
                # if tors_name not in const_tors_names:
                    # const_tors_names.append(tors_name)

        print('ts_tors_names test:', ts_tors_names, automol.zmatrix.string(ts_zma))
        print('const_tors_names test:', const_tors_names)
        ts_zma = ts_zma_max

        # (v) vary angles to decrease rad_atm to new_rad_atm to < 2.25 Angstroms
        dist_thresh = 4.25
        # dist_thresh = 4.
        ang_names = automol.zmatrix.central_angle_names(ts_zma)
        ring_angs = []
        const_angs_names = []
        for ang_name in ang_names:
            ang_atms = coo_dct[ang_name][0]
            if ((ang_atms[0] in chain_between) and (ang_atms[1] in chain_between) and
                    (ang_atms[2] in chain_between)):
                ring_angs.append(ang_name)
                const_angs_names.append(ang_name)
        dist = 1.e30
        ang_stp = numpy.pi/360.
        # ang_stp = 0.5 degrees
        counter = 0
        while ((dist > dist_thresh) and (counter < 30)):
            counter += 1
            values = automol.zmatrix.values(ts_zma)
            for ang_name in ring_angs:
                ang = values[ang_name] - ang_stp
                ts_zma = automol.zmatrix.set_values(ts_zma, {ang_name: ang})
            geo_i = automol.zmatrix.geometry(ts_zma)
            dist = automol.geom.distance(geo_i, rad_atm, new_rad_atm)

        brk_dist_name = automol.zmatrix.bond_key_from_idxs(ts_zma, brk_bnd_key)

        # Build the reactants graph
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

        ret = ts_zma, brk_dist_name, brk_bnd_key, const_tors_names, ts_tors_names, const_angs_names, rcts_gra

    return ret


def beta_scission(rct_zmas, prd_zmas):
    """ z-matrix for a beta-scission reaction
    """
    ret = None
    print('initial rct_zmas in beta-scission test:', automol.zmatrix.string(rct_zmas[0]))
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    tras, _, _ = automol.graph.reac.beta_scission(rct_gras, prd_gras)
    if tras:
        tra = tras[0]
        brk_bnd_key, = automol.graph.trans.broken_bond_keys(tra)
        ts_zma, = rct_zmas
        coo_dct = automol.zmatrix.coordinates(ts_zma)
        dist_coo_key = tuple(reversed(sorted(brk_bnd_key)))
        dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
                         if dist_coo_key in coo_keys)

        ts_name_dct = automol.zmatrix.standard_names(ts_zma)
        dist_name = ts_name_dct[dist_name]
        print('initial ts_zma in beta-scission test:', automol.zmatrix.string(ts_zma))
        ts_zma = automol.zmatrix.standard_form(ts_zma)
        tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)

        brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)

        # Build the reactants graph
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

        ret = ts_zma, dist_name, brk_bnd_key, tors_names, rcts_gra

        print('zmas in beta-scission test:', automol.zmatrix.string(ts_zma), automol.zmatrix.string(rct_zmas[0]))

    return ret


def shift_vals_from_dummy(vals, zma):
    """ Shift a set of values using remdummy
        Shift requires indices be 1-indexed
    """
    type_ = type(vals)

    dummy_idxs = automol.zmatrix.atom_indices(zma, sym='X')

    shift_vals = []
    print('vals test:', vals)
    for val in vals:
        shift = 0
        for dummy in dummy_idxs:
            if val >= dummy:
                shift += 1
        shift_vals.append(val+shift)

    shift_vals = type_(shift_vals)

    return shift_vals
