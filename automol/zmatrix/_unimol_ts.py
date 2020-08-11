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
            min_dist = 10
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

            # figure out which idx in frm_bnd_keys corresponds to the hydrogen
            symbols = automol.vmatrix.symbols(automol.zmatrix.var_(init_zma))
            dist_coo_key = tuple(reversed(sorted(frm_bnd_key)))
            for idx in dist_coo_key:
                if symbols[idx] == 'H':
                    h_idx = idx
                else:
                    a1_idx = idx

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

    # determine the backbone atoms to redefine the z-matrix entry
    _, gras = shifted_standard_zmas_graphs([rct_zma], remove_stereo=True)
    gra = functools.reduce(automol.graph.union, gras)
    xgr1, = automol.graph.connected_components(gra)
    chains_dct = automol.graph.atom_longest_chains(xgr1)

    print('rct_zma test:', automol.zmatrix.string(rct_zma))
    tors_names = automol.zmatrix.torsion_coordinate_names(rct_zma)
    print('tors_names rct test:', tors_names)
    # print('rct_zma:\n',automol.zmatrix.string(rct_zma))
    idx_found = True
    a3_idx = chains_dct[a2_idx][1]
    if a3_idx in (h_idx, a1_idx):
        idx_found = False
        a2_neighbors = _atom_neighbor_keys(xgr1)[a2_idx]
        for idx in a2_neighbors:
            # print('idx test 1:', a3_idx, a2_neighbors)
            if idx not in (h_idx, a1_idx):
                a3_idx = idx
                idx_found = True

    if not idx_found:
        a3_idx = chains_dct[a1_idx][1]
        if a3_idx in (h_idx, a2_idx):
            a1_neighbors = _atom_neighbor_keys(xgr1)[a1_idx]
            # print('idx test 2:', a3_idx, a1_neighbors)
            for idx in a1_neighbors:
                if idx not in (h_idx, a2_idx):
                    a3_idx = idx

    # idx_found = True
    # a3_idx = chains_dct[a1_idx][1]
    # if a3_idx in (h_idx, a2_idx):
        # idx_found = False
        # a1_neighbors = _atom_neighbor_keys(xgr1)[a1_idx]
        # for idx in a1_neighbors:
            # print('idx test 1:', a3_idx, a1_neighbors)
            # if idx not in (h_idx, a2_idx):
                # a3_idx = idx
                # idx_found = True

    # if not idx_found:
        # a3_idx = chains_dct[a2_idx][1]
        # if a3_idx in (h_idx, a1_idx):
            # a2_neighbors = _atom_neighbor_keys(xgr1)[a2_idx]
            # print('idx test 2:', a3_idx, a2_neighbors)
            # for idx in a2_neighbors:
                # if idx not in (h_idx, a1_idx):
                    # a3_idx = idx

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
        ts_zma, new_h_idx, shft_tors_names = automol.zmatrix.shift_row_to_end(ts_zma, h_idx, tors_names)
    else:
        new_h_idx = h_idx
    print('shft test:', shft_tors_names)
    print('ts_zma test before standard:', automol.zmatrix.string(ts_zma))

    ts_name_dct = automol.zmatrix.standard_names(ts_zma)
    dist_name = ts_name_dct[dist_name]
    ts_zma = automol.zmatrix.standard_form(ts_zma)

    print('ts_zma test after standard:', automol.zmatrix.string(ts_zma))

    pot_tors_names = []
    #for tors_name in shft_tors_names:
    for tors_name in tors_names:
        pot_tors_names.append(ts_name_dct[tors_name])

    # get full set of potential torsional coordinates
    # tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)
    print('pot_tors_names test:', pot_tors_names)

    # remove the torsional coordinates that are in the ring
    chn_end = h_idx
    for key in frm_bnd_key:
        if key != h_idx
            chn_end = key

    ring_atoms = automol.zmatrix.chain_between(ts_zma, chn_end, h_idx)
    tors_names = []
    coo_dct = automol.zmatrix.coordinates(ts_zma)
    axis_dict = {}
    bad_tors = []
    for tors_name in pot_tors_names:
        axis = coo_dct[tors_name][0][1:3]
        if axis[0] not in ring_atoms or axis[1] not in ring_atoms:
            tors_names.append(tors_name)

    # remove the torsional coordinates that would break reaction coordinate
    # gra = automol.zmatrix.graph(ts_zma, remove_stereo=True)
    # coo_dct = automol.zmatrix.coordinates(ts_zma)
    # key_mat = automol.zmatrix.key_matrix(ts_zma)
    # h_keys = key_mat[new_h_idx]
    # a1_key = h_keys[0]
    # a2_key = h_keys[1]
    # tors_names = []
    # axis_dict = {}
    # bad_tors = []
    # for tors_name in pot_tors_names:
        # axis = coo_dct[tors_name][0][1:3]
        # if axis not in axis_dict:
            # axis_dict[axis] = [tors_name]
        # else:
            # axis_dict[axis].append(tors_name)
        # if a1_key in axis and a2_key in axis:
            # bad_tors.append(tors_name)
    # for axis in axis_dict:
        # if len(axis_dict[axis]) > 1:
            # coord1 = coo_dct[axis_dict[axis][0]]
            # coord2 = coo_dct[axis_dict[axis][1]]
            # if new_h_idx in coord1[0]:
                # bad_tors.append(axis_dict[axis][0])
            # if new_h_idx in coord2[0]:
                # bad_tors.append(axis_dict[axis][1])
    # for tors_name in pot_tors_names:
        # if tors_name not in bad_tors:
            # axis = coo_dct[tors_name][0][1:3]
            # grp1 = [axis[1]] + (
                # list(automol.graph.branch_atom_keys(gra, axis[0], axis) -
                     # set(axis)))
            # grp2 = [axis[0]] + (
                # list(automol.graph.branch_atom_keys(gra, axis[1], axis) -
                     # set(axis)))
            # if not ((new_h_idx in grp1 and a1_idx in grp2) or
                    # (new_h_idx in grp2 and a1_idx in grp1)):
                 # tors_names.append(tors_name)

    frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
    brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)

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

        # print('rct_zmas', automol.zmatrix.string(rct_zmas[0]))
        # print('rct_gras', rct_gras)
        # print('atm_key_dct', atm_key_dct)
        rcts_gra = automol.graph.relabel(rct_gras[0], atm_key_dct)
        # rcts_gra = automol.graph.union_from_sequence(new_rct_gra)
    else:
        rcts_gra = automol.graph.union_from_sequence(rct_gras)

    ret = ts_zma, dist_name, frm_bnd_key, brk_bnd_key, tors_names, rcts_gra

    print('ts_zma test:', automol.zmatrix.string(ts_zma))
    print('tors_names test:', tors_names)

    import sys
    sys.exit()
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

        frm_bnd_key = shift_vals_from_dummy(frm_bnd_key, ts_zma)
        ret = ts_zma, dist_name, brk_dist_name, frm_bnd_key, tors_names

    return ret


def beta_scission(rct_zmas, prd_zmas):
    """ z-matrix for a beta-scission reaction
    """
    ret = None
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
        ts_zma = automol.zmatrix.standard_form(ts_zma)
        tors_names = automol.zmatrix.torsion_coordinate_names(ts_zma)

        brk_bnd_key = shift_vals_from_dummy(brk_bnd_key, ts_zma)

        # Build the reactants graph
        rcts_gra = automol.graph.union_from_sequence(rct_gras)        

        ret = ts_zma, dist_name, brk_bnd_key, tors_names, rcts_gra

    return ret


def shift_vals_from_dummy(vals, zma):
    """ Shift a set of values using remdummy
        Shift requires indices be 1-indexed
    """
    type_ = type(vals)

    dummy_idxs = automol.zmatrix.atom_indices(zma, sym='X')

    shift_vals = []
    for val in vals:
        shift = 0
        for dummy in dummy_idxs:
            if val >= dummy:
                shift += 1
        shift_vals.append(val+shift)

    shift_vals = type_(shift_vals)

    return shift_vals
