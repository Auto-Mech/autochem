"""
  Generate special constraint coordinates for the ZMAs for
  various kinds of transision states
"""

import automol


# UNIMOL
def hydrogen_migration(ts_zma, rct_zma, new_h_idx,
                       frm_bnd_keys, brk_bnd_keys):
    """ z-matrix for a hydrogen migration reaction

        tras in the ZMATRIX coordinate system
    """

    # Standardize the zma (maybe move out of here)
    ts_name_dct = automol.zmatrix.standard_names(ts_zma)

    # Set the initial torsional names using the reactant zma
    tors_names = automol.zmatrix.torsion_coordinate_names(rct_zma)
    pot_tors_names = []
    for tors_name in tors_names:
        pot_tors_names.append(ts_name_dct[tors_name])

    # Remove the torsional coordinates that are in the ring
    chn_end = new_h_idx
    for key in frm_bnd_keys:
        if key != new_h_idx:
            chn_end = key
    ring_atoms1 = automol.zmatrix.chain_between(ts_zma, chn_end, new_h_idx)

    chn_end = new_h_idx
    for key in brk_bnd_keys:
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
    tors_keys = [automol.zmatrix.coord_idxs(ts_zma, name) for name in tors_names]

    # Build bond names
    if len(ring_atoms1) > 3:
        bnd_keys = frozenset({ring_atoms1[1], ring_atoms1[2]})
    else:
        bnd_keys = frozenset({})

    return (bnd_keys, frozenset({}), tors_keys)


def elimination(ts_zma, rct_zma, mig_key, a1_idx):
    """  torsionas for an elimination reaction
    """

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

    return tors_names


def ring_forming_scission():
    """
    """

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

    angles = [0., 2.*numpy.pi/3, 4.*numpy.pi/3]
    # angles = [0., numpy.pi/3., 2.*numpy.pi/3, 3.*numpy.pi/3., 4.*numpy.pi/3, 5*numpy.pi/3.]
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
        #     if tors_name not in const_tors_names:
        #         const_tors_names.append(tors_name)

    ts_zma = ts_zma_max

    # (v) vary angles to decrease rad_atm to new_rad_atm to < 2.25 Ang
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

    return tors_names, const_tors_names


# BIMOL
def bimol(rct1_zma, rct2_zma,
          ts_name_dct, frm_bnd_keys, rct2_gra):
    """ tors names for hydrogen abstraction
    """

    # Get the torsional coordinates of the transition state
    rct1_tors_names = automol.zmatrix.torsion_coordinate_names(rct1_zma)
    rct2_tors_names = automol.zmatrix.torsion_coordinate_names(rct2_zma)
    tors_names = (
        tuple(map(ts_name_dct.__getitem__, rct1_tors_names)) +
        tuple(map(ts_name_dct.__getitem__, rct2_tors_names))
    )

    if 'babs2' in ts_name_dct:
        geo1 = automol.convert.zmatrix.geometry(rct1_zma)
        if not automol.geom.is_linear(geo1):
            tors_name = ts_name_dct['babs2']
            tors_names += (tors_name,)

    if 'babs3' in ts_name_dct and _include_babs3(frm_bnd_keys, rct2_gra):
        tors_name = ts_name_dct['babs3']
        tors_names += (tors_name,)


def _include_babs3(frm_bnd_keys, rct2_gra):
    """Should we include babs3?
    """
    include = False
    atm_ngbs = automol.graph.atom_neighbor_keys(rct2_gra)
    is_terminal = False
    for atm in list(frm_bnd_keys):
        if atm in atm_ngbs:
            if len(atm_ngbs[atm]) == 1:
                is_terminal = True
    if len(atm_ngbs.keys()) > 2 and is_terminal:
        include = True
    return include



