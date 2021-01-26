""" New functions that are unsorted
"""

def z_atom_closest_to(zma, idx1, idx2, chk_idxs):
    """
    """
    geo = automol.zmat.geometry(zma) 
    return def atom_closest_to(geo, idx1, idx2, chk_idxs)


def atom_closest_to(geo, idx1, idx2, chk_idxs):
    """ Check to see whether the radical atom is still closest to the bond
        formation site.
    """

    atom_closest = True

    dist1 = automol.geom.distance(geo, idx1, idx2)
    for idx in chk_idxs:
        if dist1 < automol.geom.distance(geo, idx1, idx) - 0.01:
                atom_closest = False

    return atom_closest


def calc_rxn_angle(ts_zma, frm_bnd_keys, brk_bnd_keys, rxn_class):
    """ Calculate the angle over a forming-and-breaking bond
    """

    angle = None
    if 'abstraction' in rxn_class or 'addition' in rxn_class:
        if frm_bnd_keys and brk_bnd_keys:

            ang_atms = [0, 0, 0]
            cent_atm = list(set(brk_bnd_keys) & set(frm_bnd_keys))
            if cent_atm:
                ang_atms[1] = cent_atm[0]
                for idx in brk_bnd_keys:
                    if idx != ang_atms[1]:
                        ang_atms[0] = idx
                for idx in frm_bnd_keys:
                    if idx != ang_atms[1]:
                        ang_atms[2] = idx

                geom = automol.zmatrix.geometry(ts_zma)
                angle = automol.geom.central_angle(
                    geom, *ang_atms)

    return angle

