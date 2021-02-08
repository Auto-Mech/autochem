""" random functions that I am gonna distribute to other parts of automol
"""

def rxn_angle(grxn, zma):
    """ Calculate the angle over a forming-and-breaking bond
    """

    geo = automol.zmat.geometry(zma)

    frm_bnd_keys = automol.reac.forming_bond_keys(grxn)
    brk_bnd_keys = automol.reac.breaking_bond_keys(grxn)

    angle = None
    if 'abstraction' in grxn.class_ or 'addition' in grxn.class_:
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
