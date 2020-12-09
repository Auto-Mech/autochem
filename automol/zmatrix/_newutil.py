""" new stuff
"""

import automol


def remove_dummies(zma, frm_key, brk_key, geo=None):
    """get zma and bond key idxs without dummy atoms
    """
    zgeo = automol.zmatrix.geometry(zma)
    brk_key2 = None
    if isinstance(brk_key, list):
        brk_key, brk_key2 = brk_key
    dummy_idxs = automol.geom.dummy_atom_indices(zgeo)
    for idx in dummy_idxs:
        if frm_key:
            frm1, frm2 = frm_key
            if idx < frm1:
                frm1 -= 1
            if idx < frm2:
                frm2 -= 1
            frm_key = frozenset({frm1, frm2})
        if brk_key:
            brk1, brk2 = brk_key
            if idx < brk1:
                brk1 -= 1
            if idx < brk2:
                brk2 -= 1
            brk_key = frozenset({brk1, brk2})
        if brk_key2:
            brk3, brk4 = brk_key2
            if idx < brk3:
                brk3 -= 1
            if idx < brk4:
                brk4 -= 1
            brk_key2 = frozenset({brk3, brk4})
    if not geo:
        geo = automol.geom.without_dummy_atoms(zgeo)
    
    gra = automol.geom.graph(geo)
    #if dummy_idxs:        
    #    if not geo:
    #        geo = automol.geom.without_dummy_atoms(zgeo)
    #    zma = automol.geom.zmatrix(geo, [frm_key, brk_key])
    #    atm_ord = automol.geom.zmatrix_atom_ordering(geo, [frm_key, brk_key])
    #    frm_key = frozenset({atm_ord[atm] for atm in frm_key})    
    #    brk_key = frozenset({atm_ord[atm] for atm in brk_key})    
    #return zma, frm_key, brk_key
    return gra, frm_key, brk_key, brk_key2

