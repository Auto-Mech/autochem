"""
  Util functions for ts.py
"""

import functools
import automol


def reorder_zmatrix_for_migration(zma, a_idx, h_idx):
    """ performs z-matrix reordering operations required to
        build proper z-matrices for hydrogen migrations
    """

    # initialize zmat components needed later
    symbols = automol.zmatrix.symbols(zma)

    # Get the longest chain for all the atoms
    _, gras = shifted_standard_zmas_graphs([zma], remove_stereo=True)
    gra = functools.reduce(automol.graph.union, gras)
    xgr1, = automol.graph.connected_components(gra)
    chains_dct = automol.graph.atom_longest_chains(xgr1)

    # find the longest heavy-atom chain for the forming atom
    form_chain = chains_dct[a_idx]

    # get the indices used for reordering
    # get the longest chain from the bond forming atom (not including H)
    order_idxs = [idx for idx in form_chain
                  if symbols[idx] != 'H' and idx != h_idx]

    # add all the heavy-atoms not in the chain
    for i, atom in enumerate(symbols):
        if i not in order_idxs and atom != 'H':
            order_idxs.append(i)

    # add all the hydrogens
    for i, atom in enumerate(symbols):
        if i != h_idx and atom == 'H':
            order_idxs.append(i)

    # add the migrating atoms
    order_idxs.append(h_idx)

    # get the geometry and redorder it according to the order_idxs list
    geo = [list(x) for x in automol.zmatrix.geometry(zma)]
    geo2 = tuple(tuple(geo[idx]) for idx in order_idxs)

    # convert to a z-matrix
    zma_ret = automol.geom.zmatrix(geo2)

    return zma_ret


def include_babs3(frm_bnd, rct2_gra):
    """Should we include babs3?
    """
    include = False
    atm_ngbs = automol.graph.atom_neighbor_keys(rct2_gra)
    is_terminal = False
    for atm in list(frm_bnd):
        if atm in atm_ngbs:
            if len(atm_ngbs[atm]) == 1:
                is_terminal = True
    if len(atm_ngbs.keys()) > 2 and is_terminal:
        include = True
    return include


def shifted_standard_zmas_graphs(zmas, remove_stereo=False):
    """ Generate zmas and graphs from input zmas
        shfited and in their standard form
    """
    conv = functools.partial(
        automol.convert.zmatrix.graph, remove_stereo=remove_stereo)
    gras = list(map(conv, zmas))
    shift = 0
    for idx, (zma, gra) in enumerate(zip(zmas, gras)):
        zmas[idx] = automol.zmatrix.standard_form(zma, shift=shift)
        gras[idx] = automol.graph.transform_keys(gra, lambda x: x+shift)
        shift += len(automol.graph.atoms(gra))
    zmas = tuple(zmas)
    gras = tuple(map(automol.graph.without_dummy_atoms, gras))
    return zmas, gras


def join_atom_keys(zma, atm1_key):
    """ returns available join atom keys (if available) and a boolean
    indicating whether the atoms are in a chain or not
    """
    gra = automol.convert.zmatrix.graph(zma)
    atm1_chain = (
        automol.graph.atom_longest_chains(gra)[atm1_key])
    atm1_ngb_keys = (
        automol.graph.atom_neighbor_keys(gra)[atm1_key])
    if len(atm1_chain) == 1:
        atm2_key = None
        atm3_key = None
        chain = False
    elif len(atm1_chain) == 2 and len(atm1_ngb_keys) == 1:
        atm2_key = atm1_chain[1]
        atm3_key = None
        chain = False
    elif len(atm1_chain) == 2:
        atm2_key = atm1_chain[1]
        atm3_key = sorted(atm1_ngb_keys - {atm2_key})[0]
        chain = False
    else:
        atm2_key, atm3_key = atm1_chain[1:3]
        chain = True

    return atm2_key, atm3_key, chain


def reorder_zma_for_radicals(zma, rad_idx):
    """ Creates a zmatrix where the radical atom is the first entry
        in the zmatrix
    """
    geo = automol.zmatrix.geometry(zma)
    geo_swp = automol.geom.swap_coordinates(geo, 0, rad_idx)
    zma_swp = automol.geom.zmatrix(geo_swp)

    return zma_swp
# old elim function; delete when we want
# def concerted_unimolecular_elimination2(rct_zmas, prd_zmas):
#     """ z-matrix for a concerted unimolecular elimination reaction
#     """
#     ret = None, None, None, None
#     prd_zmas, prd_gras = shifted_standard_forms_with_gaphs(
#         prd_zmas, remove_stereo=True)
#     if len(rct_zmas) == 1:
#         count = 1
#         while True:
#             rct_zmas, rct_gras = shifted_standard_forms_with_gaphs(
#                 rct_zmas, remove_stereo=True)
#             init_zma, = rct_zmas
#
#             tras = automol.graph.reac.elimination(rct_gras, prd_gras)
#             if tras is not None:
#                 min_dist = 100.
#                 frm_bnd_key = None
#                 for tra_i in tras:
#                     # Get the bond formation and breaking keys
#                     bnd_key, = automol.graph.trans.formed_bond_keys(tra_i)
#                     geo = automol.zmatrix.geometry(rct_zmas[0])
#                     dist = automol.geom.distance(geo, *list(bnd_key))
#                     if dist < min_dist:
#                         min_dist = dist
#                         frm_bnd_key = bnd_key
#                         tra = tra_i
#                 brk_keys = automol.graph.trans.broken_bond_keys(tra)
#                 brk_bnd_key1, brk_bnd_key2 = brk_keys
#                 init_zma, = rct_zmas
#
#                 # Get index for migrating atom (or bond-form atom in group)
#                 for bnd_key in (brk_bnd_key1, brk_bnd_key2):
#                     if bnd_key & frm_bnd_key:
#                         key_iter = iter(bnd_key & frm_bnd_key)
#                         mig_key = next(key_iter)
#                         a1_idx = next(key_iter)
#
#                 # Get chain for redefining the rc1_atm1_key z-matrix entries
#                 _, gras = shifted_standard_forms_with_gaphs(
#                     [init_zma], remove_stereo=True)
#                 gra = functools.reduce(automol.graph.union, gras)
#                 xgr1, = automol.graph.connected_components(gra)
#                 atm1_neighbors = _atom_neighbor_keys(xgr1)[a1_idx]
#                 for idx in atm1_neighbors:
#                     num_keys = len(_atom_neighbor_keys(xgr1)[idx])
#                     if idx != mig_key and num_keys > 1:
#                         a2_idx = idx
#                 atm2_neighbors = _atom_neighbor_keys(xgr1)[a2_idx]
#                 for idx in atm2_neighbors:
#                     if idx != mig_key:
#                         a3_idx = idx
#
#                 mig_redef_keys = (a1_idx, a2_idx, a3_idx)
#
#                 # determine if the zmatrix needs to be rebuilt by x2z
#              # determines if the hydrogen atom is used to define other atoms
#                 rebuild = False
#                 if any(idx > mig_key for idx in mig_redef_keys):
#                     rebuild = True
#
#                 # rebuild zmat and go through while loop again if needed
#              # shift order of cartesian coords & rerun x2z to get a new zmat
#                 # else go to next stage
#                 if rebuild:
#                     reord_zma = _reorder_zmatrix_for_migration(
#                         init_zma, a1_idx, mig_key)
#                     rct_zmas = [reord_zma]
#                     count += 1
#                     if count == 3:
#                         break
#                 else:
#                     rct_zma = init_zma
#                     break
#             else:
#                 return None, None, None, None
#
#         # Get the name of the coordinates for the bonds breaking
#         brk_names = []
#         for brk_key in (brk_bnd_key1, brk_bnd_key2):
#             brk_names.append(
#                 automol.zmatrix.bond_key_from_idxs(rct_zma, brk_key))
#
#         # increase the length of the bond breaking coordinates
#         for name in brk_names:
#             val = automol.zmatrix.values(ts_zma)[name]
#             ts_zma = automol.zmatrix.set_values(
#                 ts_zma, {name: val + (0.4 * ANG2BOHR)})
#
#         # determine the new coordinates
#         rct_geo = automol.zmatrix.geometry(rct_zma)
#         distance = automol.geom.distance(
#             rct_geo, mig_key, a1_idx)
#         angle = automol.geom.central_angle(
#             rct_geo, mig_key, a1_idx, a2_idx)
#         dihedral = automol.geom.dihedral_angle(
#             rct_geo, mig_key, a1_idx, a2_idx, a3_idx)
#         # Reset the keys for the migrating H atom
#         new_idxs = (a1_idx, a2_idx, a3_idx)
#         key_dct = {mig_key: new_idxs}
#         ts_zma = automol.zmatrix.set_keys(rct_zma, key_dct)
#
#         # Reset the values in the value dict
#         mig_names = automol.zmatrix.name_matrix(ts_zma)[mig_key]
#         ts_zma = automol.zmatrix.set_values(
#             ts_zma, {mig_names[0]: distance,
#                      mig_names[1]: angle,
#                      mig_names[2]: dihedral}
#         )
#
#         # standardize the ts zmat and get tors and dist coords
#         coo_dct = automol.zmatrix.coordinates(ts_zma)
#         dist_coo_key = tuple(reversed(sorted(frm_bnd_key)))
#         dist_name = next(coo_name for coo_name, coo_keys in coo_dct.items()
#                          if dist_coo_key in coo_keys)
#         ts_name_dct = automol.zmatrix.standard_names(ts_zma)
#         dist_name = ts_name_dct[dist_name]
#         ts_zma = automol.zmatrix.standard_form(ts_zma)
#
#         # Get the name of the coordinate of the other bond that is breaking
#         for brk_key in (brk_bnd_key1, brk_bnd_key2):
#             if not brk_key.intersection(frm_bnd_key):
#                 brk_dist_name = automol.zmatrix.bond_key_from_idxs(
#                     ts_zma, brk_key)
#
#         # get full set of potential torsional coordinates
#         pot_tors_names = automol.zmatrix.torsion_coordinate_names(rct_zma)
#
#       # remove the torsional coordinates that would break reaction coordinate
#         gra = automol.zmatrix.graph(ts_zma, remove_stereo=True)
#         coo_dct = automol.zmatrix.coordinates(ts_zma)
#         tors_names = []
#         for tors_name in pot_tors_names:
#             axis = coo_dct[tors_name][0][1:3]
#             grp1 = [axis[1]] + (
#                 list(automol.graph.branch_atom_keys(gra, axis[0], axis) -
#                      set(axis)))
#             grp2 = [axis[0]] + (
#                 list(automol.graph.branch_atom_keys(gra, axis[1], axis) -
#                      set(axis)))
#             if not ((mig_key in grp1 and a1_idx in grp2) or
#                     (mig_key in grp2 and a1_idx in grp1)):
#                 tors_names.append(tors_name)
#
#         ret = ts_zma, dist_name, brk_dist_name, frm_bnd_key, tors_names
#
#         return ret
