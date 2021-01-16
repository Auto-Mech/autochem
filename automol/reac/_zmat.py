""" TS z-matrices for specific reaction classes
"""
import automol.geom
import automol.graph
from automol.graph import ts


# Unimolecular reactions
def hydrogen_migration_ts_zmatrix(rxn, ts_geo):
    """ z-matrix for a hydrogen migration transition state geometry

    :param rxn: a hydrogen abstraction Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = rxn.copy()
    tsg = rxn.forward_ts_graph

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(tsg)
    geo, dummy_idx_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn.insert_dummy_atoms_(dummy_idx_dct)

    # 4. Generate a z-matrix for the geometry
    # Start the z-matrix from the forming bond ring
    rng_keys, = ts.forming_rings_atom_keys(tsg)
    frm_bnd_key, = ts.forming_bond_keys(tsg)
    brk_bnd_key, = ts.breaking_bond_keys(tsg)
    hyd_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    # First, cycle the migrating h to the front of the ring keys and, if
    # needed, reverse the ring so that the attacking atom is last
    rng_keys = automol.graph.cycle_ring_atom_key_to_front(
        rng_keys, hyd_key, end_key=att_key)
    # Now, cycle the third-to-last key to the front so that the ring order is:
    #       (atm1, atm2, attacking atom, migrating h atom, ....)
    rng_keys = automol.graph.cycle_ring_atom_key_to_front(
        rng_keys, rng_keys[-3])
    vma, row_keys = automol.graph.vmat.vmatrix(tsg, rng_keys=rng_keys)

    zma_geo = automol.geom.from_subset(geo, row_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, row_keys, dummy_idx_dct


# Bimolecular reactions
def hydrogen_abstraction_ts_zmatrix(rxn, ts_geo):
    """ z-matrix for a hydrogen abstraction transition state geometry

    :param rxn: a hydrogen abstraction Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = rxn.copy()

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))
    # Add a dummy atom over the transferring hydrogen
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    hyd_idx, = frm_bnd_key & brk_bnd_key
    lin_idxs.append(hyd_idx)

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.forward_ts_graph)
    geo, dummy_idx_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn.insert_dummy_atoms_(dummy_idx_dct)

    # 3. Generate a z-matrix for the geometry
    vma, row_keys = automol.graph.vmat.vmatrix(rxn.forward_ts_graph)

    zma_geo = automol.geom.from_subset(geo, row_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, row_keys, dummy_idx_dct
