""" TS z-matrices for specific reaction classes
"""
import automol.geom
from automol.graph import ts


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

    return zma, rxn
