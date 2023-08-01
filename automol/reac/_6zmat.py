""" TS z-matrices for specific reaction classes
"""
import copy
import automol.geom
import automol.graph
from automol.par import ReactionClass
from automol.graph import ts
from automol.zmat import distance_coordinate_name
from automol.reac._0core import Reaction
from automol.reac._0core import add_dummy_atoms
from automol.reac._1util import hydrogen_migration_atom_keys
from automol.reac._1util import ring_forming_scission_chain
from automol.reac._1util import elimination_breaking_bond_keys
from automol.reac._1util import insertion_forming_bond_keys
from automol.reac._1util import hydrogen_abstraction_atom_keys
from automol.reac._1util import hydrogen_abstraction_is_sigma
from automol.reac._1util import substitution_atom_keys


# Unimolecular reactions
# 1. Hydrogen migrations
def hydrogen_migration_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for a hydrogen migration transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    # Start the z-matrix from the forming bond ring
    rng_keys, = ts.forming_rings_atom_keys(rxn.ts_graph)
    _, hyd_key, don_key, _ = hydrogen_migration_atom_keys(rxn)
    # Cycle the migrating h to the front of the ring keys and, if
    # needed, reverse the ring so that the donating atom is last:
    #       (migrating h atom, attacking atom, ... , donating atom)
    # This ensures that the forming bond coordinate is included in the z-matrix
    # and drops the breaking bond coordinate from the z-matrix.
    rng_keys = automol.graph.cycle_ring_atom_key_to_front(
        rng_keys, hyd_key, end_key=don_key)
    vma, zma_keys = automol.graph.vmat.vmatrix(
        rxn.ts_graph, rng_keys=rng_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 2. Beta scissions
def beta_scission_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for a beta scission transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    vma, zma_keys = automol.graph.vmat.vmatrix(rxn.ts_graph)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 3. Ring-forming scissions
def ring_forming_scission_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for a ring-forming scission transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    chain_keys = ring_forming_scission_chain(rxn)
    rng_keys = list(reversed(chain_keys[1:]))
    vma, zma_keys = automol.graph.vmat.vmatrix(rxn.ts_graph,
                                               rng_keys=rng_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 4. Eliminations
def elimination_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for an elimination transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    rng_keys, = ts.forming_rings_atom_keys(rxn.ts_graph)
    frm_bnd_key, = ts.ts_forming_bond_keys(rxn.ts_graph)
    _, brk_bnd_key = elimination_breaking_bond_keys(rxn)
    # Cycle the ring keys such that the atom closest to the forming bond is the
    # beginning of the ring and the other atom is the end
    if brk_bnd_key & frm_bnd_key:
        key1, = brk_bnd_key & frm_bnd_key
        key2, = brk_bnd_key - frm_bnd_key
    else:
        path = automol.graph.shortest_path_between_groups(
            rxn.ts_graph, brk_bnd_key, frm_bnd_key)
        key1, = brk_bnd_key & set(path)
        key2, = brk_bnd_key - set(path)
    rng_keys = automol.graph.cycle_ring_atom_key_to_front(
        rng_keys, key1, end_key=key2)

    vma, zma_keys = automol.graph.vmat.vmatrix(
        rxn.ts_graph, rng_keys=rng_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# Bimolecular reactions
# 1. Hydrogen abstractions
def hydrogen_abstraction_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for a hydrogen abstraction transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))
    # Add a dummy atom over the transferring hydrogen
    att_key, hyd_key, _ = hydrogen_abstraction_atom_keys(rxn)
    lin_idxs.append(hyd_key)

    if hydrogen_abstraction_is_sigma(rxn):
        if att_key not in lin_idxs:
            lin_idxs.append(att_key)

    lin_idxs = sorted(lin_idxs)

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    tsg = rxn.ts_graph
    rct1_keys, rct2_keys = rxn.reactants_keys
    vma, zma_keys = automol.graph.vmat.vmatrix(tsg, rct1_keys)
    vma, zma_keys = automol.graph.vmat.continue_vmatrix(
        tsg, rct2_keys, vma, zma_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 2. Additions
def addition_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for an addition transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = rxn.ts_graph

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    tsg = rxn.ts_graph
    rct1_keys, rct2_keys = rxn.reactants_keys
    vma, zma_keys = automol.graph.vmat.vmatrix(tsg, rct1_keys)
    vma, zma_keys = automol.graph.vmat.continue_vmatrix(
        tsg, rct2_keys, vma, zma_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 3. Insertions
def insertion_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for an insertion transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = rxn.ts_graph

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    rng_keys, = ts.forming_rings_atom_keys(rxn.ts_graph)
    brk_bnd_key, = ts.ts_breaking_bond_keys(rxn.ts_graph)
    # Drop one of the forming bonds from the z-matrix by sorting the ring atom
    # keys to exclude it.
    _, frm_bnd_key = insertion_forming_bond_keys(rxn)
    # Cycle the ring keys such that the atom closest to the breaking bond is
    # the beginning of the ring and the other atom is the end
    if frm_bnd_key & brk_bnd_key:
        key1, = frm_bnd_key & brk_bnd_key
        key2, = frm_bnd_key - brk_bnd_key
    else:
        path = automol.graph.shortest_path_between_groups(
            rxn.ts_graph, frm_bnd_key, brk_bnd_key)
        key2, = frm_bnd_key - set(path)
    rng_keys = automol.graph.cycle_ring_atom_key_to_front(
        rng_keys, key1, end_key=key2)

    vma, zma_keys = automol.graph.vmat.vmatrix(
        rxn.ts_graph, rng_keys=rng_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


# 4. Substitutions
def substitution_ts_zmatrix(rxn: Reaction, ts_geo):
    """ z-matrix for a substitution transition state geometry

    :param rxn: a Reaction object
    :param ts_geo: a transition state geometry
    """
    rxn = copy.deepcopy(rxn)
    rxn.ts_graph = rxn.ts_graph

    # 1. Get keys to linear or near-linear atoms
    lin_idxs = list(automol.geom.linear_atoms(ts_geo))
    # Add a dummy atom over the transferring hydrogen
    _, tra_key, _ = substitution_atom_keys(rxn)
    lin_idxs.append(tra_key)

    # 2. Add dummy atoms over the linear atoms
    rcts_gra = ts.reactants_graph(rxn.ts_graph)
    geo, dummy_key_dct = automol.geom.insert_dummies_on_linear_atoms(
        ts_geo, lin_idxs=lin_idxs, gra=rcts_gra)

    # 3. Add dummy atoms to the Reaction object as well
    rxn = add_dummy_atoms(rxn, dummy_key_dct)

    # 4. Generate a z-matrix for the geometry
    tsg = rxn.ts_graph
    rct1_keys, rct2_keys = rxn.reactants_keys
    vma, zma_keys = automol.graph.vmat.vmatrix(tsg, rct1_keys)
    vma, zma_keys = automol.graph.vmat.continue_vmatrix(
        tsg, rct2_keys, vma, zma_keys)

    zma_geo = automol.geom.from_subset(geo, zma_keys)
    zma = automol.zmat.from_geometry(vma, zma_geo)

    return zma, zma_keys, dummy_key_dct


def ts_zmatrix(rxn: Reaction, ts_geo):
    """ reaction-class-specific embedding info

    :param rxn: a hydrogen migration Reaction object
    :param ts_geo: the TS geometry
    :returns: the TS z-matrix, the row keys, and the dummy index dictionary
    """
    function_dct = {
        # unimolecular
        ReactionClass.Typ.HYDROGEN_MIGRATION: hydrogen_migration_ts_zmatrix,
        ReactionClass.Typ.BETA_SCISSION: beta_scission_ts_zmatrix,
        ReactionClass.Typ.RING_FORM_SCISSION: ring_forming_scission_ts_zmatrix,
        ReactionClass.Typ.ELIMINATION: elimination_ts_zmatrix,
        # bimolecular
        ReactionClass.Typ.HYDROGEN_ABSTRACTION:
        hydrogen_abstraction_ts_zmatrix,
        ReactionClass.Typ.ADDITION: addition_ts_zmatrix,
        ReactionClass.Typ.INSERTION: insertion_ts_zmatrix,
        ReactionClass.Typ.SUBSTITUTION: substitution_ts_zmatrix,
    }

    fun_ = function_dct[rxn.class_]
    ret = fun_(rxn, ts_geo)
    return ret


# Z-Matrix coordinate functions
def zmatrix_coordinate_names(zrxn: Reaction, zma):
    """ Get the Z-matrix coordinate names for the forming and
        breaking bonds of a reaction

        It is not always guaranteed that the bond keys will be
        present in the transition state Z-Matrix. For these cases,
        the name will be returned as None.

        :param zrxn: a Reaction object
        :rtype: str
    """

    def _zma_names(zma, bnd_keys):
        _names = ()
        for keys in bnd_keys:
            try:
                name = distance_coordinate_name(zma, *keys)
            except AssertionError:
                name = None
            _names += (name,)

        return _names

    frm_bnd_keys = ts.ts_forming_bond_keys(zrxn.ts_graph)
    brk_bnd_keys = ts.ts_breaking_bond_keys(zrxn.ts_graph)

    return (_zma_names(zma, frm_bnd_keys), _zma_names(zma, brk_bnd_keys))
