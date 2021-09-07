""" Common utilities for reaction classes
"""

import itertools
import automol.graph
from automol.graph import ts
from automol.par import ReactionClass
from automol.reac._reac import reactant_graphs, product_graphs


def hydrogen_migration_atom_keys(rxn):
    """ Obtain the atoms involved in a hydrogen migration reaction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom, and
    a neighbor to the attacking atom along the chain to the donating atom
    :rtype: (int, int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    tra_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    don_key, = brk_bnd_key - frm_bnd_key

    gra = ts.reactants_graph(rxn.forward_ts_graph)
    path = automol.graph.shortest_path_between_atoms(gra, att_key, don_key)
    ngb_key = automol.graph.atom_neighbor_atom_key(
        gra, att_key, incl_atm_keys=path)
    return att_key, tra_key, don_key, ngb_key


def ring_forming_scission_atom_keys(rxn):
    """ Obtain the atoms involved in a ring-forming scission reaction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom
    :rtype: (int, int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    tra_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    don_key, = brk_bnd_key - frm_bnd_key
    return att_key, tra_key, don_key


def ring_forming_scission_chain(rxn):
    """ Obtain the chain in a ring-forming scission reaction from the donating
    attom to the attacking atom.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: atoms along the chain
    :rtype: tuple[int]
    """
    att_key, _, don_key = ring_forming_scission_atom_keys(rxn)
    gra = ts.reactants_graph(rxn.forward_ts_graph)
    path = automol.graph.shortest_path_between_atoms(gra, don_key, att_key)
    return tuple(path)


def hydrogen_abstraction_atom_keys(rxn):
    """ Obtain the atoms involved in a hydrogen abstraction, sorted in
    canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the donating atom
    :rtype: (int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    hyd_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    don_key, = brk_bnd_key - frm_bnd_key
    return att_key, hyd_key, don_key


def hydrogen_abstraction_is_sigma(rxn):
    """ Is this a sigma radical hydrogen abstraction?

    :param rxn: the reaction object
    :type rxn: Reaction
    :rtype: bool
    """
    assert rxn.class_ == ReactionClass.Typ.HYDROGEN_ABSTRACTION
    tsg = rxn.forward_ts_graph
    rct_gra = automol.graph.ts.reactants_graph(tsg)
    sig_rad_keys = automol.graph.sigma_radical_atom_keys(rct_gra)

    brk_bnd_key, = ts.breaking_bond_keys(tsg)
    frm_bnd_key, = ts.forming_bond_keys(tsg)
    rad_key, = frm_bnd_key - brk_bnd_key
    return rad_key in sig_rad_keys


def elimination_breaking_bond_keys(rxn):
    """ Obtain the breaking bonds for an elimination reaction, sorted in canonical
    order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the breaking bond keys
    :rtype: (frozenset[int], frozenset[int])
    """
    assert rxn.class_ == ReactionClass.Typ.ELIMINATION
    tsg = rxn.forward_ts_graph
    frm_bnd_key, = ts.forming_bond_keys(tsg)
    brk_bnd_keys = ts.breaking_bond_keys(tsg)
    brk_bnd_key1, brk_bnd_key2 = brk_bnd_keys

    symbs = automol.graph.atom_symbols(tsg)

    if len(frm_bnd_key | brk_bnd_key1 | brk_bnd_key2) > 3:
        # for ring_size > 3: use brk-bnd that doesn't involve atoms in frm bond
        scn_brk_bnd_key = None
        for brk_bnd_key in brk_bnd_keys:
            if not frm_bnd_key & brk_bnd_key:
                scn_brk_bnd_key = brk_bnd_key
    else:
        # if one brk bnd doesn't have H, use that, else use arbitrary brk bnd
        scn_brk_bnd_key = None
        for brk_bnd_key in brk_bnd_keys:
            brk_symbs = tuple(symbs[key] for key in brk_bnd_key1)
            if 'H' not in brk_symbs:
                scn_brk_bnd_key = brk_bnd_key
                break
        if scn_brk_bnd_key is None:
            scn_brk_bnd_key = brk_bnd_key1

    brk_bnd_key1 = scn_brk_bnd_key
    brk_bnd_key2, = brk_bnd_keys - {brk_bnd_key1}

    return brk_bnd_key1, brk_bnd_key2


def insertion_forming_bond_keys(rxn):
    """ Obtain the forming bonds for an insertion reaction, sorted in canonical
    order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the forming bond keys
    :rtype: (frozenset[int], frozenset[int])
    """

    assert rxn.class_ == ReactionClass.Typ.INSERTION
    # Choose the forming bond with the fewest neighbors, to get the terminal
    # atom if there is one.
    tsg = rxn.forward_ts_graph
    # For tricky reasons, we need to sort in descending order here.
    # The reason is that the ordering of bonds here determines the ordering of
    # atoms in the z-matrix, which means this function could otherwise yield
    # inconsistent results (i.e. the two bonds will reverse order).
    # If needed, we could add sorting that is based on the symbols of the
    # atoms instead.
    frm_bnd_keys = reversed(sorted(ts.forming_bond_keys(tsg), key=sorted))
    frm_bnd_keys = sorted(
        frm_bnd_keys, key=lambda x: automol.graph.atom_count(
            automol.graph.bond_neighborhood(tsg, x)))

    return tuple(frm_bnd_keys)


def substitution_atom_keys(rxn):
    """ Obtain the atoms involved in a substitution reaction, sorted in
        canonical order.

    :param rxn: the reaction object
    :type rxn: Reaction
    :returns: the attacking atom, the transferring atom, the leaving atom
    :rtype: (int, int, int)
    """
    frm_bnd_key, = ts.forming_bond_keys(rxn.forward_ts_graph)
    brk_bnd_key, = ts.breaking_bond_keys(rxn.forward_ts_graph)
    tra_key, = frm_bnd_key & brk_bnd_key
    att_key, = frm_bnd_key - brk_bnd_key
    lea_key, = brk_bnd_key - frm_bnd_key
    return att_key, tra_key, lea_key


# Conversion stuff
def reaction_inchis(rxn, stereo=True):
    """ Get inchis for one side of reactions
    """
    rct_ichs = tuple(automol.graph.inchi(gra, stereo=stereo)
                     for gra in reactant_graphs(rxn))
    prd_ichs = tuple(automol.graph.inchi(gra, stereo=stereo)
                     for gra in product_graphs(rxn))
    return (rct_ichs, prd_ichs)


# Get a reaction object from various identifiers
def rxn_objs_from_inchi(rct_ichs, prd_ichs, indexing='geo'):
    """ Generate obj
    """

    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    prd_geos = list(map(automol.inchi.geometry, prd_ichs))

    return rxn_objs_from_geometry(
        rct_geos, prd_geos, indexing=indexing)


def rxn_objs_from_smiles(rct_smis, prd_smis, indexing='geo'):
    """ Generate obj
    """

    # Is this adding stero? prob should?
    rct_ichs = list(map(automol.smiles.inchi, rct_smis))
    prd_ichs = list(map(automol.smiles.inchi, prd_smis))

    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    prd_geos = list(map(automol.inchi.geometry, prd_ichs))

    return rxn_objs_from_geometry(
        rct_geos, prd_geos, indexing=indexing)


def rxn_objs_from_zmatrix(rct_zmas, prd_zmas, indexing='geo'):
    """ Generate rxn obj
    """

    rct_geos = list(map(automol.zmat.geometry, rct_zmas))
    prd_geos = list(map(automol.zmat.geometry, prd_zmas))

    return rxn_objs_from_geometry(
        rct_geos, prd_geos, indexing=indexing)


def rxn_objs_from_geometry(rct_geos, prd_geos, indexing='geo'):
    """ from
    """

    # Identify the reaction based on the reactants and products
    rct_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, rct_geos)))
    prd_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, prd_geos)))

    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = automol.reac.find(rct_gras, prd_gras)

    # Obtain the reaction objects and structures to return
    rxn_objs = tuple()
    for rxn in rxns:
        std_rxn, std_rgeos, std_pgeos = (
            automol.reac.standard_keys_with_sorted_geometries(
                rxn, rct_geos, prd_geos))
        ts_geo = automol.reac.ts_geometry(std_rxn, std_rgeos, log=False)

        # Determine which geometries to store
        if indexing == 'geo':
            rxn_objs += ((std_rxn, ts_geo, std_rgeos, std_pgeos),)
        elif indexing == 'zma':
            ts_zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(
                std_rxn, ts_geo)
            std_zrxn = automol.reac.relabel_for_zmatrix(
                std_rxn, zma_keys, dummy_key_dct)
            rct_zmas = tuple(map(automol.geom.zmatrix, rct_geos))
            prd_zmas = tuple(map(automol.geom.zmatrix, prd_geos))

            rxn_objs += ((std_zrxn, ts_zma, rct_zmas, prd_zmas),)

    return rxn_objs


def assert_is_valid_reagent_graph_list(gras):
    """ Assert that a sequence of graphs has the appropriate form for reactants
    or products in a reaction object

    The sequence is appropriate if every graph is explicit and without stereo
    assignments, and none of them have overlapping atom keys.

    :param gras: the graphs
    :type gras: list

    """
    gras_str = '\n---\n'.join(map(automol.graph.string, gras))
    assert _are_all_explicit(gras), (
        "Implicit hydrogens are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_stereo_assignments(gras), (
        "Stereo assignments are not allowed here!\nGraphs:\n{}"
        .format(gras_str))
    assert _have_no_common_atom_keys(gras), (
        "Overlapping atom keys are not allowed here!\nGraphs:\n{}"
        .format(gras_str))


def _are_all_explicit(gras):
    return all(gra == automol.graph.explicit(gra) for gra in gras)


def _have_no_stereo_assignments(gras):
    return all(gra == automol.graph.without_stereo_parities(gra)
               for gra in gras)


def _have_no_common_atom_keys(gras):
    atm_keys = list(itertools.chain(*map(automol.graph.atom_keys, gras)))
    return len(atm_keys) == len(set(atm_keys))


def argsort_reagents(gras):
    """ Get indices to sort reagents by size, from largest to smallest

    :param gras: the reagent (i.e. reactant or product) graphs
    :returns: indices for sorting the reagents
    :rtype: tuple[int]
    """

    def __sort_value(args):
        _, gra = args
        val = (-automol.graph.heavy_atom_count(gra),
               -automol.graph.atom_count(gra),
               -automol.graph.electron_count(gra))
        return val

    idxs = tuple(idx for idx, gra in sorted(enumerate(gras), key=__sort_value))
    return idxs


def sort_reagents(gras):
    """ Sort reagents by size, from largest to smallest

    :param gras: the reagent (i.e. reactant or product) graphs
    :returns: the reagent graphs, in sorted order
    """
    gras = tuple(gras)
    idxs = argsort_reagents(gras)
    gras = tuple(map(gras.__getitem__, idxs))
    return gras


# if __name__ == '__main__':
#     import automol
#
#     RCT_SMIS = ['CCC(CO[O])OO']
#     PRD_SMIS = ['C=C(CC)OO', 'O[O]']
#
#     rxn_objs = automol.reac.rxn_objs_from_smiles(
#         RCT_SMIS, PRD_SMIS, indexing='zma')
#     zrxn, zma, _, _ = rxn_objs[0]
#     scan_inf = automol.reac.build_scan_info(zrxn, zma)
