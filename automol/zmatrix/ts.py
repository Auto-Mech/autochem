"""
Import all the ts builder functions
"""
import automol.inchi
import automol.graph
import automol.geom
import automol.zmatrix
from automol.zmatrix._unimol_ts import min_hyd_mig_dist
from automol.zmatrix._unimol_ts import hydrogen_migration
from automol.zmatrix._unimol_ts import min_unimolecular_elimination_dist
from automol.zmatrix._unimol_ts import concerted_unimolecular_elimination
from automol.zmatrix._unimol_ts import beta_scission
from automol.zmatrix._bimol_ts import insertion
from automol.zmatrix._bimol_ts import substitution
from automol.zmatrix._bimol_ts import addition
from automol.zmatrix._bimol_ts import hydrogen_abstraction


def zmatrix_reaction_info(ts_zma, rct_gras, prd_gras):
    """ determines reactant graph and formed/broken bonds for the TS zmatrix

    The keys used in the reactant and product graphs here are irrelevant -- the
    return value will use the keys of the TS zmatrix.

    The transformation returned will be one that aligns with the TS z-matrix
    keys, and it will be the transformation that most closely matches the TS
    z-matrix geometry, in terms of the bond distances for the bonds broken and
    formed.
    """
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    ts_gra = automol.zmatrix.connectivity_graph(ts_zma)
    ts_gra = automol.graph.without_dummy_atoms(ts_gra)

    # First, try aligning the reactant graph to the TS graph as-is, using a
    # subgraph isomoprhism
    aligned_rct_gras = _align_reactant_graphs_to_ts(ts_gra, rct_gras)

    # If that doesn't work, bond the closest unbonded atoms and try again
    if aligned_rct_gras is None:
        ts_geo = automol.zmatrix.geometry(ts_zma)

        bnd_key, _ = automol.geom.closest_unbonded_atoms(ts_geo, ts_gra)

        ts_gra = automol.graph.add_bonds(ts_gra, [bnd_key])

        aligned_rct_gras = _align_reactant_graphs_to_ts(ts_gra, rct_gras)

    # Hopefully at this point we have the correctly aligned  reactant graphs.
    # Now, we can determine the transformation (formed/broken keys). Since
    # there may be multiple transformations between these reactants and
    # products, choose the one that most closely matches the TS geometry -- the
    # one in which the formed and broken keys are as short as possible.
    if aligned_rct_gras is not None:
        tras, _, _, _ = automol.graph.reac.classify(aligned_rct_gras, prd_gras)
        tra = _find_closest_transformation(ts_zma, tras)
        rct_gra = automol.graph.union_from_sequence(aligned_rct_gras)
    else:
        tra = None
        rct_gra = None

    return tra, rct_gra


def _align_reactant_graphs_to_ts(ts_gra, rct_gras):
    ts_atm_keys = automol.graph.atom_keys(ts_gra)

    # Get the largest reactant
    rct1_gra = sorted(rct_gras, key=automol.graph.atom_count)[0]

    # See if the TS graph contains a subgraph matching the largest reacant
    rct1_iso_dct = automol.graph.full_subgraph_isomorphism(
        ts_gra, rct1_gra)

    rct_gras_ = None
    if rct1_iso_dct:
        rct1_atm_keys = set(rct1_iso_dct.keys())

        rct1_gra = automol.graph.subgraph(ts_gra, rct1_atm_keys)

        if len(rct_gras) == 1:
            rct_gras_ = [rct1_gra]
        elif len(rct_gras) == 2:
            rct2_atm_keys = ts_atm_keys - rct1_atm_keys
            rct2_gra = automol.graph.subgraph(ts_gra, rct2_atm_keys)
            rct_gras_ = [rct1_gra, rct2_gra]
        else:
            raise NotImplementedError("Doesn't handle termolecular rxns")

    return rct_gras_


def zmatrix_transformation(ts_zma, rct_gras, prd_gras,
                           reactant_connectivity=False):
    """ determines bonds broken and formed, using the TS zmatrix keys

    The keys used in the reactant and product graphs here are irrelevant -- the
    return value will use the keys of the TS zmatrix.

    The transformation returned will be one that aligns with the TS z-matrix
    keys, and it will be the transformation that most closely matches the TS
    z-matrix geometry, in terms of the bond distances for the bonds broken and
    formed.
    """
    ts_gras = None

    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    ts_gra = automol.zmatrix.connectivity_graph(ts_zma)
    ts_gra = automol.graph.without_dummy_atoms(ts_gra)
    if reactant_connectivity:
        # If the connectivity matches the reactants, we can just split up the
        # graph into components
        ts_gras = automol.graph.connected_components(ts_gra)
        assert len(ts_gras) == len(rct_gras)
        assert all(any(automol.graph.backbone_isomorphic(ts_gra, rct_gra)
                       for ts_gra in ts_gras)
                   for rct_gra in rct_gras)
    else:
        # If the connectivity doesn't match the reactants, we have to do some
        # work.
        ts_atm_keys = automol.graph.atom_keys(ts_gra)

        # Get the largest reactant
        rct1_gra = sorted(rct_gras, key=automol.graph.atom_count)[0]

        # See if the TS graph contains a subgraph matching the largest reacant
        rct1_iso_dct = automol.graph.full_subgraph_isomorphism(
            ts_gra, rct1_gra)

        if rct1_iso_dct:
            rct1_atm_keys = set(rct1_iso_dct.keys())

            rct1_gra = automol.graph.subgraph(ts_gra, rct1_atm_keys)

            if len(rct_gras) == 1:
                ts_gras = [rct1_gra]
            elif len(rct_gras) == 2:
                rct2_atm_keys = ts_atm_keys - rct1_atm_keys
                rct2_gra = automol.graph.subgraph(ts_gra, rct2_atm_keys)
                ts_gras = [rct1_gra, rct2_gra]
            else:
                raise NotImplementedError("Doesn't handle termolecular rxns")

    if ts_gras is not None:
        tras, _, _, _ = automol.graph.reac.classify(ts_gras, prd_gras)
        tra = _find_closest_transformation(ts_zma, tras)
    else:
        tra = None

    return tra


def _find_closest_transformation(ts_zma, tras):
    """ find the transformation which best matches the structure of the TS
    zmatrix (the one in which the bonds broken and formed are shortest)
    """
    min_tra = None
    min_dist_val = 1000.

    ts_geo = automol.zmatrix.geometry(ts_zma)
    for tra in tras:
        frm_bnd_keys = automol.graph.trans.formed_bond_keys(tra)
        brk_bnd_keys = automol.graph.trans.broken_bond_keys(tra)

        frm_dists = [automol.geom.distance(ts_geo, *frm_bnd_key)
                     for frm_bnd_key in frm_bnd_keys]
        brk_dists = [automol.geom.distance(ts_geo, *brk_bnd_key)
                     for brk_bnd_key in brk_bnd_keys]

        dist_val = sum(frm_dists) + sum(brk_dists)

        if dist_val < min_dist_val:
            min_dist_val = dist_val
            min_tra = tra

    return min_tra


def zmatrix_reactant_graph(ts_zma, frm_keys, brk_keys):
    """ determine the reactant graph for a TS zmatrix

    (graph keys are aligned with z-matrix rows)
    """
    rct_gra = automol.convert.zmatrix.connectivity_graph(ts_zma)
    rct_gra = automol.graph.without_dummy_atoms(rct_gra)
    rct_gra = automol.graph.add_bonds(rct_gra, brk_keys, check=False)
    rct_gra = automol.graph.remove_bonds(rct_gra, frm_keys, check=False)
    return rct_gra


def zmatrix_product_graph(ts_zma, frm_keys, brk_keys):
    """ determine the product graph for a TS zmatrix

    (graph keys are aligned with z-matrix rows)
    """
    tra = automol.graph.trans.from_data(frm_keys, brk_keys)

    rct_gra = zmatrix_reactant_graph(ts_zma, frm_keys, brk_keys)
    prd_gra = automol.graph.trans.apply(tra, rct_gra)
    return prd_gra


def zmatrix_reactants_inchi(ts_zma, frm_keys, brk_keys, remove_stereo=True):
    """ determine the InChI for the reactants of a TS zmatrix
    """
    rct_gra = zmatrix_reactant_graph(ts_zma, frm_keys, brk_keys)
    rct_ich = automol.convert.graph.inchi(rct_gra, remove_stereo=remove_stereo)
    return rct_ich


def zmatrix_products_inchi(ts_zma, frm_keys, brk_keys, remove_stereo=True):
    """ determine the InChI for the products of a TS zmatrix
    """
    prd_gra = zmatrix_product_graph(ts_zma, frm_keys, brk_keys)
    prd_ich = automol.convert.graph.inchi(prd_gra, remove_stereo=remove_stereo)
    return prd_ich


def zmatrix_reactant_inchis(ts_zma, frm_keys, brk_keys, remove_stereo=True):
    """ determine the InChIs for the reactants of a TS zmatrix
    """
    rct_ich = zmatrix_reactants_inchi(ts_zma, frm_keys, brk_keys,
                                      remove_stereo=remove_stereo)
    rct_ichs = tuple(map(automol.inchi.recalculate,
                         automol.inchi.split(rct_ich)))
    return rct_ichs


def zmatrix_product_inchis(ts_zma, frm_keys, brk_keys, remove_stereo=True):
    """ determine the InChIs for the products of a TS zmatrix
    """
    prd_ich = zmatrix_products_inchi(ts_zma, frm_keys, brk_keys,
                                     remove_stereo=remove_stereo)
    prd_ichs = tuple(map(automol.inchi.recalculate,
                         automol.inchi.split(prd_ich)))
    return prd_ichs


__all__ = [
    'min_hyd_mig_dist',
    'hydrogen_migration',
    'min_unimolecular_elimination_dist',
    'concerted_unimolecular_elimination',
    'beta_scission',
    'insertion',
    'substitution',
    'addition',
    'hydrogen_abstraction'
]


if __name__ == '__main__':
    # CH4+H
    RCT_ZMAS = [
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
          ('H', (0, 1, 2), ('R4', 'A4', 'D4'))),
         {'R1': 2.063,
          'R2': 2.063, 'A2': 1.9106,
          'R3': 2.063, 'A3': 1.9106, 'D3': 2.0943,
          'R4': 2.063, 'A4': 1.9106, 'D4': 4.1887}),
        ((('H', (None, None, None), (None, None, None)),),
         {}),
    ]
    PRD_ZMAS = [
        ((('C', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None)),
          ('H', (0, 1, None), ('R2', 'A2', None)),
          ('H', (0, 1, 2), ('R3', 'A3', 'D3'))),
         {'R1': 2.045,
          'R2': 2.045, 'A2': 2.0943,
          'R3': 2.045, 'A3': 2.0943, 'D3': 3.1415}),
        ((('H', (None, None, None), (None, None, None)),
          ('H', (0, None, None), ('R1', None, None))),
         {'R1': 1.31906}),
    ]
    TS_ZMA, DIST_NAME, FRM_KEY, BRK_KEY, TORS_NAMES = (
        hydrogen_abstraction(RCT_ZMAS, PRD_ZMAS))
    RCT_GRA = zmatrix_reactant_graph(TS_ZMA, [FRM_KEY], [BRK_KEY])
    PRD_GRA = zmatrix_product_graph(TS_ZMA, [FRM_KEY], [BRK_KEY])
    print('reactant:')
    print(automol.graph.string(RCT_GRA))
    print(zmatrix_reactants_inchi(TS_ZMA, [FRM_KEY], [BRK_KEY]))
    print(zmatrix_reactant_inchis(TS_ZMA, [FRM_KEY], [BRK_KEY]))
    print('product:')
    print(automol.graph.string(PRD_GRA))
    print(zmatrix_products_inchi(TS_ZMA, [FRM_KEY], [BRK_KEY]))
    print(zmatrix_product_inchis(TS_ZMA, [FRM_KEY], [BRK_KEY]))
