"""
Import all the ts builder functions
"""
import automol.inchi
from automol.zmatrix._unimol_ts import min_hyd_mig_dist
from automol.zmatrix._unimol_ts import hydrogen_migration
from automol.zmatrix._unimol_ts import min_unimolecular_elimination_dist
from automol.zmatrix._unimol_ts import concerted_unimolecular_elimination
from automol.zmatrix._unimol_ts import beta_scission
from automol.zmatrix._bimol_ts import insertion
from automol.zmatrix._bimol_ts import substitution
from automol.zmatrix._bimol_ts import addition
from automol.zmatrix._bimol_ts import hydrogen_abstraction


def zmatrix_reactant_graph(ts_zma, frm_keys, brk_keys):
    """ determine the reactant graph for a TS zmatrix

    (graph keys are aligned with z-matrix rows)
    """
    rct_gra = automol.convert.zmatrix.graph(ts_zma)
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


def zmatrix_reactants_inchi(ts_zma, frm_keys, brk_keys):
    """ determine the InChI for the reactants of a TS zmatrix
    """
    rct_gra = zmatrix_reactant_graph(ts_zma, frm_keys, brk_keys)
    rct_ich = automol.convert.graph.inchi(rct_gra)
    return rct_ich


def zmatrix_products_inchi(ts_zma, frm_keys, brk_keys):
    """ determine the InChI for the products of a TS zmatrix
    """
    prd_gra = zmatrix_product_graph(ts_zma, frm_keys, brk_keys)
    prd_ich = automol.convert.graph.inchi(prd_gra)
    return prd_ich


def zmatrix_reactant_inchis(ts_zma, frm_keys, brk_keys):
    """ determine the InChIs for the reactants of a TS zmatrix
    """
    rct_ich = zmatrix_reactants_inchi(ts_zma, frm_keys, brk_keys)
    rct_ichs = tuple(map(automol.inchi.recalculate,
                         automol.inchi.split(rct_ich)))
    return rct_ichs


def zmatrix_product_inchis(ts_zma, frm_keys, brk_keys):
    """ determine the InChIs for the products of a TS zmatrix
    """
    prd_ich = zmatrix_products_inchi(ts_zma, frm_keys, brk_keys)
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
