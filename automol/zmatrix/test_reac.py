""" test automol.graph
"""

import automol
from automol import graph
from automol.zmatrix._util import shifted_standard_zmas_graphs
import automol.graph.reac_old as reac_old
import automol.graph.trans_old as trans_old


# ZMA Bank
C2H6_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CC')))
C4H10_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCCC')))
C2H4_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('C=C')))
CH4_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('C')))
CH2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]')))
OH_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[OH]')))
H_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[H]')))
CH3_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[CH3]')))
H2O_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('O')))
HO2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('O[O]')))
CH2O_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('C=O')))
CH3CH2O_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CC[O]')))
CH3CH2CH2O_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCC[O]')))
H2O2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('OO')))
CH2COH_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]CO')))
CH3CH2CH2CH2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCC[CH2]')))
CH3CHCH2CH3_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CC[CH]C')))
C3H8_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCC')))
CH3CH2OO_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCO[O]')))
CH2CH2OOH_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('[CH2]COO')))
CCH2OCH2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('C1CO1')))
CCCCCH2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('CCCC[CH2]')))
CH2CCH2_ZMA = automol.convert.geom.zmatrix_x2z_old(
    automol.inchi.geometry(automol.smiles.inchi('C=C=C')))


def test__mult():
    """ test automol.mult.ts.high
        test automol.mult.ts.low
        test automol.mult.spin
    """

    rct_muls = (2, 2)
    prd_muls1 = (1, 1)
    prd_muls2 = (3, 1)
    assert automol.mult.ts.low(rct_muls, prd_muls1) == 1
    assert automol.mult.ts.high(rct_muls, prd_muls2) == 3

    mult = 3
    assert automol.mult.spin(mult) == 2


def test__trans():
    """ test automol.trans.string
        test automol.trans.from_string
    """

    ref_trans = (
        'hydrogen abstraction',
        frozenset({1, 6}),
        frozenset({0, 1})
    )

    # Test relabeling
    ref_relabel_trans = (
        'hydrogen abstraction',
        frozenset({3, 6}),
        frozenset({0, 3})
    )
    relabel_dct = {1: 3, 8: 10}
    relabel_trans = trans_old.relabel(ref_trans, relabel_dct)
    assert trans_old.from_string(relabel_trans) == ref_relabel_trans

    # Test I/O
    trans_str = trans_old.string(ref_trans)
    assert trans_old.from_string(trans_str) == ref_trans


def test__trans__is_stereo_compatible():
    """ test trans_old.is_stereo_compatible
    """
    cgr1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
             frozenset({1, 3}): (1, None)})
    cgr2 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})

    cgr1 = graph.explicit(cgr1)
    cgr2 = graph.explicit(cgr2)

    cgr1s = graph.connected_components(cgr1)
    cgr2s = graph.connected_components(cgr2)

    ref_compat_lst = (
        True, False, True, False, True, False, True, False,
        True, False, True, False, True, False, True, False
    )
    tras, _, _ = reac_old.addition(cgr1s, cgr2s)
    for tra in tras:
        assert graph.backbone_isomorphic(trans_old.apply(tra, cgr1), cgr2)

        sgr1 = graph.stereomers(cgr1)[0]
        for idx, sgr2 in enumerate(graph.stereomers(cgr2)):
            assert ref_compat_lst[idx] == (
                trans_old.is_stereo_compatible(tra, sgr1, sgr2))


def test__prod__hydrogen_abstraction():
    """ test reac_old.prod_hydrogen_abstraction
    """

    c4h10_gra = automol.zmatrix.graph(C4H10_ZMA)
    oh_gra = automol.zmatrix.graph(OH_ZMA)

    prod_gras = reac_old.prod_hydrogen_abstraction(c4h10_gra, oh_gra)

    assert len(prod_gras) == 2
    assert all(len(prod_gra) == 2 for prod_gra in prod_gras)
    assert prod_gras[0] == (
        (({0: ('C', 0, None), 1: ('C', 0, None), 3: ('H', 0, None),
           4: ('H', 0, None), 5: ('C', 0, None), 6: ('H', 0, None),
           7: ('H', 0, None), 8: ('C', 0, None), 9: ('H', 0, None),
           10: ('H', 0, None), 11: ('H', 0, None), 12: ('H', 0, None),
           13: ('H', 0, None)},
          {frozenset({8, 11}): (1, None), frozenset({1, 7}): (1, None),
           frozenset({0, 3}): (1, None), frozenset({0, 1}): (1, None),
           frozenset({0, 4}): (1, None), frozenset({1, 5}): (1, None),
           frozenset({10, 5}): (1, None), frozenset({8, 13}): (1, None),
           frozenset({1, 6}): (1, None), frozenset({8, 12}): (1, None),
           frozenset({9, 5}): (1, None), frozenset({8, 5}): (1, None)}),),
        (({0: ('O', 0, None), 1: ('H', 0, None), 5: ('H', 0, None)},
          {frozenset({0, 1}): (1, None), frozenset({0, 5}): (1, None)}),))


def test__prod__hydrogen_migration():
    """ test reac_old.prod_hydrogen migration
    """

    ccccch2_gra = automol.zmatrix.graph(CCCCCH2_ZMA)

    prod_gras = reac_old.prod_hydrogen_migration(ccccch2_gra)

    assert len(prod_gras) == 2
    assert all(len(prod_gra) == 1 for prod_gra in prod_gras)
    assert prod_gras[0] == (
        ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
          3: ('H', 0, None), 4: ('C', 0, None), 6: ('H', 0, None),
          7: ('C', 0, None), 8: ('H', 0, None), 9: ('H', 0, None),
          10: ('C', 0, None), 11: ('H', 0, None), 12: ('H', 0, None),
          13: ('H', 0, None), 14: ('H', 0, None), 15: ('H', 0, None),
          16: ('H', 0, None)},
         {frozenset({12, 7}): (1, None), frozenset({1, 4}): (1, None),
          frozenset({10, 15}): (1, None), frozenset({10, 7}): (1, None),
          frozenset({0, 3}): (1, None), frozenset({0, 1}): (1, None),
          frozenset({0, 2}): (1, None), frozenset({0, 16}): (1, None),
          frozenset({10, 13}): (1, None), frozenset({10, 14}): (1, None),
          frozenset({9, 4}): (1, None), frozenset({1, 6}): (1, None),
          frozenset({11, 7}): (1, None), frozenset({8, 4}): (1, None),
          frozenset({4, 7}): (1, None)}),)


def test__prod__addition():
    """ test reac_old.prod_addition
    """

    ch2cch2_gra = automol.zmatrix.graph(CH2CCH2_ZMA)
    h_gra = automol.zmatrix.graph(H_ZMA)

    prod_gras = reac_old.prod_addition(ch2cch2_gra, h_gra)

    assert len(prod_gras) == 2
    assert all(len(prod_gra) == 1 for prod_gra in prod_gras)
    assert prod_gras[0] == (
        ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
          3: ('H', 0, None), 4: ('X', 0, None), 5: ('C', 0, None),
          6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None)},
         {frozenset({5, 6}): (1, None), frozenset({1, 5}): (1, None),
          frozenset({0, 8}): (1, None), frozenset({0, 3}): (1, None),
          frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
          frozenset({5, 7}): (1, None)}),)


def test__prod__beta_scission():
    """ test reac_old.prod_beta_scission
    """

    ch3ch2ch2o_gra = automol.zmatrix.graph(CH3CH2CH2O_ZMA)

    prod_gras = reac_old.prod_beta_scission(ch3ch2ch2o_gra)

    assert len(prod_gras) == 2
    assert all(len(prod_gra) == 2 for prod_gra in prod_gras)
    assert prod_gras[0] == (
        ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
          3: ('H', 0, None), 4: ('H', 0, None), 6: ('H', 0, None),
          7: ('H', 0, None)},
         {frozenset({1, 7}): (1, None), frozenset({0, 3}): (1, None),
          frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
          frozenset({0, 4}): (1, None), frozenset({1, 6}): (1, None)}),
        ({5: ('C', 0, None), 8: ('O', 0, None), 9: ('H', 0, None),
          10: ('H', 0, None)},
         {frozenset({8, 5}): (1, None), frozenset({10, 5}): (1, None),
          frozenset({9, 5}): (1, None)}))


def test__prod__homolytic_scission():
    """ test reac_old.prod_homolytic_scission
    """

    c4h10_gra = automol.zmatrix.graph(C4H10_ZMA)

    prod_gras = reac_old.prod_homolytic_scission(c4h10_gra)

    assert len(prod_gras) == 4
    assert all(len(prod_gra) == 2 for prod_gra in prod_gras)
    assert prod_gras[0] == (
        ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
          3: ('H', 0, None), 4: ('H', 0, None), 5: ('C', 0, None),
          6: ('H', 0, None), 7: ('H', 0, None), 8: ('C', 0, None),
          9: ('H', 0, None), 10: ('H', 0, None), 12: ('H', 0, None),
          13: ('H', 0, None)},
         {frozenset({1, 7}): (1, None), frozenset({0, 3}): (1, None),
          frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
          frozenset({0, 4}): (1, None), frozenset({1, 5}): (1, None),
          frozenset({10, 5}): (1, None), frozenset({8, 13}): (1, None),
          frozenset({1, 6}): (1, None), frozenset({8, 12}): (1, None),
          frozenset({9, 5}): (1, None), frozenset({8, 5}): (1, None)}),
        ({11: ('H', 0, None)}, {}))


# BIMOL TS
def test__ts_hydrogen_abstraction():
    """ test zmatrix.ts.hydrogen_abstraction
    """

    rct_zmas = [CH4_ZMA, OH_ZMA]
    prd_zmas = [CH3_ZMA, H2O_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.hydrogen_abstraction(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


def test__ts_addition():
    """ test zmatrix.ts.addition
    """

    rct_zmas = [C2H4_ZMA, OH_ZMA]
    prd_zmas = [CH2COH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.addition(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


def test__ts_insertion():
    """ test zmatrix.ts.insertion
    """

    rct_zmas = [C2H6_ZMA, CH2_ZMA]
    prd_zmas = [C3H8_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.insertion(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


def test__ts_substitution():
    """ test zmatrix.ts.substitution
    """

    rct_zmas = [H2O2_ZMA, H_ZMA]
    prd_zmas = [H2O_ZMA, OH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.substitution(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])


# UNIMOL TS
def test__ts_hydrogen_migration():
    """ test zmatrix.ts.hydrogen_migration
    """

    rct_zmas = [CH3CH2OO_ZMA]
    prd_zmas = [CH2CH2OOH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.hydrogen_migration(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


def test__ts_beta_scission():
    """ test zmatrix.ts.beta_scission
    """

    rct_zmas = [CH3CH2O_ZMA]
    prd_zmas = [CH3_ZMA, CH2O_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.beta_scission(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


def test__ts_elimination():
    """ test zmatrix.ts.elimination
    """

    rct_zmas = [CH3CH2OO_ZMA]
    prd_zmas = [C2H4_ZMA, HO2_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.concerted_unimol_elimination(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


def test__ts_ring_forming_scission():
    """ test zmatrix.ts.ring_forming_scission
    """

    rct_zmas = [CH2CH2OOH_ZMA]
    prd_zmas = [CCH2OCH2_ZMA, OH_ZMA]
    rct_zmas, rct_gras = shifted_standard_zmas_graphs(
        rct_zmas, remove_stereo=True)
    prd_zmas, prd_gras = shifted_standard_zmas_graphs(
        prd_zmas, remove_stereo=True)
    rct_zmas = [rct_zmas]
    prd_zmas = [prd_zmas]

    tras, _, _ = reac_old.classify(rct_gras, prd_gras)
    rtyp = trans_old.reaction_class(tras[0])
    print('\nrtyp', rtyp)

    zma_ret = automol.zmatrix.ts.ring_forming_scission(
        rct_zmas, prd_zmas, tras)
    print('zma\n', automol.zmatrix.string(zma_ret['ts_zma']))
    print('bnd keys\n', zma_ret['bnd_keys'])
    print('const keys\n', zma_ret['const_keys'])

    zma = zma_ret['ts_zma']
    geo = automol.convert.zmatrix.geometry(zma)
    print(automol.geom.string(geo))


if __name__ == '__main__':
    # test__mult()
    # test__trans()
    # test__trans__is_stereo_compatible()
    # test__prod__hydrogen_abstraction()
    # test__prod__hydrogen_migration()
    # test__prod__addition()
    # test__prod__beta_scission()
    # test__prod__homolytic_scission()
    # BIMOL
    # test__ts_hydrogen_abstraction()
    # test__ts_addition()
    # test__ts_substitution()
    # test__ts_insertion()
    # UNIMOL
    # test__ts_hydrogen_migration()
    # test__ts_beta_scission()
    test__ts_elimination()
    # test__ts_ring_forming_scission()
