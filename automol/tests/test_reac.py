""" test automol.reac
"""

import numpy
import automol
from automol.par import ReactionClass


SUBSTITUTION_RXN_STR = """
reaction class: substitution
forward TS atoms:
  1: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 0.9, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 0, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 0.1, stereo_parity: null}
  8-9: {order: 1, stereo_parity: null}
  8-10: {order: 1, stereo_parity: null}
  8-11: {order: 1, stereo_parity: null}
  9-12: {order: 1, stereo_parity: null}
  9-13: {order: 1, stereo_parity: null}
  9-14: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7]
- [8, 9, 10, 11, 12, 13, 14]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-3: {order: 0.9, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 1, stereo_parity: null}
  1-6: {order: 1, stereo_parity: null}
  1-12: {order: 0.1, stereo_parity: null}
  2-3: {order: 1, stereo_parity: null}
  2-7: {order: 1, stereo_parity: null}
  2-8: {order: 1, stereo_parity: null}
  2-9: {order: 1, stereo_parity: null}
  3-10: {order: 1, stereo_parity: null}
  3-11: {order: 1, stereo_parity: null}
  12-13: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
- [12, 13]
"""

MIGRATION_RXN_STR = """
reaction class: hydrogen migration
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-7: {order: 1, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  5-6: {order: 0.1, stereo_parity: null}
  6-7: {order: 0.9, stereo_parity: null}
  7-8: {order: 1, stereo_parity: null}
  7-9: {order: 1, stereo_parity: null}
  8-10: {order: 1, stereo_parity: null}
  8-11: {order: 1, stereo_parity: null}
  8-12: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
backward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-6: {order: 1, stereo_parity: null}
  1-7: {order: 1, stereo_parity: null}
  1-8: {order: 1, stereo_parity: null}
  2-3: {order: 1, stereo_parity: null}
  2-9: {order: 1, stereo_parity: null}
  2-12: {order: 0.1, stereo_parity: null}
  3-5: {order: 1, stereo_parity: null}
  3-10: {order: 1, stereo_parity: null}
  3-11: {order: 1, stereo_parity: null}
  4-5: {order: 1, stereo_parity: null}
  4-12: {order: 0.9, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
"""

# ZMA Bank
C4H10_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCCC')))
OH_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[OH]')))
H_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('[H]')))
CCCCCH2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCCC[CH2]')))
CH2CCH2_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('C=C=C')))
CH3CH2CH2O_ZMA = automol.geom.zmatrix(
    automol.inchi.geometry(automol.smiles.inchi('CCC[O]')))


def test__reac__string():
    """ test reac.string
    """
    rxn_str = SUBSTITUTION_RXN_STR
    rxn = automol.reac.from_string(rxn_str)
    assert automol.reac.string(rxn).strip() == rxn_str.strip()


def test__reac__forming_bond_keys():
    """ test reac.forming_bond_keys
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert (automol.reac.forming_bond_keys(rxn) ==
            frozenset({frozenset({1, 7})}))
    assert (automol.reac.forming_bond_keys(rxn, rev=True) ==
            frozenset({frozenset({0, 11})}))


def test__reac__breaking_bond_keys():
    """ test reac.breaking_bond_keys
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert (automol.reac.breaking_bond_keys(rxn) ==
            frozenset({frozenset({0, 1})}))
    assert (automol.reac.breaking_bond_keys(rxn, rev=True) ==
            frozenset({frozenset({0, 2})}))


def test__reac__forming_rings_atom_keys():
    """ test reac.forming_rings_atom_keys
    """
    rxn = automol.reac.from_string(MIGRATION_RXN_STR)
    assert automol.reac.forming_rings_atom_keys(rxn) == (
        (0, 1, 4, 5, 6),
    )
    assert automol.reac.forming_rings_atom_keys(rxn, rev=True) == (
        (1, 2, 4, 3, 11),
    )


def test__reac__forming_rings_bond_keys():
    """ test reac.forming_rings_bond_keys
    """
    rxn = automol.reac.from_string(MIGRATION_RXN_STR)
    assert automol.reac.forming_rings_bond_keys(rxn) == (
        frozenset({frozenset({1, 4}), frozenset({0, 6}), frozenset({4, 5}),
                   frozenset({0, 1}), frozenset({5, 6})}),
    )
    assert automol.reac.forming_rings_bond_keys(rxn, rev=True) == (
        frozenset({frozenset({3, 4}), frozenset({1, 2}), frozenset({1, 11}),
                   frozenset({2, 4}), frozenset({11, 3})}),
    )


def test__reac__breaking_rings_atom_keys():
    """ test reac.breaking_rings_atom_keys
    """
    rxn = automol.reac.from_string(MIGRATION_RXN_STR)
    assert automol.reac.breaking_rings_atom_keys(rxn) == (
        (0, 1, 4, 5, 6),
    )
    assert automol.reac.breaking_rings_atom_keys(rxn, rev=True) == (
        (1, 2, 4, 3, 11),
    )


def test__reac__breaking_rings_bond_keys():
    """ test reac.breaking_rings_bond_keys
    """
    rxn = automol.reac.from_string(MIGRATION_RXN_STR)
    assert automol.reac.breaking_rings_bond_keys(rxn) == (
        frozenset({frozenset({1, 4}), frozenset({0, 6}), frozenset({4, 5}),
                   frozenset({0, 1}), frozenset({5, 6})}),
    )
    assert automol.reac.breaking_rings_bond_keys(rxn, rev=True) == (
        frozenset({frozenset({3, 4}), frozenset({1, 2}), frozenset({1, 11}),
                   frozenset({2, 4}), frozenset({11, 3})}),
    )


def test__reac__reactant_graphs():
    """ test reac.reactant_graphs
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert automol.reac.reactant_graphs(rxn) == (
        ({0: ('O', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
          3: ('X', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
          6: ('H', 0, None)},
         {frozenset({1, 4}): (1, None), frozenset({0, 1}): (1, None),
          frozenset({0, 2}): (1, None), frozenset({1, 5}): (1, None),
          frozenset({1, 6}): (1, None), frozenset({1, 3}): (0, None)}),
        ({7: ('C', 0, None), 8: ('C', 0, None), 9: ('H', 0, None),
          10: ('H', 0, None), 11: ('H', 0, None), 12: ('H', 0, None),
          13: ('H', 0, None)},
         {frozenset({8, 11}): (1, None), frozenset({10, 7}): (1, None),
          frozenset({9, 7}): (1, None), frozenset({8, 7}): (1, None),
          frozenset({8, 13}): (1, None), frozenset({8, 12}): (1, None)})
    )


def test__reac__product_graphs():
    """ test reac.product_graphs
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert automol.reac.product_graphs(rxn) == (
        ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
          3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
          6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
          9: ('H', 0, None), 10: ('H', 0, None)},
         {frozenset({1, 7}): (1, None), frozenset({10, 2}): (1, None),
          frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
          frozenset({0, 2}): (1, None), frozenset({0, 4}): (1, None),
          frozenset({0, 5}): (1, None), frozenset({8, 1}): (1, None),
          frozenset({1, 6}): (1, None), frozenset({9, 2}): (1, None)}),
        ({11: ('O', 0, None), 12: ('H', 0, None)},
         {frozenset({11, 12}): (1, None)})
    )


def test__reac__reagents_graph():
    """ test reac.reactants_graph
        test products_graph
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)

    rcts_gra = automol.reac.reactants_graph(rxn)
    prds_gra = automol.reac.products_graph(rxn)

    assert rcts_gra == (
        {0: ('O', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
         3: ('X', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
         6: ('H', 0, None), 7: ('C', 0, None), 8: ('C', 0, None),
         9: ('H', 0, None), 10: ('H', 0, None), 11: ('H', 0, None),
         12: ('H', 0, None), 13: ('H', 0, None)},
        {frozenset({1, 4}): (1, None), frozenset({8, 11}): (1, None),
         frozenset({10, 7}): (1, None), frozenset({0, 1}): (1, None),
         frozenset({0, 2}): (1, None), frozenset({9, 7}): (1, None),
         frozenset({8, 7}): (1, None), frozenset({1, 5}): (1, None),
         frozenset({8, 13}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({1, 3}): (0, None), frozenset({8, 12}): (1, None)}
    )

    assert prds_gra == (
        {0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
         3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
         6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
         9: ('H', 0, None), 10: ('H', 0, None), 11: ('O', 0, None),
         12: ('H', 0, None)},
        {frozenset({1, 7}): (1, None), frozenset({10, 2}): (1, None),
         frozenset({1, 2}): (1, None), frozenset({0, 3}): (1, None),
         frozenset({11, 12}): (1, None), frozenset({0, 2}): (1, None),
         frozenset({0, 4}): (1, None), frozenset({0, 5}): (1, None),
         frozenset({8, 1}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({9, 2}): (1, None)})

    assert prds_gra == automol.reac.reactants_graph(rxn, rev=True)


def test__reac__hydrogen_migration():
    """ test hydrogen migration functionality
    """

    rct_smis = ['CCCO[O]']
    prd_smis = ['C[CH]COO']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R2',)
    ref_constraint_dct = {'R1': 2.65}
    ref_scan_grid = (numpy.array([
       3.77945225, 3.66829189, 3.55713153, 3.44597117, 3.33481081,
       3.22365045, 3.11249009, 3.00132973, 2.89016937, 2.77900901,
       2.66784865, 2.55668829, 2.44552793, 2.33436757, 2.22320721,
       2.11204685, 2.00088649, 1.88972613]),)
    ref_update_guess = True
    ref_tors_names = {'D9'}
    ref_tors_symms = [3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac__2ts_hydrogen_migration():
    """ test hydrogen migration functionality

        EXPAND OT GET ALL OF THE STUFF NEEDED
    """

    rct_smis = ['CCC[CH2]']
    prd_smis = ['CC[CH]C']

    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    # Deal with rxn object 1
    rxn1, ts_geo1, _, _ = rxn_objs[0]
    zma1, zma_keys1, dummy_key_dct1 = automol.reac.ts_zmatrix(rxn1, ts_geo1)
    zrxn1 = automol.reac.relabel_for_zmatrix(rxn1, zma_keys1, dummy_key_dct1)
    print(zrxn1)

    bnd_keys1 = automol.reac.rotational_bond_keys(zrxn1)
    names1 = {
        automol.zmat.torsion_coordinate_name(zma1, *k) for k in bnd_keys1}
    print(names1)

    scan_name1 = automol.reac.scan_coordinate(zrxn1, zma1)
    const_names1 = automol.reac.constraint_coordinates(zrxn1, zma1)
    print(scan_name1)
    print(const_names1)

    # Deal with rxn object 2
    rxn2, ts_geo2, _, _ = rxn_objs[1]
    zma2, zma_keys2, dummy_key_dct2 = automol.reac.ts_zmatrix(rxn2, ts_geo2)
    zrxn2 = automol.reac.relabel_for_zmatrix(rxn2, zma_keys2, dummy_key_dct2)

    bnd_keys2 = automol.reac.rotational_bond_keys(zrxn2)
    names2 = {
        automol.zmat.torsion_coordinate_name(zma2, *k) for k in bnd_keys2}
    print(names2)

    scan_name2 = automol.reac.scan_coordinate(zrxn2, zma2)
    const_names2 = automol.reac.constraint_coordinates(zrxn2, zma2)
    print(scan_name2)
    print(const_names2)


def test__reac__beta_scission():
    """ test beta scission functionality
    """

    rct_smis = ['CCCO[O]']
    prd_smis = ['[O][O]', 'CC[CH2]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R8',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.89128097, 2.99303546, 3.09478994, 3.19654442, 3.29829891,
        3.40005339, 3.50180787, 3.60356236, 3.70531684, 3.80707133,
        3.90882581, 4.01058029, 4.11233478, 4.21408926]),)
    ref_update_guess = False
    ref_tors_names = {'D11', 'D5', 'D8'}
    ref_tors_symms = [3, 1, 1]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac__ring_forming_scission():
    """ test ring-forming scission functionality
    """

    rct_smis = ['[CH2]CCCOO']
    prd_smis = ['C1CCCO1', '[OH]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R13',)
    ref_constraint_dct = {'A4': 1.85, 'A7': 2.15, 'A10': 1.91,
                          'D7': 6.28, 'D10': 0.01, 'D13': 3.15}
    ref_scan_grid = ((
        2.834589188186742, 3.0235618007325247, 3.212534413278308,
        3.4015070258240905, 3.590479638369873, 3.7794522509156563,
        3.968424863461439),)
    ref_update_guess = False
    ref_tors_names = {'D14'}
    ref_tors_symms = [1]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac__elimination():
    """ test elimination functionality
    """

    rct_smis = ['CCCCO[O]']
    prd_smis = ['CCC=C', 'O[O]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R2', 'R3')
    ref_constraint_dct = None
    ref_scan_grid = (
        numpy.array([
            2.17318504, 2.49713809, 2.82109114, 3.14504419,
            3.46899724, 3.79295029, 4.11690334, 4.44085639]),
        numpy.array([3.08025358, 3.45819881, 3.83614403, 4.21408926]))
    ref_update_guess = False  # correct?
    ref_tors_names = {'D9', 'D12'}
    ref_tors_symms = [1, 3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)

    # Extra test cases:
    rxn_smis_lst = [
        # HONO elim.; use dbl bnd to make the forming bond in TS
        (['CCCON(=O)=O'], ['CCC=O', 'N(=O)O']),
        # CH2 elim.; 3-member ring in TS (forms C-C, not C-H as it should?)
        (['CCC'], ['CC', '[CH2]']),
        # H2 elim.; 3-member ring in TS
        (['C=O'], ['[C-]#[O+]', '[HH]'])
    ]
    for rct_smis, prd_smis in rxn_smis_lst:
        print('\n\nRXN ID FOR', rct_smis, prd_smis)
        rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
        print(rxn_objs)
        _check_reaction(rxn_objs[0])


def test__reac__hydrogen_abstraction():
    """ test hydrogen abstraction functionality
    """

    rct_smis = ['CCO', '[CH3]']
    prd_smis = ['[CH2]CO', 'C']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R10',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.24877409, 2.49173888, 2.73470366, 2.97766845, 3.22063324,
        3.46359803, 3.70656281, 3.9495276]),)
    ref_update_guess = False
    ref_tors_names = {'D3', 'D11', 'D6'}
    ref_tors_symms = [1, 1, 3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)

    # Extra test cases:
    rxn_smis_lst = [
        (['C(C)(C)C', '[OH]'], ['[C](C)(C)C', 'O']),
        (['C', '[H]'], ['[CH3]', '[H][H]']),
        (['C', '[OH]'], ['[CH3]', 'O']),
        (['CC', '[H]'], ['C[CH2]', '[H][H]']),
    ]
    for rct_smis, prd_smis in rxn_smis_lst:
        rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
        _check_reaction(rxn_objs[0])


def test__reac__sigma_hydrogen_abstraction():
    """ test sigma hydrogen abstraction functionality
    """

    rct_smis = ['CCO', 'C#[C]']
    prd_smis = ['CC[O]', 'C#C']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R10',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.24877409, 2.49173888, 2.73470366, 2.97766845, 3.22063324,
        3.46359803, 3.70656281, 3.9495276]),)
    ref_update_guess = False
    ref_tors_names = {'D8', 'D5'}
    ref_tors_symms = [3, 1]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac__addition():
    """ test addition functionality
    """

    rct_smis = ['CC[CH2]', '[O][O]']
    prd_smis = ['CCCO[O]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R10',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.89128097, 2.94128097, 2.99628097, 3.05678097, 3.12333097,
        3.19653597, 3.27706147, 3.36563952, 3.46307538, 3.57025482,
        3.6881522, 3.81783933, 3.96049516, 4.11741658, 4.29003014]),)
    ref_update_guess = False
    ref_tors_names = {'D7', 'D11', 'D4'}
    ref_tors_symms = [1, 1, 3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)

    # Extra test cases:
    rxn_smis_lst = [
        (['C=CCCCCCC', '[CH2]C'], ['CCC[CH]CCCCCC']),
    ]
    for rct_smis, prd_smis in rxn_smis_lst:
        rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
        _check_reaction(rxn_objs[0])


def test__reac__radrad_addition():
    """ test addition functionality
    """

    rct_smis = ['CC[CH2]', '[H]']
    prd_smis = ['CCC']

    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R10',)
    ref_constraint_dct = None
    ref_scan_grid = None  # not calling the radrad builder
    ref_update_guess = False
    ref_tors_names = {'D7', 'D4'}
    ref_tors_symms = [3, 3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac__isc_addition():
    """ test addition functionality
    """

    rct_smis = ['N#N', '[O]']
    prd_smis = ['[N-]=[N+]=O']

    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
    print(rxn_objs)


def test__reac__radrad_hydrogen_abstraction():
    """ test addition functionality
    """

    rct_smis = ['CCC', '[H]']
    prd_smis = ['CC[CH2]', '[HH]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R12',)
    ref_constraint_dct = None
    ref_scan_grid = None  # not calling the radrad builder
    ref_update_guess = False
    ref_tors_names = {'D8', 'D5'}
    ref_tors_symms = [3, 1]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def __reac__insertion():
    """ test insertion functionality
    """

    rct_smis = ['CC=C', 'O[O]']
    prd_smis = ['CCCO[O]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R3',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.05980148, 2.23617592, 2.41255035, 2.58892479, 2.76529923,
        2.94167367, 3.11804811, 3.29442255, 3.47079698, 3.64717142,
        3.82354586, 3.9999203, 4.17629474, 4.35266918, 4.52904361,
        4.70541805]),)
    ref_update_guess = False
    ref_tors_names = {'D9'}
    ref_tors_symms = [3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)

    # Extra test cases:
    rxn_smis_lst = [
        (['CC', '[CH2]'], ['CCC']),
    ]
    for rct_smis, prd_smis in rxn_smis_lst:
        rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
        _check_reaction(rxn_objs[0])


def test__reac__substitution():
    """ test substitution functionality
    """

    rct_smis = ['CO', '[CH2]C']
    prd_smis = ['CCC', '[OH]']
    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)

    ref_scan_names = ('R7',)
    ref_constraint_dct = None
    ref_scan_grid = (numpy.array([
        2.91017823, 3.1136872, 3.31719617, 3.52070514, 3.7242141,
        3.92772307, 4.13123204, 4.334741, 4.53824997, 4.74175894,
        4.94526791, 5.14877687, 5.35228584, 5.55579481]),)
    ref_scan_grid = ()
    ref_update_guess = False
    ref_tors_names = {'D8', 'D3', 'D11'}
    ref_tors_symms = [1, 1, 3]

    _check_reaction(rxn_objs[0],
                    ref_scan_names, ref_constraint_dct,
                    ref_scan_grid, ref_update_guess,
                    ref_tors_names, ref_tors_symms)


def test__reac_util():
    """ test if the internal converter in the reac.util functions work
    """

    rct_smis = ['CC', '[H]']
    prd_smis = ['C[CH2]', '[HH]']

    rxn_objs = automol.reac.rxn_objs_from_smiles(
        rct_smis, prd_smis)
    rxn, geo, _, _ = rxn_objs[0]
    _, zma_keys1, dummy_key_dct1 = automol.reac.ts_zmatrix(rxn, geo)
    zrxn1 = automol.reac.relabel_for_zmatrix(rxn, zma_keys1, dummy_key_dct1)

    zrxn_objs = automol.reac.rxn_objs_from_smiles(
        rct_smis, prd_smis, indexing='zma')
    zrxn2, _, _, _ = zrxn_objs[0]

    assert zrxn1 == zrxn2


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


def test__stereo():
    """ test stereo functionality
    """

    # example 1
    rct_smis = ['FC=C(C(O)F)C(O)F', '[OH]']
    prd_smis = ['FC(O)[C](C(O)F)C(O)F']

    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    assert len(srxns) == 16
    print("Complete stereo expansion for the reaction:")
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_ichs = list(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = list(map(automol.graph.stereo_inchi, prd_gras))
        print(rct_ichs)
        print(prd_ichs)
        print()

    # Assign reactant and product stereo from geometries.
    srxn = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)
    # Note that the original stereo assignments from the product geometries
    # could be inconsistent with the reactant stereo assignments.
    print('Consistent?', automol.reac.is_stereo_consistent(srxn))

    # Add stereo from geometries and expand stereo possibilities consistent
    # with the reactants.
    srxns = automol.reac.expand_product_stereo(srxn)
    print(len(srxns))
    assert len(srxns) == 2
    print("Product expansion for reactant geometry stereo assignments:")
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_ichs = list(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = list(map(automol.graph.stereo_inchi, prd_gras))
        print(rct_ichs)
        print(prd_ichs)
        print()

    # example 2
    rct_smis = ['FC=CC=CF', '[OH]']
    prd_smis = ['FC=C[CH]C(O)F']
    print("Reaction:", rct_smis, "=>", prd_smis)

    rxn_objs = automol.reac.rxn_objs_from_smiles(rct_smis, prd_smis)
    rxn, _, rct_geos, prd_geos = rxn_objs[0]

    # Complete stereo expansion for the reaction
    srxns = automol.reac.expand_stereo(rxn)
    print(len(srxns))
    assert len(srxns) == 16
    print("Complete stereo expansion for the reaction:")
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_ichs = list(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = list(map(automol.graph.stereo_inchi, prd_gras))
        print(rct_ichs)
        print(prd_ichs)
        print()

    # Assign reactant and product stereo from geometries.
    srxn = automol.reac.add_stereo_from_geometries(rxn, rct_geos, prd_geos)
    # Note that the original stereo assignments from the product geometries
    # could be inconsistent with the reactant stereo assignments.
    print('Consistent?', automol.reac.is_stereo_consistent(srxn))
    # Add stereo from geometries and expand stereo possibilities consistent
    # with the reactants.
    srxns = automol.reac.expand_product_stereo(srxn)
    print(len(srxns))
    assert len(srxns) == 4
    print("Product expansion for reactant geometry stereo assignments:")
    for srxn in srxns:
        rct_gras = automol.reac.reactant_graphs(srxn)
        prd_gras = automol.reac.product_graphs(srxn)
        rct_ichs = list(map(automol.graph.stereo_inchi, rct_gras))
        prd_ichs = list(map(automol.graph.stereo_inchi, prd_gras))
        print(rct_ichs)
        print(prd_ichs)
        print()


def test__prod__hydrogen_migration():
    """ test hydrogen migration product enumeration
    """
    rct_gras = _gras_for_prod_tests(['C=CCC[CH2]'])
    rclass = ReactionClass.Typ.HYDROGEN_MIGRATION
    nprods = 7
    _check_products(rct_gras, rclass, nprods)


def __prod__homolytic_scission():
    """ test homolytic scission product enumeration
    """
    rct_gras = _gras_for_prod_tests(['CCCl'])
    rclass = ReactionClass.Typ.HOMOLYT_SCISSION
    nprods = 1
    _check_products(rct_gras, rclass, nprods)
    # check fails because some reactions ID'd as a beta scission


def test__prod__beta_scission():
    """ test beta scission product enumeration
    """
    rct_gras = _gras_for_prod_tests(['C=C[CH]CC'])
    rclass = ReactionClass.Typ.BETA_SCISSION
    nprods = 1
    _check_products(rct_gras, rclass, nprods)


def test__prod__ring_forming_scission():
    """ test ring-forming scission product enumeration
    """
    rct_gras = _gras_for_prod_tests(['CC(OO)CC(OO)C[CH2]'])
    rclass = ReactionClass.Typ.RING_FORM_SCISSION
    nprods = 2
    _check_products(rct_gras, rclass, nprods)


def test__prod__elimination():
    """ test elimination product enumeration
    """
    rct_gras = _gras_for_prod_tests(['CCCO[O]'])
    rclass = ReactionClass.Typ.ELIMINATION
    nprods = 8
    _check_products(rct_gras, rclass, nprods)


def test__prod__hydrogen_abstraction():
    """ test hydrogen abstraction product enumeration
    """
    rct_gras = _gras_for_prod_tests(['CC(=O)C', '[CH3]'])
    rclass = ReactionClass.Typ.HYDROGEN_ABSTRACTION
    nprods = 1
    _check_products(rct_gras, rclass, nprods)


def test__prod__addition():
    """ test addition product enumeration
    """
    rct_gras = _gras_for_prod_tests(['C=CC=C', '[CH3]'])
    rclass = ReactionClass.Typ.ADDITION
    nprods = 2
    _check_products(rct_gras, rclass, nprods)


def test__prod__insertion():
    """ test insertion product enumeration
    """

    rct_gras = _gras_for_prod_tests(['CC=C', 'O[O]'])
    rclass = ReactionClass.Typ.INSERTION
    nprods = 4
    _check_products(rct_gras, rclass, nprods)


# Utility functions for building information
def _gras_for_prod_tests(rct_smis):
    """ Get reactant graphs from smiles
    """

    rct_ichs = list(map(automol.smiles.inchi, rct_smis))
    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    rct_gras = tuple(map(automol.geom.connectivity_graph, rct_geos))
    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)

    return rct_gras


# Checker functions for assessing if tests output correct information
def _check_reaction(rxn_obj,
                    ref_scan_names=None, ref_constraint_dct=None,
                    ref_scan_grid=None, ref_update_guess=None,
                    ref_tors_names=None, ref_tors_symms=None):
    """ Check if all of the information for reactions is correct
    """

    # Unpack the reaction object
    rxn, geo, _, _ = rxn_obj

    # Build Reaction object aligned to z-matrix keys
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # Get scan information
    scan_info = automol.reac.build_scan_info(zrxn, zma)
    scan_names, constraint_dct, scan_grid, update_guess = scan_info
    print('scan grid', scan_grid)
    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.zmat.geometry_with_conversion_info(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)

    # Get torsion information
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    tors_names = {automol.zmat.torsion_coordinate_name(zma, *k)
                  for k in bnd_keys}

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    tors_symms = [automol.reac.rotational_symmetry_number(grxn, *a)
                  for a in axes]

    # Check that the information is correct, requested
    if ref_scan_names is not None:
        assert scan_names == ref_scan_names
    if ref_constraint_dct is not None:
        assert set(constraint_dct.keys()) == set(ref_constraint_dct.keys())
    if ref_scan_grid is not None:
        for rgrd, grd in zip(ref_scan_grid, scan_grid):
            # assert numpy.allclose(rgrd, grd)  correct?
            if rxn.class_ != 'elimination':
                assert numpy.allclose(rgrd, grd)
            else:
                for sub_rgrd, sub_grd in zip(rgrd, grd):
                    assert numpy.allclose(sub_rgrd, sub_grd)
    if ref_update_guess is not None:
        assert update_guess == ref_update_guess
    if ref_tors_names is not None:
        assert tors_names == ref_tors_names
    if ref_tors_names is not None:
        assert tors_symms == ref_tors_symms

    # Check that zrxn -> grxn -> zrxn conversion holds
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn


def _check_products(rct_gras, rxn_class_typ, num_rxns):
    """ Check the products
    """

    # Enumerate all possible reactions, but select the insertions
    rxns = [r for r in automol.reac.enumerate_reactions(rct_gras)
            if r.class_ == rxn_class_typ]
    print('PRODUCTS FOR {}'.format(rxn_class_typ))
    print('num prods\n', len(rxns))

    assert rxns
    # assert len(rxns) == num_rxns

    # Verify the enumerated reactions with the classifier
    for rxn in rxns:
        rct_gras_ = automol.reac.reactant_graphs(rxn)
        prd_gras_ = automol.reac.product_graphs(rxn)
        for gra in prd_gras_:
            print(automol.geom.string(
                automol.inchi.geometry(automol.graph.inchi(gra))))
            print('')
        print('\n\n')
        assert rct_gras_ == rct_gras
        rxns_ = automol.reac.find(rct_gras_, prd_gras_)
        assert any(r.class_ == rxn_class_typ for r in rxns_)


if __name__ == '__main__':
    # test__reac__hydrogen_migration()
    # test__reac__beta_scission()
    # test__reac__ring_forming_scission()
    # test__reac__hydrogen_abstraction()
    # test__reac__sigma_hydrogen_abstraction()
    # test__reac__addition()
    # test__reac__radrad_addition()
    # test__reac__isc_addition()
    # test__reac__radrad_hydrogen_abstraction()
    # test__reac__insertion()
    # test__reac__substitution()
    # test__prod__homolytic_scission()
    test__prod__beta_scission()
    # test__prod__ring_forming_scission()
