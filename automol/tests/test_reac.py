""" test automol.reac
"""

import sys
import automol


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


def test__reac__reactants_graph():
    """ test reac.reactants_graph
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert automol.reac.reactants_graph(rxn) == (
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


def test__reac__products_graph():
    """ test reac.product_graphs
    """
    rxn = automol.reac.from_string(SUBSTITUTION_RXN_STR)
    assert automol.reac.products_graph(rxn) == (
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


def test__reac__hydrogen_migration():
    """ test hydrogen migration functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['C[CH]COO']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R5'
    assert const_names == ('R4',)
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__beta_scission():
    """ test beta scission functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['[O][O]', 'CC[CH2]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D8', 'D11', 'D5'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R8'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [3, 1, 1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__ring_forming_scission():
    """ test ring-forming scission functionality
    """
    rct_smis = ['[CH2]CCCOO']
    prd_smis = ['C1CCCO1', '[OH]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D14'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R13'
    assert const_names == ('A4', 'A7', 'A10', 'D7', 'D10', 'D13')
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__elimination():
    """ test elimination functionality
    """
    rct_smis = ['CCCO[O]']
    prd_smis = ['CC=C', 'O[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R2'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__hydrogen_abstraction():
    """ test hydrogen abstraction functionality
    """
    rct_smis = ['CCO', '[CH3]']
    prd_smis = ['[CH2]CO', 'C']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D11', 'D3', 'D6'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R10'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    print(rxn)
    print(zrxn)
    print(grxn)

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [1, 1, 3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)

    # Extra test cases:
    rxn_smis_lst = [
        (['C(C)(C)C', '[OH]'], ['[C](C)(C)C', 'O']),
        (['C', '[H]'], ['[CH3]', '[H][H]']),
        (['C', '[OH]'], ['[CH3]', 'O']),
        (['[O]O', 'CCC=C[CH]CCCCC'], ['O=O', 'CCCC=CCCCCC']),
    ]
    for rct_smis, prd_smis in rxn_smis_lst:
        rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
        geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

        # reaction object aligned to z-matrix keys
        # (for getting torsion coordinate names)
        zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
        zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

        # You can also do this to determine linear atoms from zmatrix:
        # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
        bnd_keys = automol.reac.rotational_bond_keys(zrxn)
        names = {automol.zmat.torsion_coordinate_name(zma, *k)
                 for k in bnd_keys}
        print(automol.zmat.string(zma, one_indexed=False))
        print(names)

        scan_name = automol.reac.scan_coordinate(zrxn, zma)
        const_names = automol.reac.constraint_coordinates(zrxn, zma)
        print(scan_name)
        print(const_names)

        # graph aligned to geometry keys
        # (for getting rotational groups and symmetry numbers)
        geo, _ = automol.convert.zmat.geometry(zma)
        grxn = automol.reac.relabel_for_geometry(zrxn)
        print(automol.geom.string(geo))

        gbnd_keys = automol.reac.rotational_bond_keys(grxn)

        axes = sorted(map(sorted, gbnd_keys))
        for axis in axes:
            print('axis:', axis)
            groups = automol.reac.rotational_groups(grxn, *axis)
            print('\tgroup 1:', groups[0])
            print('\tgroup 2:', groups[1])
            sym_num = automol.reac.rotational_symmetry_number(grxn, *axis)
            print('\tsymmetry number:', sym_num)


def test__reac__sigma_hydrogen_abstraction():
    """ test sigma hydrogen abstraction functionality
    """
    rct_smis = ['CCO', 'C#[C]']
    prd_smis = ['CC[O]', 'C#C']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k)
             for k in bnd_keys}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, _ = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)

    axes = sorted(map(sorted, gbnd_keys))
    for axis in axes:
        print('axis:', axis)
        groups = automol.reac.rotational_groups(grxn, *axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        sym_num = automol.reac.rotational_symmetry_number(grxn, *axis)
        print('\tsymmetry number:', sym_num)


def test__reac__addition():
    """ test addition functionality
    """
    rct_smis = ['CC[CH2]', '[O][O]']
    prd_smis = ['CCCO[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D11', 'D8', 'D5'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R10'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [3, 1, 1]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__insertion():
    """ test insertion functionality
    """
    rct_smis = ['CC=C', 'O[O]']
    prd_smis = ['CCCO[O]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R3'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__reac__substitution():
    """ test substitution functionality
    """
    rct_smis = ['CO', '[CH2]C']
    prd_smis = ['CCC', '[OH]']

    rxn, rct_geos, _ = _from_smiles(rct_smis, prd_smis)
    geo = automol.reac.ts_geometry(rxn, rct_geos, log=False)

    # reaction object aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.reac.ts_zmatrix(rxn, geo)
    zrxn = automol.reac.relabel_for_zmatrix(rxn, zma_keys, dummy_key_dct)

    # You can also do this to determine linear atoms from zmatrix:
    # bnd_keys = automol.reac.rotational_bond_keys(zrxn, zma=zma)
    bnd_keys = automol.reac.rotational_bond_keys(zrxn)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D3', 'D8', 'D11'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    scan_name = automol.reac.scan_coordinate(zrxn, zma)
    const_names = automol.reac.constraint_coordinates(zrxn, zma)
    assert scan_name == 'R7'
    assert const_names == ()
    print(scan_name)
    print(const_names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.string(geo))

    print(rxn)
    print(zrxn)
    print(grxn)

    # Check that the reaction object can be converted back, if needed
    old_zrxn = zrxn
    zrxn = automol.reac.insert_dummy_atoms(grxn, gdummy_key_dct)
    assert zrxn == old_zrxn

    gbnd_keys = automol.reac.rotational_bond_keys(grxn)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.reac.rotational_groups(grxn, *a) for a in axes]
    sym_nums = [
        automol.reac.rotational_symmetry_number(grxn, *a) for a in axes]
    assert sym_nums == [1, 1, 3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def test__species__demo():
    """ doesn't really belong here, but demonstrates equivalent functionality
    for species
    """
    ich = automol.smiles.inchi('CC#CC#CCCCC#CC')
    geo = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo)

    # graph aligned to z-matrix keys
    # (for getting torsion coordinate names)
    zma, zma_keys, dummy_key_dct = automol.convert.geom.zmatrix(geo)
    gra = automol.graph.relabel_for_zmatrix(gra, zma_keys, dummy_key_dct)

    lin_keys = sorted(
        automol.graph.dummy_atoms_neighbor_atom_key(gra).values())
    bnd_keys = automol.graph.rotational_bond_keys(gra, lin_keys=lin_keys)
    names = {automol.zmat.torsion_coordinate_name(zma, *k) for k in bnd_keys}
    assert names == {'D9', 'D12', 'D15', 'D26'}
    print(automol.zmat.string(zma, one_indexed=False))
    print(names)

    # graph aligned to geometry keys
    # (for getting rotational groups and symmetry numbers)
    geo, gdummy_key_dct = automol.convert.zmat.geometry(zma)
    ggra = automol.graph.relabel_for_geometry(gra)

    # Check that the geometry graph can be converted back, if needed
    zgra = automol.graph.insert_dummy_atoms(ggra, gdummy_key_dct)
    assert zgra == gra

    lin_keys = sorted(gdummy_key_dct.keys())
    gbnd_keys = automol.graph.rotational_bond_keys(ggra, lin_keys=lin_keys)
    assert len(gbnd_keys) == len(bnd_keys)

    axes = sorted(map(sorted, gbnd_keys))
    groups_lst = [automol.graph.rotational_groups(ggra, *a) for a in axes]
    sym_nums = [
        automol.graph.rotational_symmetry_number(ggra, *a, lin_keys=lin_keys)
        for a in axes]
    assert sym_nums == [3, 1, 1, 3]
    for axis, groups, sym_num in zip(axes, groups_lst, sym_nums):
        print('axis:', axis)
        print('\tgroup 1:', groups[0])
        print('\tgroup 2:', groups[1])
        print('\tsymmetry number:', sym_num)


def _from_smiles(rct_smis, prd_smis):
    rct_ichs = list(map(automol.smiles.inchi, rct_smis))
    prd_ichs = list(map(automol.smiles.inchi, prd_smis))

    rct_geos = list(map(automol.inchi.geometry, rct_ichs))
    prd_geos = list(map(automol.inchi.geometry, prd_ichs))

    rct_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, rct_geos)))
    prd_gras = list(map(automol.graph.without_stereo_parities,
                        map(automol.geom.graph, prd_geos)))

    rct_gras, _ = automol.graph.standard_keys_for_sequence(rct_gras)
    prd_gras, _ = automol.graph.standard_keys_for_sequence(prd_gras)

    rxns = automol.reac.find(rct_gras, prd_gras)
    rxn = rxns[0]

    rxn, rct_geos, prd_geos = (
        automol.reac.standard_keys_with_sorted_geometries(
            rxn, rct_geos, prd_geos))
    return rxn, rct_geos, prd_geos


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


# def test__trans():
#     """ test automol.trans.string
#         test automol.trans.from_string
#     """
#
#     ref_trans = (
#         'hydrogen abstraction',
#         frozenset({1, 6}),
#         frozenset({0, 1})
#     )
#
#     # Test relabeling
#     ref_relabel_trans = (
#         'hydrogen abstraction',
#         frozenset({3, 6}),
#         frozenset({0, 3})
#     )
#     relabel_dct = {1: 3, 8: 10}
#     relabel_trans = automol.graph.trans.relabel(ref_trans, relabel_dct)
#     assert automol.graph.trans.from_string(relabel_trans) == ref_relabel_trans
#
#     # Test I/O
#     trans_str = automol.graph.trans.string(ref_trans)
#     assert automol.graph.trans.from_string(trans_str) == ref_trans
#
#     # apply code for rct gra test
#     # ref_gra = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
#     #             3: ('C', 0, None), 4: ('C', 0, None), 5: ('C', 0, None),
#     #             6: ('C', 0, None), 8: ('H', 0, None), 9: ('H', 0, None),
#     #             10: ('H', 0, None), 11: ('H', 0, None), 12: ('H', 0, None),
#     #             13: ('H', 0, None), 14: ('H', 0, None), 15: ('H', 0, None),
#     #             16: ('H', 0, None), 17: ('H', 0, None), 18: ('H', 0, None),
#     #             19: ('H', 0, None), 20: ('H', 0, None), 21: ('H', 0, None),
#     #             22: ('H', 0, None), 7: ('O', 0, None)},
#     #            {frozenset({4, 6}): (1, None), frozenset({21, 6}): (1, None),
#     #             frozenset({0, 2}): (1, None), frozenset({2, 4}): (1, None),
#     #             frozenset({5, 6}): (1, None), frozenset({17, 4}): (1, None),
#     #             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None),
#     #             frozenset({20, 6}): (1, None), frozenset({0, 10}): (1, None),
#     #             frozenset({1, 12}): (1, None), frozenset({2, 14}): (1, None),
#     #             frozenset({18, 5}): (1, None), frozenset({1, 13}): (1, None),
#     #             frozenset({0, 8}): (1, None), frozenset({0, 9}): (1, None),
#     #             frozenset({3, 15}): (1, None), frozenset({1, 11}): (1, None),
#     #             frozenset({19, 5}): (1, None), frozenset({16, 3}): (1, None),
#     #             frozenset({22, 7}): (1, None)})
#
#     # rct_atm_keys_lst = automol.graph.connected_components_atom_keys(rct_gra)
#     # print(rct_atm_keys_lst)
#
#     # # this is how we can get the product graph
#     # prd_gra = automol.graph.trans.apply(tra, rct_gra)
#     # prd_atm_keys_lst = automol.graph.connected_components_atom_keys(prd_gra)
#     # print(prd_atm_keys_lst)


# def test__trans__is_stereo_compatible():
#     """ test graph.trans.is_stereo_compatible
#     """
#     cgr1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
#              3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
#              6: ('O', 1, None)},
#             {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
#              frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
#              frozenset({1, 3}): (1, None)})
#     cgr2 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
#              3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
#              6: ('O', 1, None)},
#             {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
#              frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
#              frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
#
#     cgr1 = automol.graph.explicit(cgr1)
#     cgr2 = automol.graph.explicit(cgr2)
#
#     cgr1s = automol.graph.connected_components(cgr1)
#     cgr2s = automol.graph.connected_components(cgr2)
#
#     ref_compat_lst = (
#         True, False, True, False, True, False, True, False,
#         True, False, True, False, True, False, True, False
#     )
#     tras, _, _ = graph.reac.addition(cgr1s, cgr2s)
#     for tra in tras:
#         assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)
#
#         sgr1 = graph.stereomers(cgr1)[0]
#         for idx, sgr2 in enumerate(graph.stereomers(cgr2)):
#             assert ref_compat_lst[idx] == (
#                 graph.trans.is_stereo_compatible(tra, sgr1, sgr2))


def test__prod__hydrogen_abstraction():
    """ test graph.reac.prod_hydrogen_abstraction
    """

    c4h10_gra = automol.zmat.graph(C4H10_ZMA)
    oh_gra = automol.zmat.graph(OH_ZMA)

    prod_gras = automol.reac.prod_hydrogen_abstraction(c4h10_gra, oh_gra)

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
    """ test graph.reac.prod_hydrogen migration
    """

    ccccch2_gra = automol.zmat.graph(CCCCCH2_ZMA)

    prod_gras = automol.reac.prod_hydrogen_migration(ccccch2_gra)

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


# def test__prod__addition():
#     """ test graph.reac.prod_addition
#     """
# 
#     ch2cch2_gra = automol.zmat.graph(CH2CCH2_ZMA)
#     h_gra = automol.zmat.graph(H_ZMA)
# 
#     prod_gras = automol.reac.prod_addition(ch2cch2_gra, h_gra)
# 
#     assert len(prod_gras) == 2
#     assert all(len(prod_gra) == 1 for prod_gra in prod_gras)
#     assert prod_gras[0] == (
#         ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
#           3: ('H', 0, None), 4: ('X', 0, None), 5: ('C', 0, None),
#           6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None)},
#          {frozenset({5, 6}): (1, None), frozenset({1, 5}): (1, None),
#           frozenset({0, 8}): (1, None), frozenset({0, 3}): (1, None),
#           frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
#           frozenset({5, 7}): (1, None)}),)


def test__prod__beta_scission():
    """ test graph.reac.prod_beta_scission
    """

    ch3ch2ch2o_gra = automol.zmat.graph(CH3CH2CH2O_ZMA)

    prod_gras = automol.reac.prod_beta_scission(ch3ch2ch2o_gra)

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
    """ test graph.reac.prod_homolytic_scission
    """

    c4h10_gra = automol.zmat.graph(C4H10_ZMA)

    prod_gras = automol.reac.prod_homolytic_scission(c4h10_gra)

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


if __name__ == '__main__':
    # test__reac__string()
    # test__reac__hydrogen_migration()
    # test__reac__ring_forming_scission()
    test__reac__hydrogen_abstraction()
    # test__reac__insertion()
    # test__reac__substitution()
    # test__reac__sigma_hydrogen_abstraction()
    # test__reac__hydrogen_migration()
    # test__reac__beta_scission()
    # test__reac__ring_forming_scission()
    # test__reac__elimination()
    # test__reac__addition()
    # test__reac__insertion()
    test__reac__substitution()
