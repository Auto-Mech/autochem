""" test automol.graph.ts
"""
import numpy
from automol import graph


# CCO[C@H](O[O])C => C[CH]O[C@H](OO)C
C4H9O3_TSG = ({0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
               2: ('C', 0, None, None, True), 3: ('C', 0, True, True, None),
               4: ('O', 0, None, None, None), 5: ('O', 0, None, None, None),
               6: ('O', 0, None, None, None), 7: ('H', 0, None, None, None),
               8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
               10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
               12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
               14: ('H', 0, None, None, None), 15: ('H', 0, None, None, None)},
              {frozenset({4, 6}): (1, None, None, None),
               frozenset({2, 13}): (0.9, None, None, None),
               frozenset({3, 15}): (1, None, None, None),
               frozenset({4, 13}): (0.1, None, None, None),
               frozenset({1, 11}): (1, None, None, None),
               frozenset({3, 6}): (1, None, None, None),
               frozenset({0, 2}): (1, None, None, None),
               frozenset({2, 5}): (1, None, None, None),
               frozenset({1, 12}): (1, None, None, None),
               frozenset({2, 14}): (1, None, None, None),
               frozenset({3, 5}): (1, None, None, None),
               frozenset({1, 3}): (1, None, None, None),
               frozenset({0, 7}): (1, None, None, None),
               frozenset({1, 10}): (1, None, None, None),
               frozenset({0, 8}): (1, None, None, None),
               frozenset({0, 9}): (1, None, None, None)})


def test__from_data():
    """ test getters
    """
    tsg = graph.from_data(
        atm_symb_dct=graph.atom_symbols(C4H9O3_TSG),
        bnd_keys=graph.bond_keys(C4H9O3_TSG),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C4H9O3_TSG)),
        atm_ste_par_dct=graph.atom_stereo_parities(C4H9O3_TSG),
        atm_prd_ste_par_dct=graph.ts_atom_product_stereo_parities(C4H9O3_TSG),
        atm_ts_ste_par_dct=graph.ts_atom_fleeting_stereo_parities(C4H9O3_TSG),
        bnd_ord_dct=graph.bond_orders(C4H9O3_TSG),
        bnd_ste_par_dct=graph.bond_stereo_parities(C4H9O3_TSG),
        bnd_prd_ste_par_dct=graph.ts_bond_product_stereo_parities(C4H9O3_TSG),
        bnd_ts_ste_par_dct=graph.ts_bond_fleeting_stereo_parities(C4H9O3_TSG),
    )
    assert tsg == C4H9O3_TSG


def test__getters():
    """ test getters
    """
    orig_tsg = C4H9O3_TSG
    tsg = graph.ts_graph(
        gra=graph.ts.reactants_graph(orig_tsg),
        frm_bnd_keys=graph.ts.forming_bond_keys(orig_tsg),
        brk_bnd_keys=graph.ts.breaking_bond_keys(orig_tsg),
        atm_prd_ste_par_dct=graph.ts.atom_product_stereo_parities(orig_tsg),
        atm_ts_ste_par_dct=graph.ts.atom_fleeting_stereo_parities(orig_tsg),
        bnd_prd_ste_par_dct=graph.ts.bond_product_stereo_parities(orig_tsg),
        bnd_ts_ste_par_dct=graph.ts.bond_fleeting_stereo_parities(orig_tsg),
    )
    assert tsg == orig_tsg


def test__setters():
    """ test graph setters
    """
    atm_symbs = numpy.array(list('CHON'))
    bnd_ords = numpy.arange(1, 4)
    atm_imp_hyd_vlcs = numpy.arange(0, 4)
    pars = numpy.array([None, True, False])

    print("\nTesting setters for a TS graph...")
    orig_gra = C4H9O3_TSG

    atm_keys = graph.atom_keys(orig_gra)
    bnd_keys = graph.bond_keys(orig_gra)
    natms = len(atm_keys)
    nbnds = len(bnd_keys)

    # atom symbols
    orig_atm_symb_dct = graph.atom_symbols(orig_gra)
    atm_symb_dct = dict(
        zip(atm_keys, numpy.random.choice(atm_symbs, size=natms)))
    gra = graph.set_atom_symbols(orig_gra, atm_symb_dct)
    print(atm_symb_dct)
    assert atm_symb_dct == graph.atom_symbols(gra)
    assert orig_gra == graph.set_atom_symbols(gra, orig_atm_symb_dct)

    # bond orders
    orig_bnd_ord_dct = graph.bond_orders(orig_gra)
    bnd_ord_dct = dict(
        zip(bnd_keys, numpy.random.choice(bnd_ords, size=nbnds)))
    gra = graph.set_bond_orders(orig_gra, bnd_ord_dct)
    print(bnd_ord_dct)
    assert bnd_ord_dct == graph.bond_orders(gra)
    assert orig_gra == graph.set_bond_orders(gra, orig_bnd_ord_dct)

    # atom implicit hydrogen valences
    orig_atm_imp_hyd_vlc_dct = graph.atom_implicit_hydrogen_valences(orig_gra)
    atm_imp_hyd_vlc_dct = dict(
        zip(atm_keys, numpy.random.choice(atm_imp_hyd_vlcs, size=natms)))
    gra = graph.set_atom_implicit_hydrogen_valences(
        orig_gra, atm_imp_hyd_vlc_dct)
    print(atm_imp_hyd_vlc_dct)
    assert atm_imp_hyd_vlc_dct == graph.atom_implicit_hydrogen_valences(gra)
    assert orig_gra == graph.set_atom_implicit_hydrogen_valences(
        gra, orig_atm_imp_hyd_vlc_dct)

    # atom stereo parities
    orig_atm_par_dct = graph.atom_stereo_parities(orig_gra)
    atm_par_dct = dict(zip(atm_keys, numpy.random.choice(pars, size=natms)))
    gra = graph.set_atom_stereo_parities(orig_gra, atm_par_dct)
    print(atm_par_dct)
    assert atm_par_dct == graph.atom_stereo_parities(gra)
    assert orig_gra == graph.set_atom_stereo_parities(gra, orig_atm_par_dct)

    # bond stereo parities
    orig_bnd_par_dct = graph.bond_stereo_parities(orig_gra)
    bnd_par_dct = dict(zip(bnd_keys, numpy.random.choice(pars, size=nbnds)))
    gra = graph.set_bond_stereo_parities(orig_gra, bnd_par_dct)
    print(bnd_par_dct)
    assert bnd_par_dct == graph.bond_stereo_parities(gra)
    assert orig_gra == graph.set_bond_stereo_parities(gra, orig_bnd_par_dct)

    # TS atom product stereo parities
    orig_atm_par_dct = graph.ts_atom_product_stereo_parities(orig_gra)
    atm_par_dct = dict(zip(atm_keys, numpy.random.choice(pars, size=natms)))
    gra = graph.ts_set_atom_product_stereo_parities(orig_gra, atm_par_dct)
    print(atm_par_dct)
    assert atm_par_dct == graph.ts_atom_product_stereo_parities(gra)
    assert orig_gra == graph.ts_set_atom_product_stereo_parities(
        gra, orig_atm_par_dct)

    # TS bond product stereo parities
    orig_bnd_par_dct = graph.ts_bond_product_stereo_parities(orig_gra)
    bnd_par_dct = dict(zip(bnd_keys, numpy.random.choice(pars, size=nbnds)))
    gra = graph.ts_set_bond_product_stereo_parities(orig_gra, bnd_par_dct)
    print(bnd_par_dct)
    assert bnd_par_dct == graph.ts_bond_product_stereo_parities(gra)
    assert orig_gra == graph.ts_set_bond_product_stereo_parities(
        gra, orig_bnd_par_dct)

    # TS atom fleeting stereo parities
    orig_atm_par_dct = graph.ts_atom_fleeting_stereo_parities(orig_gra)
    atm_par_dct = dict(zip(atm_keys, numpy.random.choice(pars, size=natms)))
    gra = graph.ts_set_atom_fleeting_stereo_parities(orig_gra, atm_par_dct)
    print(atm_par_dct)
    assert atm_par_dct == graph.ts_atom_fleeting_stereo_parities(gra)
    assert orig_gra == graph.ts_set_atom_fleeting_stereo_parities(
        gra, orig_atm_par_dct)

    # TS bond fleeting stereo parities
    orig_bnd_par_dct = graph.ts_bond_fleeting_stereo_parities(orig_gra)
    bnd_par_dct = dict(zip(bnd_keys, numpy.random.choice(pars, size=nbnds)))
    gra = graph.ts_set_bond_fleeting_stereo_parities(orig_gra, bnd_par_dct)
    print(bnd_par_dct)
    assert bnd_par_dct == graph.ts_bond_fleeting_stereo_parities(gra)
    assert orig_gra == graph.ts_set_bond_fleeting_stereo_parities(
        gra, orig_bnd_par_dct)


def test__string():
    """ test graph.string and graph.from_string
    """
    gra = C4H9O3_TSG
    assert gra == graph.from_string(graph.string(gra))


def test__atoms_neighbor_atom_keys():
    """ test graph.atoms_neighbor_atom_keys
    """
    assert graph.atoms_neighbor_atom_keys(C4H9O3_TSG) == {
        0: frozenset({2, 7, 8, 9}),
        1: frozenset({3, 10, 11, 12}),
        2: frozenset({0, 5, 13, 14}),
        3: frozenset({1, 5, 6, 15}),
        4: frozenset({6, 13}),
        5: frozenset({2, 3}),
        6: frozenset({3, 4}),
        7: frozenset({0}),
        8: frozenset({0}),
        9: frozenset({0}),
        10: frozenset({1}),
        11: frozenset({1}),
        12: frozenset({1}),
        13: frozenset({2, 4}),
        14: frozenset({2}),
        15: frozenset({3}),
    }
    assert graph.atoms_neighbor_atom_keys(C4H9O3_TSG, ts_=False) == {
        0: frozenset({2, 7, 8, 9}),
        1: frozenset({3, 10, 11, 12}),
        2: frozenset({0, 5, 13, 14}),
        3: frozenset({1, 5, 6, 15}),
        4: frozenset({6}),
        5: frozenset({2, 3}),
        6: frozenset({3, 4}),
        7: frozenset({0}),
        8: frozenset({0}),
        9: frozenset({0}),
        10: frozenset({1}),
        11: frozenset({1}),
        12: frozenset({1}),
        13: frozenset({2}),
        14: frozenset({2}),
        15: frozenset({3}),
    }


def test__atoms_bond_keys():
    """ test graph.atoms_neighbor_atom_keys
    """
    assert graph.atoms_bond_keys(C4H9O3_TSG) == {
        0: frozenset({frozenset({0, 7}), frozenset({0, 2}),
                      frozenset({0, 8}), frozenset({0, 9})}),
        1: frozenset({frozenset({1, 11}), frozenset({1, 10}),
                      frozenset({1, 12}), frozenset({1, 3})}),
        2: frozenset({frozenset({2, 5}), frozenset({2, 14}),
                      frozenset({2, 13}), frozenset({0, 2})}),
        3: frozenset({frozenset({3, 6}), frozenset({1, 3}),
                      frozenset({3, 5}), frozenset({3, 15})}),
        4: frozenset({frozenset({4, 6}), frozenset({4, 13})}),
        5: frozenset({frozenset({2, 5}), frozenset({3, 5})}),
        6: frozenset({frozenset({3, 6}), frozenset({4, 6})}),
        7: frozenset({frozenset({0, 7})}),
        8: frozenset({frozenset({0, 8})}),
        9: frozenset({frozenset({0, 9})}),
        10: frozenset({frozenset({1, 10})}),
        11: frozenset({frozenset({1, 11})}),
        12: frozenset({frozenset({1, 12})}),
        13: frozenset({frozenset({2, 13}), frozenset({4, 13})}),
        14: frozenset({frozenset({2, 14})}),
        15: frozenset({frozenset({3, 15})}),
    }

    assert graph.atoms_bond_keys(C4H9O3_TSG, ts_=False) == {
        0: frozenset({frozenset({0, 7}), frozenset({0, 2}),
                      frozenset({0, 8}), frozenset({0, 9})}),
        1: frozenset({frozenset({1, 11}), frozenset({1, 10}),
                      frozenset({1, 12}), frozenset({1, 3})}),
        2: frozenset({frozenset({2, 5}), frozenset({2, 14}),
                      frozenset({2, 13}), frozenset({0, 2})}),
        3: frozenset({frozenset({3, 6}), frozenset({1, 3}),
                      frozenset({3, 5}), frozenset({3, 15})}),
        4: frozenset({frozenset({4, 6})}),
        5: frozenset({frozenset({2, 5}), frozenset({3, 5})}),
        6: frozenset({frozenset({3, 6}), frozenset({4, 6})}),
        7: frozenset({frozenset({0, 7})}),
        8: frozenset({frozenset({0, 8})}),
        9: frozenset({frozenset({0, 9})}),
        10: frozenset({frozenset({1, 10})}),
        11: frozenset({frozenset({1, 11})}),
        12: frozenset({frozenset({1, 12})}),
        13: frozenset({frozenset({2, 13})}),
        14: frozenset({frozenset({2, 14})}),
        15: frozenset({frozenset({3, 15})}),
    }


# # bond properties
def test__bonds_neighbor_atom_keys():
    """ test graph.bonds_neighbor_atom_keys
    """

    assert graph.bonds_neighbor_atom_keys(C4H9O3_TSG) == {
        frozenset({2, 13}): frozenset({0, 4, 5, 14}),
        frozenset({4, 6}): frozenset({3, 13}),
        frozenset({3, 15}): frozenset({1, 5, 6}),
        frozenset({4, 13}): frozenset({2, 6}),
        frozenset({1, 11}): frozenset({3, 10, 12}),
        frozenset({1, 10}): frozenset({3, 11, 12}),
        frozenset({3, 6}): frozenset({1, 4, 5, 15}),
        frozenset({0, 2}): frozenset({5, 7, 8, 9, 13, 14}),
        frozenset({2, 5}): frozenset({0, 3, 13, 14}),
        frozenset({1, 12}): frozenset({3, 10, 11}),
        frozenset({2, 14}): frozenset({0, 5, 13}),
        frozenset({3, 5}): frozenset({1, 2, 6, 15}),
        frozenset({0, 7}): frozenset({2, 8, 9}),
        frozenset({1, 3}): frozenset({5, 6, 10, 11, 12, 15}),
        frozenset({0, 8}): frozenset({2, 7, 9}),
        frozenset({0, 9}): frozenset({2, 7, 8})
    }

    assert graph.bonds_neighbor_atom_keys(C4H9O3_TSG, ts_=False) == {
        frozenset({2, 13}): frozenset({0, 5, 14}),
        frozenset({4, 6}): frozenset({3}),
        frozenset({3, 15}): frozenset({1, 5, 6}),
        frozenset({1, 11}): frozenset({3, 10, 12}),
        frozenset({1, 10}): frozenset({3, 11, 12}),
        frozenset({3, 6}): frozenset({1, 4, 5, 15}),
        frozenset({0, 2}): frozenset({5, 7, 8, 9, 13, 14}),
        frozenset({2, 5}): frozenset({0, 3, 13, 14}),
        frozenset({1, 12}): frozenset({3, 10, 11}),
        frozenset({2, 14}): frozenset({0, 5, 13}),
        frozenset({3, 5}): frozenset({1, 2, 6, 15}),
        frozenset({0, 7}): frozenset({2, 8, 9}),
        frozenset({1, 3}): frozenset({5, 6, 10, 11, 12, 15}),
        frozenset({0, 8}): frozenset({2, 7, 9}),
        frozenset({0, 9}): frozenset({2, 7, 8})
    }


if __name__ == '__main__':
    test__getters()
    test__setters()
    test__string()
    test__atoms_neighbor_atom_keys()
    test__atoms_bond_keys()
    test__bonds_neighbor_atom_keys()