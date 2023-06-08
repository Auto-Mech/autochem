""" test automol.graph.ts
"""
import numpy
from automol import graph


# CCO[C@H](O[O])C => C[CH]O[C@H](OO)C
C4H9O3_TSG = (
    {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
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

# CCOCC + [OH] => C[CH]OCC + O
#  *
# [* marks a fleeting TS stereosite]
C4H11O2_TSG = ({0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
                2: ('C', 0, None, None, True), 3: ('C', 0, None, None, None),
                4: ('O', 0, None, None, None), 5: ('H', 0, None, None, None),
                6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
                8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
                10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
                12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
                14: ('H', 0, None, None, None), 15: ('O', 0, None, None, None),
                16: ('H', 0, None, None, None)},
               {frozenset({0, 5}): (1, None, None, None),
                frozenset({3, 14}): (1, None, None, None),
                frozenset({11, 15}): (0.1, None, None, None),
                frozenset({2, 11}): (0.9, None, None, None),
                frozenset({1, 10}): (1, None, None, None),
                frozenset({3, 4}): (1, None, None, None),
                frozenset({1, 9}): (1, None, None, None),
                frozenset({0, 6}): (1, None, None, None),
                frozenset({0, 2}): (1, None, None, None),
                frozenset({3, 13}): (1, None, None, None),
                frozenset({2, 12}): (1, None, None, None),
                frozenset({15, 16}): (1, None, None, None),
                frozenset({2, 4}): (1, None, None, None),
                frozenset({1, 8}): (1, None, None, None),
                frozenset({0, 7}): (1, None, None, None),
                frozenset({1, 3}): (1, None, None, None)})

# C=C(O[O])OO => [CH]=C(OO)OO
#  *
# [* marks a fleeting TS stereosite]
C2H3O4_TSG = ({0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
               2: ('O', 0, None, None, None), 3: ('O', 0, None, None, None),
               4: ('O', 0, None, None, None), 5: ('O', 0, None, None, None),
               6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
               8: ('H', 0, None, None, None)},
              {frozenset({2, 8}): (1, None, None, None),
               frozenset({1, 4}): (1, None, None, None),
               frozenset({0, 6}): (0.9, None, None, None),
               frozenset({0, 1}): (1, None, None, True),
               frozenset({3, 6}): (0.1, None, None, None),
               frozenset({2, 4}): (1, None, None, None),
               frozenset({1, 5}): (1, None, None, None),
               frozenset({3, 5}): (1, None, None, None),
               frozenset({0, 7}): (1, None, None, None)})

# F/C=C([C@@H](F)O)\[C@H](F)O + [OH] => F[C@H]([C]([C@@H](F)O)[C@H](F)O)O
C4H5F3O2_TSG = (
    {0: ('C', 0, None, False, None), 1: ('C', 0, None, None, None),
     2: ('C', 0, False, False, None), 3: ('C', 0, True, True, None),
     4: ('F', 0, None, None, None), 5: ('F', 0, None, None, None),
     6: ('F', 0, None, None, None), 7: ('O', 0, None, None, None),
     8: ('O', 0, None, None, None), 9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
     14: ('O', 0, None, None, None), 15: ('H', 0, None, None, None)},
    {frozenset({7, 12}): (1, None, None, None),
     frozenset({2, 10}): (1, None, None, None),
     frozenset({1, 2}): (1, None, None, None),
     frozenset({0, 1}): (1, False, None, None),
     frozenset({3, 6}): (1, None, None, None),
     frozenset({2, 5}): (1, None, None, None),
     frozenset({0, 4}): (1, None, None, None),
     frozenset({3, 8}): (1, None, None, None),
     frozenset({0, 14}): (0.1, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({8, 13}): (1, None, None, None),
     frozenset({14, 15}): (1, None, None, None),
     frozenset({3, 11}): (1, None, None, None),
     frozenset({2, 7}): (1, None, None, None),
     frozenset({0, 9}): (1, None, None, None)})

# F[C@H]([C]([C@@H](F)O)[C@H](F)O)O => F/C=C([C@@H](F)O)\[C@H](F)O + [OH]
C4H5F3O2_REV_TSG = (
    {0: ('C', 0, False, None, None),
     1: ('C', 0, None, None, None),
     2: ('C', 0, False, False, None),
     3: ('C', 0, True, True, None),
     4: ('F', 0, None, None, None),
     5: ('F', 0, None, None, None),
     6: ('F', 0, None, None, None),
     7: ('O', 0, None, None, None),
     8: ('O', 0, None, None, None),
     9: ('H', 0, None, None, None),
     10: ('H', 0, None, None, None),
     11: ('H', 0, None, None, None),
     12: ('H', 0, None, None, None),
     13: ('H', 0, None, None, None),
     14: ('O', 0, None, None, None),
     15: ('H', 0, None, None, None)},
    {frozenset({7, 12}): (1, None, None, None),
     frozenset({2, 10}): (1, None, None, None),
     frozenset({1, 2}): (1, None, None, None),
     frozenset({0, 1}): (1, None, False, None),
     frozenset({3, 6}): (1, None, None, None),
     frozenset({2, 5}): (1, None, None, None),
     frozenset({0, 4}): (1, None, None, None),
     frozenset({3, 8}): (1, None, None, None),
     frozenset({0, 14}): (0.9, None, None, None),
     frozenset({1, 3}): (1, None, None, None),
     frozenset({8, 13}): (1, None, None, None),
     frozenset({14, 15}): (1, None, None, None),
     frozenset({3, 11}): (1, None, None, None),
     frozenset({2, 7}): (1, None, None, None),
     frozenset({0, 9}): (1, None, None, None)})

# # CC(OO)[CH]C => CC(O[O])CC
# #        *
# # [* marks a fleeting TS stereosite]
# C4H9O2_TSG = (
#     {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
#      2: ('C', 0, None, None, None), 3: ('C', 0, None, None, None),
#      4: ('O', 0, None, None, None), 5: ('O', 0, None, None, None),
#      6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
#      8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
#      10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
#      12: ('H', 0, None, None, None), 13: ('H', 0, None, None, None),
#      14: ('H', 0, None, None, None)},
#     {frozenset({1, 9}): (1, None, None, None),
#      frozenset({0, 6}): (1, None, None, None),
#      frozenset({2, 3}): (1, None, None, None),
#      frozenset({1, 11}): (1, None, None, None),
#      frozenset({4, 5}): (1, None, None, None),
#      frozenset({0, 2}): (1, None, None, None),
#      frozenset({3, 13}): (1, None, None, None),
#      frozenset({2, 12}): (1, None, None, None),
#      frozenset({2, 14}): (0.1, None, None, None),
#      frozenset({3, 5}): (1, None, None, None),
#      frozenset({4, 14}): (0.9, None, None, None),
#      frozenset({1, 3}): (1, None, None, None),
#      frozenset({0, 7}): (1, None, None, None),
#      frozenset({1, 10}): (1, None, None, None),
#      frozenset({0, 8}): (1, None, None, None)})

# # CC(O)C + [OH] => CC(O)O + [CH3]
# #  *
# # [* marks a fleeting TS stereosite]
# C3H9O2_TSG = (
#     {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
#      2: ('C', 0, None, None, True), 3: ('O', 0, None, None, None),
#      4: ('H', 0, None, None, None), 5: ('H', 0, None, None, None),
#      6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
#      8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
#      10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
#      12: ('O', 0, None, None, None), 13: ('H', 0, None, None, None)},
#     {frozenset({1, 9}): (1, None, None, None),
#      frozenset({1, 7}): (1, None, None, None),
#      frozenset({0, 6}): (1, None, None, None),
#      frozenset({2, 3}): (1, None, None, None),
#      frozenset({1, 2}): (0.9, None, None, None),
#      frozenset({2, 10}): (1, None, None, None),
#      frozenset({0, 2}): (1, None, None, None),
#      frozenset({0, 4}): (1, None, None, None),
#      frozenset({0, 5}): (1, None, None, None),
#      frozenset({2, 12}): (0.1, None, None, None),
#      frozenset({1, 8}): (1, None, None, None),
#      frozenset({12, 13}): (1, None, None, None),
#      frozenset({3, 11}): (1, None, None, None)})

# # CC(O)C + [CH3] => CC(O)C + [CH3]
# #  *
# # [* marks a fleeting TS stereosite]
# C4H11O_TSG = (
#     {0: ('C', 0, None, None, None), 1: ('C', 0, None, None, None),
#      2: ('C', 0, None, None, True), 3: ('O', 0, None, None, None),
#      4: ('H', 0, None, None, None), 5: ('H', 0, None, None, None),
#      6: ('H', 0, None, None, None), 7: ('H', 0, None, None, None),
#      8: ('H', 0, None, None, None), 9: ('H', 0, None, None, None),
#      10: ('H', 0, None, None, None), 11: ('H', 0, None, None, None),
#      12: ('C', 0, None, None, None), 13: ('H', 0, None, None, None),
#      14: ('H', 0, None, None, None), 15: ('H', 0, None, None, None)},
#     {frozenset({1, 9}): (1, None, None, None),
#      frozenset({1, 7}): (1, None, None, None),
#      frozenset({0, 6}): (1, None, None, None),
#      frozenset({2, 3}): (1, None, None, None),
#      frozenset({1, 2}): (0.9, None, None, None),
#      frozenset({2, 10}): (1, None, None, None),
#      frozenset({12, 14}): (1, None, None, None),
#      frozenset({0, 2}): (1, None, None, None),
#      frozenset({12, 15}): (1, None, None, None),
#      frozenset({0, 4}): (1, None, None, None),
#      frozenset({0, 5}): (1, None, None, None),
#      frozenset({2, 12}): (0.1, None, None, None),
#      frozenset({1, 8}): (1, None, None, None),
#      frozenset({12, 13}): (1, None, None, None),
#      frozenset({3, 11}): (1, None, None, None)})


def test__from_data():
    """ test getters
    """
    tsg = graph.from_data(
        atm_symb_dct=graph.atom_symbols(C4H9O3_TSG),
        bnd_keys=graph.bond_keys(C4H9O3_TSG),
        atm_imp_hyd_dct=graph.atom_implicit_hydrogens(C4H9O3_TSG),
        atm_ste_par_dct=graph.atom_stereo_parities(C4H9O3_TSG),
        atm_prd_ste_par_dct=graph.atom_stereo_parities(
            C4H9O3_TSG, ts_select='P'),
        atm_ts_ste_par_dct=graph.atom_stereo_parities(
            C4H9O3_TSG, ts_select='T'),
        bnd_ord_dct=graph.bond_orders(C4H9O3_TSG),
        bnd_ste_par_dct=graph.bond_stereo_parities(C4H9O3_TSG),
        bnd_prd_ste_par_dct=graph.bond_stereo_parities(
            C4H9O3_TSG, ts_select='P'),
        bnd_ts_ste_par_dct=graph.bond_stereo_parities(
            C4H9O3_TSG, ts_select='T'),
    )
    assert tsg == C4H9O3_TSG


def test__getters():
    """ test getters
    """
    ref_tsg = C4H9O3_TSG
    tsg = graph.ts_graph(
        gra=graph.ts.reactants_graph(ref_tsg),
        frm_bnd_keys=graph.ts.forming_bond_keys(ref_tsg),
        brk_bnd_keys=graph.ts.breaking_bond_keys(ref_tsg),
        atm_prd_ste_par_dct=graph.atom_stereo_parities(ref_tsg, ts_select='P'),
        atm_ts_ste_par_dct=graph.atom_stereo_parities(ref_tsg, ts_select='T'),
        bnd_prd_ste_par_dct=graph.bond_stereo_parities(ref_tsg, ts_select='P'),
        bnd_ts_ste_par_dct=graph.bond_stereo_parities(ref_tsg, ts_select='T'),
    )
    assert tsg == ref_tsg


def test__setters():
    """ test graph setters
    """
    atm_symbs = numpy.array(list('CHON'))
    bnd_ords = numpy.arange(1, 4)
    atm_imp_hyds = numpy.arange(0, 4)
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
    orig_atm_imp_hyd_dct = graph.atom_implicit_hydrogens(orig_gra)
    atm_imp_hyd_dct = dict(
        zip(atm_keys, numpy.random.choice(atm_imp_hyds, size=natms)))
    gra = graph.set_atom_implicit_hydrogens(orig_gra, atm_imp_hyd_dct)
    print(atm_imp_hyd_dct)
    assert atm_imp_hyd_dct == graph.atom_implicit_hydrogens(gra)
    assert orig_gra == graph.set_atom_implicit_hydrogens(
        gra, orig_atm_imp_hyd_dct)

    for ts_select in ('R', 'P', 'T'):
        # atom stereo parities
        orig_atm_par_dct = graph.atom_stereo_parities(
            orig_gra, ts_select=ts_select)
        atm_par_dct = dict(
            zip(atm_keys, numpy.random.choice(pars, size=natms)))
        gra = graph.set_atom_stereo_parities(
            orig_gra, atm_par_dct, ts_select=ts_select)
        print(atm_par_dct)
        assert (atm_par_dct ==
                graph.atom_stereo_parities(gra, ts_select=ts_select))
        assert orig_gra == graph.set_atom_stereo_parities(
            gra, orig_atm_par_dct, ts_select=ts_select)

        # bond stereo parities
        orig_bnd_par_dct = graph.bond_stereo_parities(
            orig_gra, ts_select=ts_select)
        bnd_par_dct = dict(
            zip(bnd_keys, numpy.random.choice(pars, size=nbnds)))
        gra = graph.set_bond_stereo_parities(
            orig_gra, bnd_par_dct, ts_select=ts_select)
        print(bnd_par_dct)
        assert bnd_par_dct == graph.bond_stereo_parities(
            gra, ts_select=ts_select)
        assert orig_gra == graph.set_bond_stereo_parities(
            gra, orig_bnd_par_dct, ts_select=ts_select)


def test__reverse():
    """ test graph.ts.reverse
    """
    assert graph.ts.reverse(C4H5F3O2_TSG) == C4H5F3O2_REV_TSG
    assert graph.ts.reverse(C4H5F3O2_REV_TSG) == C4H5F3O2_TSG


def test__string():
    """ test graph.string and graph.from_string
    """
    gra = C4H9O3_TSG
    assert gra == graph.from_string(graph.string(gra))


def test__has_stereo():
    """ test graph.has_stereo functions
    """
    # CCOCC + [OH] => C[CH]OCC + O
    #  *
    # [* marks a fleeting TS stereosite]
    assert not graph.has_stereo(C4H11O2_TSG, ts_selects='R')
    assert not graph.has_stereo(C4H11O2_TSG, ts_selects='RP')
    assert graph.has_stereo(C4H11O2_TSG, ts_selects='RPT')
    assert graph.has_stereo(C4H11O2_TSG, ts_selects='T')
    assert not graph.has_atom_stereo(C4H11O2_TSG, ts_selects='R')
    assert not graph.has_atom_stereo(C4H11O2_TSG, ts_selects='RP')
    assert graph.has_atom_stereo(C4H11O2_TSG, ts_selects='T')
    assert graph.has_atom_stereo(C4H11O2_TSG, ts_selects='RPT')
    assert not graph.has_bond_stereo(C4H11O2_TSG, ts_selects='R')
    assert not graph.has_bond_stereo(C4H11O2_TSG, ts_selects='RPT')
    # C=C(O[O])OO => [CH]=C(OO)OO
    #  *
    # [* marks a fleeting TS stereosite]
    assert not graph.has_stereo(C2H3O4_TSG, ts_selects='R')
    assert not graph.has_stereo(C2H3O4_TSG, ts_selects='RP')
    assert graph.has_stereo(C2H3O4_TSG, ts_selects='RPT')
    assert graph.has_stereo(C2H3O4_TSG, ts_selects='T')
    assert not graph.has_atom_stereo(C2H3O4_TSG, ts_selects='R')
    assert not graph.has_atom_stereo(C2H3O4_TSG, ts_selects='RPT')
    assert not graph.has_bond_stereo(C2H3O4_TSG, ts_selects='R')
    assert not graph.has_bond_stereo(C2H3O4_TSG, ts_selects='RP')
    assert graph.has_bond_stereo(C2H3O4_TSG, ts_selects='T')
    assert graph.has_bond_stereo(C2H3O4_TSG, ts_selects='RPT')


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
    test__reverse()
    test__string()
    test__atoms_neighbor_atom_keys()
    test__atoms_bond_keys()
    test__bonds_neighbor_atom_keys()
    test__has_stereo()
