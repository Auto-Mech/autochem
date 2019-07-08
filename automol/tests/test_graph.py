""" test automol.graph
"""
import numpy
import automol
from automol import graph

C8H13O_CGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)})
C8H13O_RGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({1, 4}): (2, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (2, None), frozenset({5, 7}): (1, None)})
C8H13O_SGR = (
    {0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
    {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
     frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)})


C3H3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({2, 0}): (1, None)})
C3H3_RGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (2, None),
      frozenset({2, 0}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (2, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
     {frozenset({0, 1}): (2, None), frozenset({1, 2}): (1, None),
      frozenset({2, 0}): (1, None)}),
)

C2_CGR = ({0: ('C', 0, None), 1: ('C', 0, None)},
          {frozenset({0, 1}): (1, None)})
C2_RGRS = (
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (1, None)}),
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (2, None)}),
    ({0: ('C', 0, None), 1: ('C', 0, None)},
     {frozenset({0, 1}): (3, None)}),
)

CH2FH2H_CGR_IMP = (
    {0: ('F', 0, None), 1: ('C', 2, None), 2: ('H', 1, None),
     3: ('H', 0, None)},
    {frozenset({0, 1}): (1, None)})
CH2FH2H_CGR_EXP = (
    {0: ('F', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
     3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None), frozenset({2, 6}): (1, None)})

C2H2CL2F2_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('F', 0, None),
     3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
     frozenset({1, 5}): (1, None)})
C2H2CL2F2_SGRS = (
    ({0: ('C', 1, False), 1: ('C', 1, False), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, True), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, False), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, True), 2: ('F', 0, None),
      3: ('Cl', 0, None), 4: ('F', 0, None), 5: ('Cl', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None),
      frozenset({1, 5}): (1, None)})
)

C3H3CL2F3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
     3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
     6: ('F', 0, None), 7: ('F', 0, None)},
    {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
     frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
     frozenset({2, 7}): (1, None)})
C3H3CL2F3_SGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, False), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, True), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, False), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, False), 1: ('C', 1, True), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, False), 2: ('C', 1, True),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
    ({0: ('C', 1, True), 1: ('C', 1, True), 2: ('C', 1, False),
      3: ('Cl', 0, None), 4: ('Cl', 0, None), 5: ('F', 0, None),
      6: ('F', 0, None), 7: ('F', 0, None)},
     {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
      frozenset({0, 5}): (1, None), frozenset({2, 4}): (1, None),
      frozenset({1, 3}): (1, None), frozenset({1, 6}): (1, None),
      frozenset({2, 7}): (1, None)}),
)

C3H5N3_CGR = (
    {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
     3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
    {frozenset({1, 4}): (1, None), frozenset({1, 2}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({0, 2}): (1, None),
     frozenset({2, 5}): (1, None)})
C3H5N3_SGRS = (
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, None)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, False)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, False), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, True)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, False)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, False), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, True)}),
    ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 0, None),
      3: ('N', 1, None), 4: ('N', 1, None), 5: ('N', 1, None)},
     {frozenset({1, 4}): (1, True), frozenset({1, 2}): (1, None),
      frozenset({0, 3}): (1, True), frozenset({0, 2}): (1, None),
      frozenset({2, 5}): (1, None)}),
)

C8H13O_SGRS = (
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, False), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, False), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, False), frozenset({5, 7}): (1, None)}),
    ({0: ('C', 3, None), 1: ('C', 2, None), 2: ('C', 3, None),
      3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
      6: ('C', 1, True), 7: ('C', 1, True), 8: ('O', 0, None)},
     {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None),
      frozenset({0, 3}): (1, None), frozenset({2, 6}): (1, None),
      frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None),
      frozenset({3, 5}): (1, True), frozenset({5, 7}): (1, None)}),
)


def test__from_data():
    """ test getters
    """
    cgr = automol.graph.from_data(
        atm_sym_dct=graph.atom_symbols(C8H13O_CGR),
        bnd_keys=graph.bond_keys(C8H13O_CGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_CGR)),
    )
    assert cgr == C8H13O_CGR

    rgr = automol.graph.from_data(
        atm_sym_dct=graph.atom_symbols(C8H13O_RGR),
        bnd_keys=graph.bond_keys(C8H13O_RGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_RGR)),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGR),
    )
    assert rgr == C8H13O_RGR

    sgr = automol.graph.from_data(
        atm_sym_dct=graph.atom_symbols(C8H13O_SGR),
        bnd_keys=graph.bond_keys(C8H13O_SGR),
        atm_imp_hyd_vlc_dct=(
            graph.atom_implicit_hydrogen_valences(C8H13O_SGR)),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGR),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGR)
    )
    assert sgr == C8H13O_SGR


def test__set_atom_implicit_hydrogen_valences():
    """ test graph.set_atom_implicit_hydrogen_valences
    """
    atm_keys = graph.atom_keys(C8H13O_CGR)
    cgr = graph.set_atom_implicit_hydrogen_valences(
        C8H13O_CGR, {atm_key: 0 for atm_key in atm_keys})

    assert cgr == automol.graph.from_data(
        graph.atom_symbols(C8H13O_CGR), graph.bond_keys(C8H13O_CGR))


def test__without_bond_orders():
    """ test graph.without_bond_orders
    """
    assert C8H13O_CGR == graph.without_bond_orders(C8H13O_RGR)


def test__without_stereo_parities():
    """ test graph.without_stereo_parities
    """
    assert C8H13O_CGR == graph.without_stereo_parities(C8H13O_SGR)


# graph theory library
# # atom properties
def test__atom_neighbor_keys():
    """ test graph.atom_neighbor_keys
    """
    assert graph.atom_neighbor_keys(C8H13O_CGR) == {
        0: frozenset({3}),
        1: frozenset({4}),
        2: frozenset({6}),
        3: frozenset({0, 5}),
        4: frozenset({1, 6}),
        5: frozenset({3, 7}),
        6: frozenset({2, 4, 7}),
        7: frozenset({8, 5, 6}),
        8: frozenset({7})
    }


def test__atom_bond_keys():
    """ test graph.atom_neighbor_keys
    """
    assert graph.atom_bond_keys(C8H13O_CGR) == {
        0: frozenset({frozenset({0, 3})}),
        1: frozenset({frozenset({1, 4})}),
        2: frozenset({frozenset({2, 6})}),
        3: frozenset({frozenset({3, 5}), frozenset({0, 3})}),
        4: frozenset({frozenset({1, 4}), frozenset({4, 6})}),
        5: frozenset({frozenset({3, 5}), frozenset({5, 7})}),
        6: frozenset({frozenset({6, 7}), frozenset({4, 6}),
                      frozenset({2, 6})}),
        7: frozenset({frozenset({6, 7}), frozenset({5, 7}),
                      frozenset({8, 7})}),
        8: frozenset({frozenset({8, 7})})
    }


# # bond properties
def test__bond_neighbor_keys():
    """ test graph.bond_neighbor_keys
    """
    assert graph.bond_neighbor_keys(C8H13O_CGR) == {
        frozenset({1, 4}): frozenset({frozenset({4, 6})}),
        frozenset({4, 6}): frozenset({frozenset({6, 7}), frozenset({1, 4}),
                                      frozenset({2, 6})}),
        frozenset({2, 6}): frozenset({frozenset({6, 7}), frozenset({4, 6})}),
        frozenset({0, 3}): frozenset({frozenset({3, 5})}),
        frozenset({6, 7}): frozenset({frozenset({4, 6}), frozenset({8, 7}),
                                      frozenset({5, 7}), frozenset({2, 6})}),
        frozenset({8, 7}): frozenset({frozenset({6, 7}), frozenset({5, 7})}),
        frozenset({3, 5}): frozenset({frozenset({5, 7}), frozenset({0, 3})}),
        frozenset({5, 7}): frozenset({frozenset({6, 7}), frozenset({3, 5}),
                                      frozenset({8, 7})})
    }


# # other properties
def test__branch():
    """ test graph.branch
    """
    assert graph.branch(C8H13O_CGR, 6, frozenset({6, 4})) == (
        {1: ('C', 2, None), 4: ('C', 1, None), 6: ('C', 1, None)},
        {frozenset({1, 4}): (1, None), frozenset({4, 6}): (1, None)}
    )


def test__rings():
    """ test graph.rings
    """
    c5h5n5o_cgr = (
        {0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
         3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
         6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
         9: ('N', 1, None), 10: ('O', 1, None)},
        {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
         frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
         frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
         frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
         frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})

    assert graph.rings(c5h5n5o_cgr) == (
        ({0: ('C', 1, None), 1: ('C', 0, None), 3: ('C', 0, None),
          6: ('N', 0, None), 7: ('N', 0, None)},
         {frozenset({0, 6}): (1, None), frozenset({3, 7}): (1, None),
          frozenset({0, 7}): (1, None), frozenset({1, 6}): (1, None),
          frozenset({1, 3}): (1, None)}),
        ({1: ('C', 0, None), 2: ('C', 0, None), 3: ('C', 0, None),
          4: ('C', 0, None), 8: ('N', 0, None), 9: ('N', 1, None)},
         {frozenset({8, 2}): (1, None), frozenset({9, 3}): (1, None),
          frozenset({1, 2}): (1, None), frozenset({9, 4}): (1, None),
          frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    )


def test__connected_components():
    """ test graph.connected_components
    """
    gra1 = C3H3_CGR
    gra2 = C2_CGR
    gra1_natms = automol.formula.atom_count(graph.formula(C3H3_CGR))
    gra2 = graph.transform_keys(gra2, lambda x: x + gra1_natms)

    gra1 = gra1
    gra2 = gra2
    gra = graph.union(gra1, gra2)
    cmp_gras = graph.connected_components(gra)
    assert cmp_gras in [(gra1, gra2), (gra2, gra1)]


def test__subgraph():
    """ test graph.subgraph
    """
    assert graph.subgraph(C3H3_CGR, (1, 2)) == (
        {1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({1, 2}): (1, None)})


def test__bond_induced_subgraph():
    """ test graph.bond_induced_subgraph
    """
    assert graph.bond_induced_subgraph(
        C3H3_CGR, [frozenset({0, 1}), frozenset({1, 2})]) == (
            {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None)})


# # transformations
def test__relabel():
    """ test graph.relabel
    """
    assert graph.relabel(C3H3_CGR, {0: 10, 1: 11, 2: 12}) == (
        {10: ('C', 1, None), 11: ('C', 1, None), 12: ('C', 1, None)},
        {frozenset({10, 11}): (1, None), frozenset({11, 12}): (1, None),
         frozenset({12, 10}): (1, None)})


def test__remove_atoms():
    """ test graph.remove_atoms
    """
    assert graph.remove_atoms(C3H3_CGR, (0,)) == (
        {1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({1, 2}): (1, None)})


def test__remove_bonds():
    """ test graph.remove_bonds
    """
    assert graph.remove_bonds(C3H3_CGR, [frozenset({1, 2})]) == (
        {0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None)},
        {frozenset({0, 1}): (1, None), frozenset({2, 0}): (1, None)})


# implicit/explicit hydrogen functions
# # atom properties
def test__atom_explicit_hydrogen_valences():
    """ test graph.atom_explicit_hydrogen_valences
    """
    assert graph.atom_explicit_hydrogen_valences(CH2FH2H_CGR_EXP) == {
        0: 0, 1: 2, 2: 1, 3: 0, 4: 0, 5: 0, 6: 0
    }


def test__atom_explicit_hydrogen_keys():
    """ test graph.atom_explicit_hydrogen_keys
    """
    assert graph.atom_explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == {
        0: frozenset(),
        1: frozenset({4, 5}),
        2: frozenset({6}),
        3: frozenset(),
        4: frozenset(),
        5: frozenset(),
        6: frozenset()
    }


# # other properties
def test__backbone_keys():
    """ test graph.backbone_keys
    """
    assert graph.backbone_keys(CH2FH2H_CGR_EXP) == frozenset({0, 1, 2, 3})


def test__explicit_hydrogen_keys():
    """ test graph.explicit_hydrogen_keys
    """
    assert graph.explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == frozenset(
        {4, 5, 6})


def test__explicit():
    """ test graph.explicit
    """
    assert CH2FH2H_CGR_EXP == graph.explicit(CH2FH2H_CGR_IMP)


def test__implicit():
    """ test graph.implicit
    """
    assert CH2FH2H_CGR_IMP == graph.implicit(graph.explicit(CH2FH2H_CGR_IMP))


# # comparisons
def test__backbone_isomorphic():
    """ test graph.backbone_isomorphic
    """
    assert graph.backbone_isomorphic(CH2FH2H_CGR_IMP, CH2FH2H_CGR_EXP)

    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphic(cgr, cgr_pmt)


def test__backbone_isomorphism():
    """ test graph.backbone_isomorphism
    """
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphism(cgr, cgr_pmt) == pmt_dct


def test__backbone_unique():
    """ test graph.backbone_unique
    """
    assert graph.backbone_unique(C3H3_RGRS) == C3H3_RGRS[:2]


# chemistry library
def test__atom_element_valences():
    """ test graph.atom_element_valences
    """
    assert graph.atom_element_valences(C8H13O_CGR) == {
        0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 2}


def test__atom_lone_pair_counts():
    """ test graph.atom_lone_pair_counts
    """
    assert graph.atom_lone_pair_counts(C8H13O_CGR) == {
        0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 2}


def test__atom_bond_valences():
    """ test graph.atom_bond_valences
    """
    assert graph.atom_bond_valences(C8H13O_CGR) == {
        0: 4, 1: 3, 2: 4, 3: 3, 4: 3, 5: 3, 6: 4, 7: 4, 8: 1}


def test__atom_unsaturated_valences():
    """ test graph.atom_unsaturated_valences
    """
    assert graph.atom_unsaturated_valences(C8H13O_CGR) == {
        0: 0, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1, 6: 0, 7: 0, 8: 1}


def test__unsaturated_atom_keys():
    """ test graph.unsaturated_atom_keys
    """
    assert graph.unsaturated_atom_keys(C8H13O_CGR) == frozenset(
        {1, 3, 4, 5, 8})


def test__maximum_spin_multiplicity():
    """ test graph.maximum_spin_multiplicity
    """
    assert graph.maximum_spin_multiplicity(C2_CGR) == 7


def test__possible_spin_multiplicities():
    """ test graph.possible_spin_multiplicities
    """
    assert graph.possible_spin_multiplicities(C2_CGR) == (1, 3, 5, 7)


# miscellaneous
def test__bond_symmetry_numbers():
    """ test graph.bond_symmetry_numbers
    """
    assert graph.bond_symmetry_numbers(C8H13O_CGR) == {
        frozenset({1, 4}): 1, frozenset({4, 6}): 1, frozenset({2, 6}): 3,
        frozenset({0, 3}): 3, frozenset({6, 7}): 1, frozenset({8, 7}): 1,
        frozenset({3, 5}): 1, frozenset({5, 7}): 1}


# resonance graph library
# # atom properties
def test__resonance_dominant_atom_hybridizations():
    """ test graph.resonance_dominant_atom_hybridizations
    """
    assert graph.resonance_dominant_atom_hybridizations(C3H3_CGR) == {
        0: 2, 1: 2, 2: 2}
    assert graph.resonance_dominant_atom_hybridizations(C8H13O_CGR) == {
        0: 3, 1: 2, 2: 3, 3: 2, 4: 2, 5: 2, 6: 3, 7: 3, 8: 3}

    cgr = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('O', 0, None),
            3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
            6: ('X', 0, None)},
           {frozenset({1, 4}): (1, None), frozenset({1, 2}): (1, None),
            frozenset({0, 3}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({2, 5}): (1, None)})
    print(graph.resonance_dominant_atom_hybridizations(cgr))


def test__resonance_dominant_atom_centered_cumulene_keys():
    """ test graph.resonance_dominant_atom_centered_cumulene_keys
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None),
            6: ('C', 0, None)},
           {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
            frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    assert (graph.resonance_dominant_atom_centered_cumulene_keys(cgr) ==
            frozenset({(frozenset({1, 4}), 5)}))


def test__resonance_dominant_bond_centered_cumulene_keys():
    """ test graph.resonance_dominant_bond_centered_cumulene_keys
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 2, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 1, None), 5: ('C', 0, None)},
           {frozenset({4, 5}): (1, None), frozenset({0, 2}): (1, None),
            frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({1, 3}): (1, None)})
    assert (graph.resonance_dominant_bond_centered_cumulene_keys(cgr) ==
            frozenset({(frozenset({1, 4}), frozenset({3, 5}))}))


def test__resonance_dominant_radical_atom_keys():
    """ test graph.resonance_dominant_radical_atom_keys
    """
    assert graph.resonance_dominant_radical_atom_keys(C3H3_CGR) == frozenset(
        {0, 1, 2})
    assert graph.resonance_dominant_radical_atom_keys(C8H13O_CGR) == frozenset(
        {8})


# # bond properties
def test__resonance_dominant_bond_orders():
    """ test graph.resonance_dominant_bond_orders
    """
    assert graph.resonance_dominant_bond_orders(C3H3_CGR) == {
        frozenset({0, 1}): frozenset({1, 2}),
        frozenset({0, 2}): frozenset({1, 2}),
        frozenset({1, 2}): frozenset({1, 2})
    }


# # transformations
def test__resonances():
    """ test graph.resonances
    """
    assert graph.resonances(C3H3_CGR) == C3H3_RGRS


def test__subresonances():
    """ test graph.subresonances
    """
    assert graph.subresonances(C2_RGRS[1]) == C2_RGRS[1:]


def test__dominant_resonances():
    """ test graph.dominant_resonances
    """
    assert graph.dominant_resonances(C3H3_CGR) == C3H3_RGRS[1:]


def test__dominant_resonance():
    """ test graph.dominant_resonance
    """
    assert graph.dominant_resonance(C3H3_CGR) == C3H3_RGRS[1]


# stereo graph library
def test__stereogenic_atom_keys():
    """ test graph.stereogenic_atom_keys
    """
    assert graph.stereogenic_atom_keys(C8H13O_CGR) == frozenset({6, 7})
    assert graph.stereogenic_atom_keys(C3H3CL2F3_CGR) == frozenset({1, 2})

    cgr = ({0: ('C', 2, None), 1: ('C', 3, None), 2: ('C', 1, None),
            3: ('O', 1, None)},
           {frozenset({0, 2}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({1, 2}): (1, None)})
    assert graph.stereogenic_atom_keys(cgr) == frozenset({2})


def test__stereogenic_bond_keys():
    """ test graph.stereogenic_bond_keys
    """
    print(graph.stereogenic_bond_keys(C8H13O_CGR))
    print(graph.stereogenic_bond_keys(C3H5N3_CGR))
    assert graph.stereogenic_bond_keys(C8H13O_CGR) == frozenset(
        {frozenset({3, 5})})
    assert graph.stereogenic_bond_keys(C3H5N3_CGR) == frozenset(
        {frozenset({1, 4}), frozenset({0, 3})})


def test__stereomers():
    """ test graph.stereomers
    """
    assert graph.stereomers(C2H2CL2F2_CGR) == C2H2CL2F2_SGRS
    assert graph.stereomers(C3H3CL2F3_CGR) == C3H3CL2F3_SGRS
    assert graph.stereomers(C3H5N3_CGR) == C3H5N3_SGRS
    assert graph.stereomers(C8H13O_CGR) == C8H13O_SGRS


def test__atom_stereo_coordinates():
    """ test graph.atom_stereo_coordinates
    """
    for sgr in C2H2CL2F2_SGRS:
        sgr = graph.explicit(sgr)
        cgr = graph.without_stereo_parities(sgr)
        atm_xyz_dct = graph.atom_stereo_coordinates(sgr)
        assert graph.set_stereo_from_atom_coordinates(cgr, atm_xyz_dct) == sgr

    for sgr in C3H5N3_SGRS:
        sgr = graph.explicit(sgr)
        cgr = graph.without_stereo_parities(sgr)
        atm_xyz_dct = graph.atom_stereo_coordinates(sgr)
        assert graph.set_stereo_from_atom_coordinates(cgr, atm_xyz_dct) == sgr

    for sgr in C3H3CL2F3_SGRS:
        sgr = graph.explicit(sgr)
        cgr = graph.without_stereo_parities(sgr)
        atm_xyz_dct = graph.atom_stereo_coordinates(sgr)
        assert graph.set_stereo_from_atom_coordinates(cgr, atm_xyz_dct) == sgr

    for sgr in C8H13O_SGRS:
        sgr = graph.explicit(sgr)
        cgr = graph.without_stereo_parities(sgr)
        atm_xyz_dct = graph.atom_stereo_coordinates(sgr)
        assert graph.set_stereo_from_atom_coordinates(cgr, atm_xyz_dct) == sgr


def test__trans__hydrogen_migration():
    """ test graph.trans.hydrogen_migration
    """
    cgr1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('C', 1, None), 5: ('O', 0, None)},
            {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
             frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
             frozenset({0, 1}): (1, None)})
    cgr2 = ({0: ('C', 2, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('C', 0, None), 5: ('O', 0, None)},
            {frozenset({3, 4}): (1, None), frozenset({2, 3}): (1, None),
             frozenset({1, 2}): (1, None), frozenset({4, 5}): (1, None),
             frozenset({0, 1}): (1, None)})

    cgr1 = graph.explicit(cgr1)
    cgr2 = graph.explicit(cgr2)

    tra = graph.trans.hydrogen_migration(cgr1, cgr2)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)

    tra = graph.trans.hydrogen_migration(cgr2, cgr1)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr2), cgr1)


def test__trans__beta_scission():
    """ test graph.trans.beta_scission
    """
    cgr1 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({3, 6}): (1, None), frozenset({2, 4}): (1, None),
             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    cgr2 = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 1, None),
             3: ('C', 1, None), 4: ('F', 0, None), 5: ('F', 0, None),
             6: ('O', 1, None)},
            {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
             frozenset({1, 3}): (1, None)})

    cgr1 = graph.explicit(cgr1)
    cgr2 = graph.explicit(cgr2)

    tra = graph.trans.beta_scission(cgr1, cgr2)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)


def test__trans__addition():
    """ test graph.trans.addition
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

    tra = graph.trans.addition(cgr1, cgr2)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)


def test__trans__hydrogen_abstraction():
    """ test graph.trans.hydrogen_abstraction
    """
    cgr1 = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
             3: ('C', 2, None), 4: ('C', 1, None), 5: ('C', 2, None),
             6: ('C', 2, None), 7: ('O', 1, None)},
            {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})
    cgr2 = ({0: ('C', 2, None), 1: ('C', 3, None), 2: ('C', 1, None),
             3: ('C', 2, None), 4: ('C', 1, None), 5: ('C', 2, None),
             6: ('C', 2, None), 7: ('O', 2, None)},
            {frozenset({4, 6}): (1, None), frozenset({0, 2}): (1, None),
             frozenset({2, 4}): (1, None), frozenset({5, 6}): (1, None),
             frozenset({3, 5}): (1, None), frozenset({1, 3}): (1, None)})

    cgr1 = graph.explicit(cgr1)
    cgr2 = graph.explicit(cgr2)

    tra = graph.trans.hydrogen_abstraction(cgr1, cgr2)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)

    tra = graph.trans.hydrogen_abstraction(cgr2, cgr1)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr2), cgr1)


def test__trans__is_stereo_compatible():
    """ test graph.trans.is_stereo_compatible
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

    ich1 = graph.inchi(cgr1)
    ich2 = graph.inchi(cgr2)
    smi1 = automol.inchi.smiles(ich1)
    smi2 = automol.inchi.smiles(ich2)
    print(smi1)
    print(smi2)

    tra = graph.trans.addition(cgr1, cgr2)
    assert graph.backbone_isomorphic(graph.trans.apply(tra, cgr1), cgr2)

    sgr1 = graph.stereomers(cgr1)[0]
    print(sgr1)
    for sgr2 in graph.stereomers(cgr2):
        print(sgr2)
        print(graph.trans.is_stereo_compatible(tra, sgr1, sgr2))


if __name__ == '__main__':
    # test__from_data()
    # test__set_atom_implicit_hydrogen_valences()
    # test__without_bond_orders()
    # test__without_stereo_parities()
    # test__atom_explicit_hydrogen_valences()
    # test__atom_explicit_hydrogen_keys()
    # test__explicit()
    # test__backbone_keys()
    # test__explicit_hydrogen_keys()
    # test__stereomers()
    # test__atom_stereo_coordinates()
    # test__connected_components()
    # test__unsaturated_atom_keys()
    # test__resonance_dominant_radical_atom_keys()
    # test__remove_bonds()
    # test__trans__hydrogen_migration()
    # test__trans__beta_scission()
    # test__trans__addition()
    # test__trans__hydrogen_abstraction()
    # test__trans__is_stereo_compatible()
    # test__resonance_dominant_atom_centered_cumulene_keys()
    # test__resonance_dominant_bond_centered_cumulene_keys()
    # test__stereogenic_bond_keys()
    # test__trans__beta_scission()
    test__resonance_dominant_atom_hybridizations()
