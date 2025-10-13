""" test rotors
"""

from automol import chi as chi_
from automol import geom, reac, smiles
from automol.data import rotor, tors

# Species Z-Matrix
C3H7OH_ZMA = geom.zmatrix(smiles.geometry("CCCO"))
C4H5OH_ZMA = geom.zmatrix(smiles.geometry("CCC#CO"))
C6H13OH_ZMA = geom.zmatrix(smiles.geometry("CCC(C)CCO"))
C7H16O2_ZMA = geom.zmatrix(
    chi_.geometry("InChI=1S/C7H16O2/c1-3-5-7(9,4-2)6-8/h8-9H,3-6H2,1-2H3/t7-/m1/s1")
)
H2_ZMA = geom.zmatrix(
    smiles.geometry("[H][H]")
)  # make sure things don't break for a rotor-less structure

C4H7O3_TOR_STR = """
D5:
  axis1: 2
  group1: 1-3-4
  axis2: 5
  group2: 6-7-8-9-10-11-12-13-14
  symmetry: 1
D6:
  axis1: 5
  group1: 1-2-3-4
  axis2: 6
  group2: 7-8-9-10-11-12-13-14
  symmetry: 1
D9:
  axis1: 6
  group1: 1-2-3-4-5-8-9-13-14
  axis2: 7
  group2: 10-11-12
  symmetry: 3
D12:
  axis1: 6
  group1: 1-2-3-4-5-7-9-10-11-12
  axis2: 8
  group2: 13-14
  symmetry: 1
D13:
  axis1: 13
  group1: 14
  axis2: 8
  group2: 1-2-3-4-5-6-7-9-10-11-12
  symmetry: 1
"""
C4H7O3_ZMA = (
    ("C", (None, None, None), (None, None, None), (None, None, None)),
    ("C", (0, None, None), ("R1", None, None), (2.480088, None, None)),
    ("H", (0, 1, None), ("R2", "A2", None), (2.039487, 2.103694, None)),
    ("H", (0, 1, 2), ("R3", "A3", "D3"), (2.052919, 2.108722, 3.147005)),
    ("O", (1, 0, 2), ("R4", "A4", "D4"), (2.497653, 2.286147, 3.278218)),
    ("C", (4, 1, 0), ("R5", "A5", "D5"), (2.71753, 2.051854, 3.874355)),
    ("C", (5, 4, 1), ("R6", "A6", "D6"), (2.849875, 1.968001, 0.904697)),
    ("O", (5, 4, 6), ("R7", "A7", "D7"), (2.625782, 1.934588, 4.187407)),
    ("H", (5, 4, 6), ("R8", "A8", "D8"), (2.056744, 1.797242, 2.134375)),
    ("H", (6, 5, 4), ("R9", "A9", "D9"), (2.056056, 1.909271, 0.924466)),
    ("H", (6, 5, 9), ("R10", "A10", "D10"), (2.057218, 1.9289, 4.190752)),
    ("H", (6, 5, 9), ("R11", "A11", "D11"), (2.055208, 1.91078, 2.090482)),
    ("O", (7, 5, 4), ("R12", "A12", "D12"), (2.683853, 1.893624, 5.081201)),
    ("H", (12, 7, 5), ("R13", "A13", "D13"), (1.821823, 1.777218, 1.560426)),
)


# Transition state ZRXN object
C2H5OH_CH3_ZMA = (
    ("O", (None, None, None), (None, None, None), (None, None, None)),
    ("C", (0, None, None), ("R1", None, None), (2.684898, None, None)),
    ("H", (0, 1, None), ("R2", "A2", None), (1.839547, 1.8705, None)),
    ("C", (1, 0, 2), ("R3", "A3", "D3"), (2.852996, 1.92687, 5.212206)),
    ("H", (1, 0, 3), ("R4", "A4", "D4"), (2.064439, 1.895091, 2.135135)),
    ("H", (1, 0, 3), ("R5", "A5", "D5"), (2.069338, 1.88882, 4.191772)),
    ("H", (3, 1, 0), ("R6", "A6", "D6"), (2.072917, 1.961565, 0.946677)),
    ("H", (3, 1, 6), ("R7", "A7", "D7"), (2.083232, 1.933421, 4.303179)),
    ("H", (3, 1, 6), ("R8", "A8", "D8"), (2.059843, 1.946002, 2.14978)),
    ("X", (6, 3, 1), ("R9", "A9", "D9"), (1.889726, 1.570796, 0.0)),
    ("C", (6, 9, 3), ("R10", "A10", "D10"), (3.105266, 1.626053, 3.400042)),
    ("H", (10, 6, 9), ("R11", "A11", "D11"), (2.034492, 1.382484, 4.516064)),
    ("H", (10, 6, 11), ("R12", "A12", "D12"), (2.04136, 1.155219, 2.567991)),
    ("H", (10, 6, 11), ("R13", "A13", "D13"), (2.11299, 2.777445, 4.005115)),
)


C2H5OH_CH3_ZRXN_STR = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: X, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 1, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  4-7: {order: 0.9, stereo_parity: null}
  4-8: {order: 1, stereo_parity: null}
  4-9: {order: 1, stereo_parity: null}
  7-10: {order: 0, stereo_parity: null}
  7-11: {order: 0.1, stereo_parity: null}
  11-12: {order: 1, stereo_parity: null}
  11-13: {order: 1, stereo_parity: null}
  11-14: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
- [11, 12, 13, 14]
backward TS atoms:
  1: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogens: 0, stereo_parity: null}
  8: {symbol: O, implicit_hydrogens: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogens: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 0.1, stereo_parity: null}
  6-7: {order: 1, stereo_parity: null}
  6-9: {order: 1, stereo_parity: null}
  6-10: {order: 1, stereo_parity: null}
  7-8: {order: 1, stereo_parity: null}
  7-11: {order: 1, stereo_parity: null}
  7-12: {order: 1, stereo_parity: null}
  8-13: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3, 4, 5]
- [6, 7, 8, 9, 10, 11, 12, 13]
"""
C2H5OH_CH3_ZRXN = reac.from_old_string(C2H5OH_CH3_ZRXN_STR)


def test__rotor():
    """test a simple rotor without dummy atoms"""

    rotors = rotor.rotors_from_zmatrix(C3H7OH_ZMA)

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom", sort=True)
    assert len(tors_lst) == 3

    rotor_names1 = rotor.rotors_torsion_names(rotors)
    rotor_names2 = rotor.rotors_torsion_names(rotors, flat=True)
    assert rotor_names1 == (("D5",), ("D8",), ("D11",))
    assert rotor_names2 == ("D5", "D8", "D11")

    rotor_axes1 = rotor.rotors_torsion_axes(rotors)
    rotor_axes2 = rotor.rotors_torsion_axes(rotors, flat=True)
    rotor_axes3 = rotor.rotors_torsion_axes(rotors, "geom")
    assert rotor_axes1 == (((0, 1),), ((1, 5),), ((5, 8),))
    assert rotor_axes2 == ((0, 1), (1, 5), (5, 8))
    assert rotor_axes3 == (((0, 1),), ((1, 5),), ((5, 8),))

    rotor_groups1 = rotor.rotors_torsion_groups(rotors)
    rotor_groups2 = rotor.rotors_torsion_groups(rotors, flat=True)
    rotor_groups3 = rotor.rotors_torsion_groups(rotors, "geom")
    assert rotor_groups1 == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),),
        (((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),),
        (((0, 1, 2, 3, 4, 6, 7, 9, 10), (11,)),),
    )
    assert rotor_groups2 == (
        ((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),
        ((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),
        ((0, 1, 2, 3, 4, 6, 7, 9, 10), (11,)),
    )
    assert rotor_groups3 == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11)),),
        (((0, 2, 3, 4, 6, 7), (8, 9, 10, 11)),),
        (((0, 1, 2, 3, 4, 6, 7, 9, 10), (11,)),),
    )

    rotor_symms1 = rotor.rotors_torsion_symmetries(rotors)
    rotor_symms2 = rotor.rotors_torsion_symmetries(rotors, flat=True)
    assert rotor_symms1 == ((3,), (1,), (1,))
    assert rotor_symms2 == (3, 1, 1)

    rotor_grids1 = rotor.rotors_torsion_grids(rotors)
    rotor_grids2 = rotor.rotors_torsion_grids(rotors, flat=True)
    assert [list(map(len, r)) for r in rotor_grids1] == [[4], [12], [12]]
    assert list(map(len, rotor_grids2)) == [4, 12, 12]


def test__rotor_empty():
    """test a rotor with dummy atoms"""

    rotors = rotor.rotors_from_zmatrix(H2_ZMA)

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom", sort=True)
    assert len(tors_lst) == 0

    rotor_names = rotor.rotors_torsion_names(rotors)
    assert rotor_names == ()

    rotor_axes1 = rotor.rotors_torsion_axes(rotors)
    rotor_axes2 = rotor.rotors_torsion_axes(rotors, "geom")
    assert rotor_axes1 == ()
    assert rotor_axes2 == ()

    rotor_groups1 = rotor.rotors_torsion_groups(rotors)
    rotor_groups2 = rotor.rotors_torsion_groups(rotors, "geom")
    assert rotor_groups1 == ()
    assert rotor_groups2 == ()

    rotor_symms = rotor.rotors_torsion_symmetries(rotors)
    assert rotor_symms == ()


def test__rotor_with_dummy_atoms():
    """test a rotor with dummy atoms"""

    rotors = rotor.rotors_from_zmatrix(C4H5OH_ZMA)

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom", sort=True)
    assert len(tors_lst) == 2

    rotor_names = rotor.rotors_torsion_names(rotors)
    assert rotor_names == (("D5",), ("D12",))

    rotor_axes1 = rotor.rotors_torsion_axes(rotors)
    rotor_axes2 = rotor.rotors_torsion_axes(rotors, "geom")
    assert rotor_axes1 == (((0, 1),), ((9, 11),))
    assert rotor_axes2 == (((0, 1),), ((8, 9),))

    rotor_groups1 = rotor.rotors_torsion_groups(rotors)
    rotor_groups2 = rotor.rotors_torsion_groups(rotors, "geom")
    assert rotor_groups1 == (
        (((2, 3, 4), (5, 6, 7, 9, 11, 12)),),
        (((0, 1, 2, 3, 4, 5, 6, 7), (12,)),),
    )
    assert rotor_groups2 == (
        (((2, 3, 4), (5, 6, 7, 8, 9, 10)),),
        (((0, 1, 2, 3, 4, 5, 6, 7), (10,)),),
    )

    rotor_symms = rotor.rotors_torsion_symmetries(rotors)
    assert rotor_symms == ((3,), (1,))


def test__rotor_multidimensional():
    """test a multi-dimensional hindered rotor"""
    # Handle splitting multirotors when only one rotor of dim <= 4 exists
    rotors = rotor.rotors_from_zmatrix(C6H13OH_ZMA, multi=True)

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom")
    assert len(tors_lst) == 3

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom", flat=True)
    assert len(tors_lst) == 6

    tors_lst = rotor.rotors_torsions(rotors, key_typ="geom", sort=True)
    assert len(tors_lst) == 6

    assert rotor.rotors_torsion_names(rotors) == (
        ("D8", "D14", "D17", "D20"),
        ("D5",),
        ("D11",),
    )
    assert rotor.rotors_torsion_axes(rotors) == (
        ((1, 5), (5, 9), (9, 14), (14, 17)),
        ((0, 1),),
        ((5, 8),),
    )
    assert rotor.rotors_torsion_groups(rotors) == (
        (
            ((0, 2, 3, 4, 6, 7), (8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
            ((0, 1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13), (14, 15, 16, 17, 18, 19, 20)),
            ((0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16), (17, 18, 19, 20)),
            ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18, 19), (20,)),
        ),
        (((2, 3, 4), (5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),),
        (((0, 1, 2, 3, 4, 6, 7, 9, 10, 14, 15, 16, 17, 18, 19, 20), (11, 12, 13)),),
    )
    assert rotor.rotors_torsion_symmetries(rotors) == ((1, 1, 1, 1), (3,), (3,))
    assert rotor.rotors_dimensions(rotors) == (4, 1, 1)

    # Hand splitting multirotors when there is a rotor of dim >= 4
    rotors = rotor.rotors_from_zmatrix(C7H16O2_ZMA, multi=True)
    assert rotor.rotors_torsion_names(rotors) == (
        ("D5",),
        ("D8",),
        ("D11",),
        ("D14",),
        ("D17",),
        ("D20",),
        ("D23",),
        ("D24",),
    )


def test__rotor_multidimensional_custom_grouping():
    """test various ways to do the build"""

    tors_names_lst = (("D5",), ("D11",), ("D8", "D14", "D17", "D20"))
    rotors = rotor.rotors_from_zmatrix(C6H13OH_ZMA, tor_names_lst=tors_names_lst)
    assert rotor.rotors_torsion_names(rotors) == tors_names_lst
    assert rotor.rotors_torsion_axes(rotors) == (
        ((0, 1),),
        ((5, 8),),
        ((1, 5), (5, 9), (9, 14), (14, 17)),
    )
    assert rotor.rotors_torsion_symmetries(rotors) == ((3,), (3,), (1, 1, 1, 1))

    tors_names_lst = (("D5",), ("D11",))
    rotors = rotor.rotors_from_zmatrix(C6H13OH_ZMA, tor_names_lst=tors_names_lst)
    assert rotor.rotors_torsion_names(rotors) == (("D5",), ("D11",))
    assert rotor.rotors_torsion_axes(rotors) == (((0, 1),), ((5, 8),))
    assert rotor.rotors_torsion_symmetries(rotors) == ((3,), (3,))

    tors_names_lst = (("D8", "D14", "D17", "D20"),)
    rotors = rotor.rotors_from_zmatrix(C6H13OH_ZMA, tor_names_lst=tors_names_lst)
    assert rotor.rotors_torsion_names(rotors) == (("D8", "D14", "D17", "D20"),)
    assert rotor.rotors_torsion_axes(rotors) == (((1, 5), (5, 9), (9, 14), (14, 17)),)
    assert rotor.rotors_torsion_symmetries(rotors) == ((1, 1, 1, 1),)


def test__rotor_for_ts():
    """build rotors for a transition state"""

    ts_zgra = reac.ts_graph(C2H5OH_CH3_ZRXN)
    rotors = rotor.rotors_from_zmatrix(C2H5OH_CH3_ZMA, gra=ts_zgra)

    print(rotor.rotors_torsion_names(rotors))
    assert rotor.rotors_torsion_names(rotors) == (("D3",), ("D6",), ("D11",))
    assert rotor.rotors_torsion_axes(rotors) == (((0, 1),), ((1, 3),), ((6, 10),))
    assert rotor.rotors_torsion_groups(rotors) == (
        (((2,), (3, 4, 5, 6, 7, 8, 10, 11, 12, 13)),),
        (((0, 2, 4, 5), (6, 7, 8, 10, 11, 12, 13)),),
        (((0, 1, 2, 3, 4, 5, 7, 8), (11, 12, 13)),),
    )
    assert rotor.rotors_torsion_symmetries(rotors) == ((1,), (1,), (3,))
    assert rotor.rotors_dimensions(rotors) == (1, 1, 1)


def test__torsion_list_string():
    """Test conversion to and from a torsion list string

    (Not a complete serialization -- only serializes the flattened list of torsions)
    """

    rotors = rotor.rotors_from_zmatrix(C3H7OH_ZMA)

    tors_lst = rotor.rotors_torsions(rotors, sort=True)
    assert tors_lst == tors.torsions_from_string(tors.torsions_string(tors_lst))


def test__consistency():
    """check that torsions and z-matrices come out consistent when initializing from data"""
    zma = C4H7O3_ZMA
    tor_lst = tors.torsions_from_string(C4H7O3_TOR_STR)
    rotors = rotor.rotors_from_data(zma=zma, tor_lst=tor_lst)

    tor_coos = rotor.rotors_torsion_coordinates(rotors)
    ref_tor_coos = (
        ((0, 1, 4, 5),),
        ((1, 4, 5, 6),),
        ((4, 5, 6, 9),),
        ((4, 5, 7, 12),),
        ((5, 7, 12, 13),),
    )
    assert tor_coos == ref_tor_coos, f"{tor_coos} != {ref_tor_coos}"

    tor_grps = rotor.rotors_torsion_groups(rotors)
    ref_tor_grps = (
        (((0, 2, 3), (5, 6, 7, 8, 9, 10, 11, 12, 13)),),
        (((0, 1, 2, 3), (6, 7, 8, 9, 10, 11, 12, 13)),),
        (((0, 1, 2, 3, 4, 7, 8, 12, 13), (9, 10, 11)),),
        (((0, 1, 2, 3, 4, 6, 8, 9, 10, 11), (12, 13)),),
        (((0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11), (13,)),),
    )
    assert tor_grps == ref_tor_grps, f"{tor_grps} != {ref_tor_grps}"


def test__dummy_replacement():
    """Check that dummy atoms get replaced upon conversion to geometry."""
    zma = C2H5OH_CH3_ZMA
    zgra = reac.ts_graph(C2H5OH_CH3_ZRXN)
    rotors = rotor.rotors_from_zmatrix(zma, gra=zgra)

    tor_coos = rotor.rotors_torsion_coordinates(rotors, key_typ="geom")
    ref_tor_coos = (((2, 0, 1, 3),), ((0, 1, 3, 6),), ((3, 6, 9, 10),))
    assert tor_coos == ref_tor_coos, f"{tor_coos} != {ref_tor_coos}"


if __name__ == "__main__":
    # test__rotor()
    # test__rotor_with_dummy_atoms()
    # test__rotor_multidimensional()
    # test__torsion_list_string()
    # test__rotor_for_ts()
    # test__consistency()
    test__dummy_replacement()
