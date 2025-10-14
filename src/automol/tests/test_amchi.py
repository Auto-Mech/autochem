""" test automol.amchi
"""

import automol
import numpy
import pytest
from automol import amchi

# Note: Many of these are not canonical AMChIs. I just copied the InChIs over
# and relabeled them.
AR_CHI = "AMChI=1/Ar"
CH4O_CH_CHI = "AMChI=1/CH4O.CH/c1-2;/h2H,1H3;h1H"
CH2O2_CHI = "AMChI=1/CH2O2/c2-1-3/h1-2H/q+1/p+1"
C2H6O_CHI = "AMChI=1/C2H6O/c1-2-3/h3H,2H2,1H3/i2D/t2-/m1/s1"
C2H6O_CHI_NO_STEREO = "AMChI=1/C2H6O/c1-2-3/h3H,2H2,1H3/i2D"

C4H10ZN_CHI = "AMChI=1/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;"
C4H10ZN_CHIS = (
    "AMChI=1/C2H5/c1-2/h1H2,2H3",
    "AMChI=1/C2H5/c1-2/h1H2,2H3",
    "AMChI=1/Zn",
)

C4H5F2O_CHI = "AMChI=1/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1+,4-2+;"
C4H5F2O_CHI_NO_STEREO = "AMChI=1/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H"
C4H5F2O_CHIS = ("AMChI=1/C4H4F2/c5-3-1-2-4-6/h1-4H/b3-1+,4-2+", "AMChI=1/HO/h1H")

C2H2F2_CHI = "AMChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1+"
C2H2F2_CHI_NO_STEREO = "AMChI=1/C2H2F2/c3-1-2-4/h1-2H"
C2H2F2_CHI_STEREO_UNKNOWN = "InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?"

C3H3CL2F3_CHI = "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1"

C2H4F2O2_CHI = "AMChI=1/C2H4F2O2/c3-1(5)2(4)6/h1-2,5-6H/t1-,2-/m0/s1"

C8H13O_CHI = (
    "AMChI=1/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/" "b5-3-,6-4-/t8-/m0/s1"
)
C8H13O_CHI_NO_STEREO = "AMChI=1/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3"

C3H8FNO2_CHI = "InChI=1S/C3H8FNO2/c1-3(4,2-5)7-6/h6H,2,5H2,1H3"
C10H14CLFO_CHI = (
    "AMChI=1/C10H14ClFO/c1-7(9(6-12)10(13)5-11)8-3-2-4-8" "/h2-4,7,9-10,13H,5-6H2,1H3"
)


def test__from_data():
    """test getters"""
    assert AR_CHI == amchi.from_data(
        fml_lyr=amchi.formula_layer(AR_CHI),
    )

    assert CH2O2_CHI == amchi.from_data(
        fml_lyr=amchi.formula_layer(CH2O2_CHI),
        main_lyr_dct=amchi.main_layers(CH2O2_CHI),
        char_lyr_dct=amchi.charge_layers(CH2O2_CHI),
    )

    assert C2H6O_CHI == amchi.from_data(
        fml_lyr=amchi.formula_layer(C2H6O_CHI),
        main_lyr_dct=amchi.main_layers(C2H6O_CHI),
        iso_lyr_dct=amchi.isotope_layers(C2H6O_CHI),
    )

    assert C2H2F2_CHI == amchi.from_data(
        fml_lyr=amchi.formula_layer(C2H2F2_CHI),
        main_lyr_dct=amchi.main_layers(C2H2F2_CHI),
        ste_lyr_dct=amchi.stereo_layers(C2H2F2_CHI),
    )

    assert C8H13O_CHI == amchi.from_data(
        fml_lyr=amchi.formula_layer(C8H13O_CHI),
        main_lyr_dct=amchi.main_layers(C8H13O_CHI),
        ste_lyr_dct=amchi.stereo_layers(C8H13O_CHI),
    )


def test__version():
    """amchi.version"""
    assert amchi.version(C2H2F2_CHI) == "1"


def test__formula_string():
    """amchi.formula_string"""

    assert amchi.formula_layer(AR_CHI) == "Ar"
    assert amchi.formula_layer(CH4O_CH_CHI) == "CH4O.CH"
    assert amchi.formula_layer(CH2O2_CHI) == "CH2O2"
    assert amchi.formula_layer(C2H6O_CHI) == "C2H6O"


def test__formula():
    """amchi.formula"""
    assert amchi.formula(AR_CHI) == {"Ar": 1}
    assert amchi.formula(CH4O_CH_CHI) == {"C": 2, "H": 5, "O": 1}
    assert amchi.formula(CH2O2_CHI) == {"C": 1, "H": 2, "O": 2}
    assert amchi.formula(C2H6O_CHI) == {"C": 2, "H": 6, "O": 1}


def test__has_stereo():
    """test amchi.has_stereo"""
    assert amchi.has_stereo(C2H6O_CHI)
    assert amchi.has_stereo(C4H5F2O_CHI)
    assert amchi.has_stereo(C2H2F2_CHI)
    assert amchi.has_stereo(C8H13O_CHI)
    assert not amchi.has_stereo(AR_CHI)
    assert not amchi.has_stereo(C2H6O_CHI_NO_STEREO)


def test__has_mobile_hydrogens():
    """test amchi.has_mobile_hydrogens"""
    ich1 = "InChI=1S/CHO2/c2-1-3/h(H,2,3)"
    ich2 = "InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1,7H,(H3,6,8,9,10)"
    assert not amchi.has_mobile_hydrogens(C2H6O_CHI)
    assert not amchi.has_mobile_hydrogens(C4H5F2O_CHI)
    assert not amchi.has_mobile_hydrogens(C2H2F2_CHI)
    assert not amchi.has_mobile_hydrogens(C8H13O_CHI)
    assert not amchi.has_mobile_hydrogens(AR_CHI)
    assert not amchi.has_mobile_hydrogens(C2H6O_CHI_NO_STEREO)
    assert amchi.has_mobile_hydrogens(ich1)
    assert amchi.has_mobile_hydrogens(ich2)


def test__symbols():
    """test amchi.symbols"""
    symb_dct = amchi.symbols(C3H8FNO2_CHI)
    print(symb_dct)
    assert symb_dct == {0: "C", 1: "C", 2: "C", 3: "F", 4: "N", 5: "O", 6: "O"}

    symb_dct = amchi.symbols(C10H14CLFO_CHI)
    print(symb_dct)
    assert symb_dct == {
        0: "C",
        1: "C",
        2: "C",
        3: "C",
        4: "C",
        5: "C",
        6: "C",
        7: "C",
        8: "C",
        9: "C",
        10: "Cl",
        11: "F",
        12: "O",
    }


def test__bonds():
    """test amchi.bonds"""
    print(C3H8FNO2_CHI)
    bnds = amchi.bonds(C3H8FNO2_CHI, one_indexed=True)
    print(bnds)
    assert bnds == {
        frozenset({1, 3}),
        frozenset({3, 4}),
        frozenset({3, 2}),
        frozenset({2, 5}),
        frozenset({3, 7}),
        frozenset({7, 6}),
    }

    print(C10H14CLFO_CHI)
    bnds = amchi.bonds(C10H14CLFO_CHI, one_indexed=True)
    print(bnds)
    assert bnds == {
        frozenset({1, 7}),
        frozenset({7, 9}),
        frozenset({7, 8}),
        frozenset({8, 3}),
        frozenset({3, 2}),
        frozenset({2, 4}),
        frozenset({4, 8}),
        frozenset({9, 6}),
        frozenset({9, 10}),
        frozenset({6, 12}),
        frozenset({10, 13}),
        frozenset({10, 5}),
        frozenset({5, 11}),
    }

    print(C3H3CL2F3_CHI)
    bnds = amchi.bonds(C3H3CL2F3_CHI, one_indexed=True)
    print(bnds)
    assert bnds == {
        frozenset({4, 2}),
        frozenset({2, 7}),
        frozenset({2, 1}),
        frozenset({1, 6}),
        frozenset({1, 3}),
        frozenset({3, 5}),
        frozenset({3, 8}),
    }


def test__hydrogen_valences():
    """test amchi.hydrogen_valences"""
    nhyd_dct = amchi.hydrogen_valences(C3H8FNO2_CHI, one_indexed=True)
    print(nhyd_dct)
    assert nhyd_dct == {1: 3, 2: 2, 3: 0, 4: 0, 5: 2, 6: 1, 7: 0}

    nhyd_dct = amchi.hydrogen_valences(C10H14CLFO_CHI, one_indexed=True)
    print(nhyd_dct)
    assert nhyd_dct == {
        1: 3,
        2: 1,
        3: 1,
        4: 1,
        5: 2,
        6: 2,
        7: 1,
        8: 0,
        9: 1,
        10: 1,
        11: 0,
        12: 0,
        13: 1,
    }


def test__atom_stereo_parities():
    """test amchi.atom_stereo_parities"""
    atm_par_dct = amchi.atom_stereo_parities(C2H6O_CHI)
    print(atm_par_dct)
    assert atm_par_dct == {}

    atm_par_dct = amchi.atom_isotope_stereo_parities(C2H6O_CHI)
    print(atm_par_dct)
    assert atm_par_dct == {1: False}

    atm_par_dct = amchi.atom_stereo_parities(C2H4F2O2_CHI)
    print(atm_par_dct)
    assert atm_par_dct == {0: False, 1: False}


def test__bond_stereo_parities():
    """test amchi.bond_stereo_parities"""
    bnd_par_dct = amchi.bond_stereo_parities(C8H13O_CHI)
    assert bnd_par_dct == {frozenset({2, 4}): False, frozenset({3, 5}): False}


def test__is_inverted_enantiomer():
    """test amchi.is_inverted_enantiomer"""
    chi1 = "AMChI=1/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m1/s1"
    chi2 = "AMChI=1/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3/t4-/m0/s1"
    chi3 = "AMChI=1/C4H8O2/c1-3-4(2)6-5/h3-5H,1H2,2H3"
    assert amchi.is_inverted_enantiomer(chi1) is True
    assert amchi.is_inverted_enantiomer(chi2) is False
    assert amchi.is_inverted_enantiomer(chi3) is None


def test__join():
    """test amchi.join"""
    assert amchi.join(amchi.split(C4H10ZN_CHI)) == C4H10ZN_CHI
    assert amchi.join(amchi.split(C4H5F2O_CHI)) == C4H5F2O_CHI

    chi = (
        "AMChI=1/C3H7O4.C2H5FO/c1-3(7-5)2-6-4;1-2(3)4/"
        "h3,5H,2H2,1H3;2,4H,1H3/t3-;2-/m01/s1"
    )
    print(amchi.join(amchi.split(chi)))
    assert amchi.join(amchi.split(chi)) == chi


@pytest.mark.parametrize(
    "chi",
    [
        "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2-,3+",
        "AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-",
        "AMChI=1/C10H14ClFO/c1-8(9(5-12)10(13)6-11)7-3-2-4-7/"
        "h2-4,8-10,13H,5-6H2,1H3",
        "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m0/s1",
        "AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m0/s1",
    ],
)
def test__graph(chi):
    """test amchi.graph"""
    print("---")
    print(chi)
    gra = amchi.graph(chi)

    natms = len(automol.graph.atom_keys(gra))

    for _ in range(5):
        print()
        pmt = list(map(int, numpy.random.permutation(natms)))
        pmt_gra = automol.graph.relabel(gra, dict(enumerate(pmt)))
        pmt_chi = automol.graph.amchi(pmt_gra)
        print(automol.graph.string(pmt_gra))
        print(chi)
        print(pmt_chi)
        assert pmt_chi == chi


@pytest.mark.parametrize(
    "chi0",
    [
        "AMChI=1/CH3.H2O/h1H3;1H2",
        "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2-,3+",
        "AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-",
        "AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m0/s1",
        "AMChI=1/C5H6FO/c6-4-2-1-3-5-7/h1-5,7H/b2-1-,3-1-,4-2-,5-3+",
        "AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m0/s1",
        "AMChI=1/C10H14ClFO/c1-8(9(5-12)10(13)6-11)7-3-2-4-7/"
        "h2-4,8-10,13H,5-6H2,1H3",
    ],
)
def test__smiles(chi0):
    """test amchi.smiles"""
    smi = amchi.smiles(chi0)
    print(smi)
    chi = automol.smiles.amchi(smi)
    print(chi)
    print(chi0)
    assert chi == chi0
    print()


def test__racemic():
    """test amchi.racemic"""
    chi = "AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m0/s1"
    assert (
        automol.amchi.racemic(chi)
        == "AMChI=1/C6H11O/c1-3-4-6-5(2)7-6/h5-6H,1,3-4H2,2H3/t5-,6+/m2/s3"
    )

    assert automol.amchi.racemic(AR_CHI) == AR_CHI
    assert automol.amchi.racemic(C8H13O_CHI_NO_STEREO) == C8H13O_CHI_NO_STEREO


if __name__ == "__main__":
    # test__bonds()
    # test__atom_stereo_parities()
    # test__bond_stereo_parities()
    # test__is_inverted_enantiomer()
    # test__graph()
    test__smiles("AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t1-,2-,3+/m0/s1")
    # test__racemic()
    # test__graph()
    # test__smiles()
    # test__join()
