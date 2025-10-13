""" test automol type conversions
"""

import os
import pandas
import numpy
import automol

PATH = os.path.dirname(os.path.realpath(__file__))

# PROBLEM CASES
# (failing with inchi -> geom conversion sometimes):
# InChI=1/C7H11/c1-3-5-7-6-4-2/h1,3-4,6H,5,7H2,2H3/b3-1+,6-4+
# InChI=1S/C5H5O/c1-2-3-4-5-6/h1-5H/b2-1?,4-3+
# (failing with graph -> inchi conversion sometimes):
# InChI=1S/C4H8O4/c1-4(8-6)2-7-3(4)5/h3,5-6H,2H2,1H3/t3-,4+/m0/s1
# InChI=1S/C5H10O/c1-4-3-5(2)6-4/h4-5H,3H2,1-2H3/t4-,5-/m1/s1
# InChI=1S/C6H12O3/c1-2-5-3-6(9-5)4-8-7/h5-7H,2-4H2,1H3/t5-,6+/m1/s1
# InChI=1S/C6H12O3/c1-4-3-6(8-4)5(2)9-7/h4-7H,3H2,1-2H3/t4-,5-,6-/m0/s1
# InChI=1S/C4H8O3/c1-3-4(7-5)2-6-3/h3-5H,2H2,1H3/t3-,4+/m1/s1
# InChI=1S/C6H12O3/c1-5-4-6(9-5)2-3-8-7/h5-7H,2-4H2,1H3/t5-,6-/m1/s1
# InChI=1S/C5H10O3/c1-2-4-5(8-6)3-7-4/h4-6H,2-3H2,1H3/t4-,5+/m1/s1
# InChI=1S/C5H10O3/c1-4-2-5(8-4)3-7-6/h4-6H,2-3H2,1H3/t4-,5-/m1/s1
# InChI=1S/C5H9O/c1-4-3-5(2)6-4/h4-5H,1,3H2,2H3/t4-,5-/m1/s1
# InChI=1S/C5H10O3/c1-4-5(2,8-6)3-7-4/h4,6H,3H2,1-2H3/t4-,5+/m1/s1
# InChI=1S/C7H14O3/c1-5-3-7(9-5)4-6(2)10-8/h5-8H,3-4H2,1-2H3/t5-,6+,7-/m0/s1
# InChI=1S/C7H14O3/c1-2-6-5-7(10-6)3-4-9-8/h6-8H,2-5H2,1H3/t6-,7+/m1/s1
# InChI=1S/C6H12O3/c1-3-5-6(9-7)4(2)8-5/h4-7H,3H2,1-2H3/t4-,5-,6-/m1/s1
# InChI=1S/C5H10O3/c1-4-5(2-7-4)3-8-6/h4-6H,2-3H2,1H3/t4-,5+/m0/s1
# InChI=1S/C5H10O3/c1-3-5(8-6)4(2)7-3/h3-6H,1-2H3/t3-,4+,5-
# InChI=1S/C7H14O/c1-3-6-5-7(4-2)8-6/h6-7H,3-5H2,1-2H3/t6-,7+
# InChI=1S/C6H12O/c1-3-6-4-5(2)7-6/h5-6H,3-4H2,1-2H3/t5-,6+/m1/s1
# InChI=1S/C7H14O3/c1-2-3-4-6-7(10-8)5-9-6/h6-8H,2-5H2,1H3/t6-,7+/m0/s1


def load_pandas_csv_string_file(path_lst, file_name, path=PATH):
    """Read a file with numpy"""
    file_path = os.path.join(path, *path_lst, file_name)
    file_df = pandas.read_csv(file_path, quotechar="'")

    return file_df


def load_numpy_string_file(path_lst, file_name, path=PATH):
    """Read a file with numpy"""
    file_path = os.path.join(path, *path_lst, file_name)
    file_lst = list(numpy.loadtxt(file_path, dtype=str))

    return file_lst


# InChIs
BS_DF = load_pandas_csv_string_file(["data"], "badspecies.csv", path=PATH)
ICHS_NO_STEREO = load_numpy_string_file(
    ["data"], "heptane_inchis_no_stereo.txt", path=PATH
)
ICHS_WITH_STEREO = load_numpy_string_file(
    ["data"], "heptane_inchis_with_stereo.txt", path=PATH
)
# Use NSAMP = None to test everything
NSAMP = 10
# NSAMP = None

# Geometries
C2H6_H_GEO = (
    ("C", (1.5794198675747746, 0.2458382511130598, -0.0)),
    ("C", (-1.1217055232753832, -0.6604252657772527, -0.0)),
    ("H", (2.9000776470644833, -1.3544800976831528, 7.180959276739747e-05)),
    ("H", (1.9865651407569145, 1.3905530796920174, 1.6742255375628683)),
    ("H", (1.9866142736361765, 1.390449144755117, -1.6742878985250085)),
    ("H", (-1.7333607372068203, -1.636440463684339, -1.7110355127606616)),
    ("H", (-1.7333909728248273, -1.6363724335438226, 1.7110600792002923)),
    ("H", (-2.5967029046096086, 1.4251728623104047, -1.5117809003662624e-05)),
    ("H", (-3.556090402338792, 2.9086399961389326, -1.5117809003662624e-05)),
)

# Z-Matrices
C5H8O_ZMA = (
    ("C", (None, None, None), (None, None, None), (None, None, None)),
    ("C", (0, None, None), ("R1", None, None), (2.894126135733367, None, None)),
    (
        "H",
        (0, 1, None),
        ("R2", "A2", None),
        (2.1002551784038714, 1.9218361502726833, None),
    ),
    (
        "H",
        (0, 1, 2),
        ("R3", "A3", "D3"),
        (2.1008326855365818, 1.932519243247544, 4.182194537966691),
    ),
    (
        "H",
        (0, 1, 2),
        ("R4", "A4", "D4"),
        (2.0985632786105715, 1.9267905467443245, 2.078543469899091),
    ),
    (
        "C",
        (1, 0, 2),
        ("R5", "A5", "D5"),
        (2.7809944542976193, 1.9090367411091194, 1.0264175778482927),
    ),
    (
        "C",
        (1, 0, 5),
        ("R6", "A6", "D6"),
        (2.905905476629344, 1.9462477957117656, 2.1246205025559037),
    ),
    (
        "H",
        (1, 0, 5),
        ("R7", "A7", "D7"),
        (2.1005624282579265, 1.8988979609840009, 4.223022677525844),
    ),
    (
        "X",
        (5, 1, 0),
        ("R8", "A8", "D8"),
        (1.8897261254578288, 1.5707963267948968, 4.216836292856281),
    ),
    (
        "C",
        (5, 8, 1),
        ("R9", "A9", "D9"),
        (2.2778505841014964, 1.5732722619955628, 3.1519457859137696),
    ),
    (
        "X",
        (9, 5, 8),
        ("R10", "A10", "D10"),
        (1.8897261254578286, 1.56832039159423, 0.0),
    ),
    (
        "H",
        (9, 10, 5),
        ("R11", "A11", "D11"),
        (2.0003442808863467, 1.57550991099349, 3.14859478950736),
    ),
    (
        "O",
        (6, 1, 0),
        ("R12", "A12", "D12"),
        (2.65076899334649, 1.9387190313618887, 1.0262708014428483),
    ),
    (
        "H",
        (6, 1, 12),
        ("R13", "A13", "D13"),
        (2.1058345184525726, 1.9323237957467607, 2.129177885999989),
    ),
    (
        "H",
        (6, 1, 12),
        ("R14", "A14", "D14"),
        (2.1010240316411886, 1.9207088798352128, 4.1894956154070275),
    ),
    (
        "H",
        (12, 6, 1),
        ("R15", "A15", "D15"),
        (1.8758293656194, 1.8624105681328567, 1.2477273554765336),
    ),
)

HOOH_ZMA_C2 = (
    ("O", (None, None, None), (None, None, None), (None, None, None)),
    ("O", (0, None, None), ("R1", None, None), (2.747759350307364, None, None)),
    (
        "H",
        (0, 1, None),
        ("R2", "A2", None),
        (1.8445107263893656, 1.6890062361073546, None),
    ),
    (
        "H",
        (1, 0, 2),
        ("R3", "A3", "D3"),
        (1.8445105807874629, 1.6890041579036852, 2.2578840834994196),
    ),
)
HOOH_ZMA_CS = (
    ("O", (None, None, None), (None, None, None), (None, None, None)),
    ("O", (0, None, None), ("R1", None, None), (2.747759350307364, None, None)),
    (
        "H",
        (0, 1, None),
        ("R2", "A2", None),
        (1.8445107263893656, 1.6890062361073546, None),
    ),
    (
        "H",
        (1, 0, 2),
        ("R3", "A3", "D3"),
        (1.8445105807874629, 1.6890041579036852, 0.0000),
    ),
)


def test__geom__with_stereo():
    """test geom conversions"""
    print("geom with stereo")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        ref_geo = automol.chi.geometry(ref_ich)
        ref_chi = automol.geom.chi(ref_geo)

        # Test ChI conversion
        print(ref_chi, flush=True)
        geo = automol.chi.geometry(ref_chi)
        chi = automol.geom.chi(geo)
        assert chi == ref_chi
        assert automol.geom.formula(geo) == automol.chi.formula(chi)

        # Test z-matrix conversion
        zma = automol.geom.zmatrix(ref_geo)
        geo = automol.zmat.geometry(zma)
        chi = automol.geom.chi(geo)
        assert chi == ref_chi

        # Test symmetry factor calculations
        sym_num = automol.geom.external_symmetry_factor(geo)
        print("symmetry number:", sym_num)

        end_sym_fac = automol.symm.end_group_symmetry_factor(geo)
        print("end-group symmetry factor:", end_sym_fac)


def test__graph__with_stereo():
    """test graph conversions"""
    print("graph with stereo")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = list(numpy.random.choice(ref_ichs, NSAMP))

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        ref_geo = automol.chi.geometry(ref_ich)
        ref_chi = automol.geom.chi(ref_geo)

        print(ref_chi, flush=True)
        gra = automol.chi.graph(ref_chi)
        chi = automol.graph.chi(gra, stereo=True)
        assert chi == ref_chi

        assert automol.graph.formula(gra) == automol.chi.formula(chi)

        # If this is a good InChI, test RDKit molecule conversion as well
        if not automol.graph.inchi_is_bad(gra, ref_ich):
            print("Good InChI! Testing RDKit molecule...")
            mol = automol.extern.rdkit_.from_graph(gra, stereo=True)
            ich = automol.extern.rdkit_.to_inchi(mol)
            print(ich, flush=True)
            assert ich == ref_ich

            gra_out = automol.extern.rdkit_.to_graph(mol)
            print(gra_out, flush=True)
            assert gra == gra_out


def test__smiles__with_stereo():
    """test smiles conversions"""
    print("smiles from geom")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        ref_geo = automol.chi.geometry(ref_ich)
        ref_chi = automol.geom.chi(ref_geo)

        smi = automol.chi.smiles(ref_chi)
        chi = automol.smiles.chi(smi)
        assert chi == ref_chi, f"{chi} != {ref_chi}"


def test__graph__misc():
    """test graph conversions"""

    ref_ich = "InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H"
    print(ref_ich)
    gra = automol.chi.graph(ref_ich)
    ich = automol.graph.chi(gra)
    assert ich == ref_ich

    ref_ich = "InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1-,4-2-;"
    print(ref_ich)
    gra = automol.chi.graph(ref_ich)
    ich = automol.graph.chi(gra, stereo=True)
    assert ich == ref_ich

    ref_ich = "InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+"
    print(ref_ich)
    geo = automol.chi.geometry(ref_ich)
    gra = automol.geom.graph(geo)
    ich = automol.graph.chi(gra, stereo=True)
    assert ich == ref_ich

    ref_ich = "InChI=1S/C2H4O/c1-2-3/h2-3H,1H2"
    print(ref_ich)
    ref_conn_gra = (
        {
            0: ("C", 0, None),
            1: ("C", 0, None),
            2: ("O", 0, None),
            3: ("H", 0, None),
            4: ("H", 0, None),
            5: ("H", 0, None),
            6: ("H", 0, None),
        },
        {
            frozenset({0, 1}): (1, None),
            frozenset({0, 3}): (1, None),
            frozenset({0, 4}): (1, None),
            frozenset({1, 2}): (1, None),
            frozenset({1, 5}): (1, None),
            frozenset({2, 6}): (1, None),
        },
    )
    conn_gra = automol.geom.graph_without_stereo(automol.chi.geometry(ref_ich))
    assert conn_gra == ref_conn_gra

    ich = "InChI=1S/C2HCl/c1-2-3/h1H"
    print(ich)
    gra = automol.chi.graph(ich)
    geo = automol.graph.geometry(gra)
    assert ich == automol.geom.inchi(geo)

    # Test stereo
    def randomize_atom_ordering(geo):
        """randomize atom ordering in a geometry"""
        natms = automol.geom.count(geo)
        ord_dct = dict(enumerate(numpy.random.permutation(natms)))
        return automol.geom.reorder(geo, ord_dct)

    # smi = 'FC=C(C=CC=CF)C=CC=CF'
    # smi = 'FC=CC=CC=CF'
    # ich = automol.smiles.chi('CC([O])=CCO')
    ich = automol.smiles.chi("CC([O])=CCO.O")
    geo = automol.chi.geometry(ich)
    geo = randomize_atom_ordering(geo)
    gra = automol.geom.graph(geo)
    ich = automol.graph.chi(gra, stereo=True)
    print(ich)
    assert ich in (
        "AMChI=1/C4H7O2.H2O/c1-4(5)2-3-6;/h2,6H,3H2,1H3;1H2/b4-2-;",
        "AMChI=1/C4H7O2.H2O/c1-4(5)2-3-6;/h2,6H,3H2,1H3;1H2/b4-2+;",
    )

    # EXTRA TEST CASE 1
    gra = automol.smiles.graph(r"[H]/N=C\C(\C=N\[H])=N/[H]")

    natms = len(automol.graph.atoms(gra))
    for _ in range(10):
        pmt_dct = dict(enumerate(map(int, numpy.random.permutation(natms))))
        pmt_gra = automol.graph.relabel(gra, pmt_dct)
        pmt_geo = automol.graph.geometry(pmt_gra)
        pmt_gra_ = automol.geom.graph(pmt_geo)
        print(automol.graph.string(pmt_gra))
        print()
        print(automol.graph.string(pmt_gra_))
        assert pmt_gra == pmt_gra_
        print("------------------------------------")

    # EXTRA TEST CASE 2
    gra = automol.smiles.graph(r"F[C@H]([C@H](Cl)F)[C@@H](Cl)F")

    natms = len(automol.graph.atoms(gra))
    for _ in range(10):
        pmt_dct = dict(enumerate(map(int, numpy.random.permutation(natms))))
        pmt_gra = automol.graph.relabel(gra, pmt_dct)
        pmt_geo = automol.graph.geometry(pmt_gra)
        pmt_gra_ = automol.geom.graph(pmt_geo)
        print(automol.graph.string(pmt_gra))
        print()
        print(automol.graph.string(pmt_gra_))
        assert pmt_gra == pmt_gra_
        print("------------------------------------")

    # test graph => geom conversion with non-standard keys
    gra = (
        {
            5: ("C", 0, None),
            6: ("C", 0, None),
            7: ("H", 0, None),
            8: ("H", 0, None),
            9: ("H", 0, None),
            10: ("C", 0, None),
            11: ("H", 0, None),
            12: ("H", 0, None),
            13: ("C", 0, None),
            14: ("H", 0, None),
            15: ("H", 0, None),
            16: ("C", 0, None),
            17: ("H", 0, None),
            18: ("C", 0, None),
            19: ("H", 0, None),
            20: ("C", 0, None),
            21: ("H", 0, None),
            22: ("C", 0, None),
            23: ("H", 0, None),
            24: ("H", 0, None),
            25: ("C", 0, None),
            26: ("H", 0, None),
            27: ("H", 0, None),
            28: ("H", 0, None),
            29: ("H", 0, None),
            30: ("H", 0, None),
        },
        {
            frozenset({10, 15}): (1, None),
            frozenset({25, 28}): (1, None),
            frozenset({25, 22}): (1, None),
            frozenset({10, 13}): (1, None),
            frozenset({11, 6}): (1, None),
            frozenset({9, 5}): (1, None),
            frozenset({24, 20}): (1, None),
            frozenset({25, 30}): (1, None),
            frozenset({5, 6}): (1, None),
            frozenset({16, 13}): (1, True),
            frozenset({27, 22}): (1, None),
            frozenset({12, 6}): (1, None),
            frozenset({8, 5}): (1, None),
            frozenset({25, 29}): (1, None),
            frozenset({18, 20}): (1, None),
            frozenset({16, 18}): (1, False),
            frozenset({17, 13}): (1, None),
            frozenset({26, 22}): (1, None),
            frozenset({5, 7}): (1, None),
            frozenset({18, 21}): (1, None),
            frozenset({10, 6}): (1, None),
            frozenset({10, 14}): (1, None),
            frozenset({16, 19}): (1, None),
            frozenset({20, 22}): (1, None),
            frozenset({20, 23}): (1, None),
        },
    )
    geo = automol.graph.geometry(gra)
    gra_ = automol.geom.graph(geo)
    iso = automol.graph.isomorphism(gra, gra_, stereo=True)
    print(iso)
    assert iso

    # EXTRA TEST CASE 3
    smi = "c1ccc2c(c1)Cc1c2cccc1"
    ich = automol.smiles.inchi(smi)
    gra = automol.smiles.graph(smi)
    chi = automol.graph.chi(gra)
    assert ich == chi

    # EXTRA TEST CASE 3
    smi = "[CH]=C=C"
    ich = automol.smiles.inchi(smi)
    gra = automol.smiles.graph(smi)
    chi = automol.graph.chi(gra)
    assert ich == chi

    # EXTRA TEST CASE 4
    gra = (
        {
            0: ("C", 1, None),
            2: ("C", 0, None),
            4: ("C", 0, None),
            5: ("H", 0, None),
            6: ("H", 0, None),
        },
        {
            frozenset({4, 5}): (1, None),
            frozenset({2, 4}): (2, None),
            frozenset({0, 2}): (2, None),
            frozenset({4, 6}): (1, None),
        },
    )
    assert automol.graph.chi(gra) == "InChI=1S/C3H3/c1-3-2/h1H,2H2"

    # EXTRA TEST CASE 5
    ref_ich = "InChI=1S/C4H6O2/c1-3-4(2,5-3)6-3/h1-2H3"
    geo = automol.chi.geometry(ref_ich)
    zma = automol.geom.zmatrix(geo)
    geo = automol.zmat.geometry(zma)
    ich = automol.geom.inchi(geo, stereo=False)
    assert ich == ref_ich, f"{ich} != {ref_ich}"


def test__inchi_geometry():
    """test automol.chi.geometry"""
    ref_ich = "InChI=1S/H2S/h1H2"
    ich = automol.geom.chi(automol.chi.geometry(ref_ich))
    print(ich)
    assert ich == ref_ich

    ref_ich = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    ich = automol.geom.chi(automol.chi.geometry(ref_ich))
    print(ich)
    assert ich == ref_ich

    ref_ich = "InChI=1S/Ar"
    ich = automol.geom.chi(automol.chi.geometry(ref_ich))
    print(ich)
    assert ich == ref_ich

    ref_ich = "InChI=1S/Cl2/c1-2"
    ich = automol.geom.chi(automol.chi.geometry(ref_ich))
    print(ich)
    assert ich == ref_ich

    ich = "InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3"
    geo = automol.chi.geometry(ich)
    ich = automol.geom.chi(geo, stereo=True)
    print(automol.geom.chi(geo, stereo=True))
    print(ich)
    assert ich in (
        "AMChI=1/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5+",
        "AMChI=1/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5-",
    )

    ich = "InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3"
    geo = automol.chi.geometry(ich)
    ich = automol.geom.chi(geo, stereo=True)
    print(ich)
    assert ich in (
        "AMChI=1/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5+",
        "AMChI=1/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5-",
    )

    # Extra test case for broken InChI conversion
    ich = automol.smiles.chi("CC([O])=CCO")
    geo = automol.chi.geometry(ich)
    ich = automol.geom.chi(geo, stereo=True)
    print(ich)
    assert ich in (
        "AMChI=1/C4H7O2/c1-4(5)2-3-6/h2,6H,3H2,1H3/b4-2+",
        "AMChI=1/C4H7O2/c1-4(5)2-3-6/h2,6H,3H2,1H3/b4-2-",
    )

    # Extra test case for broken InChI conversion
    ich = automol.smiles.chi("CC([O])=CCO.[OH]")
    geo = automol.chi.geometry(ich)
    ich = automol.geom.chi(geo, stereo=True)
    print(ich)
    assert ich in (
        "AMChI=1/C4H7O2.HO/c1-4(5)2-3-6;/h2,6H,3H2,1H3;1H/b4-2+;",
        "AMChI=1/C4H7O2.HO/c1-4(5)2-3-6;/h2,6H,3H2,1H3;1H/b4-2-;",
    )


def test__multiple_rings():
    """test graph => inchi conversion for multiple rings"""
    ref_chi = (
        "AMChI=1/C17H19NO3/c1-18-7-6-17-11-9-2-3-10(19)12(11)21-16"
        "(17)14(20)5-4-13(17)15(18)8-9/h2-5,13-16,19-20H,6-8H2,1H3/"
        "t13-,14-,15-,16+,17-,18-/m1/s1"
    )
    gra = automol.chi.graph(ref_chi)
    print(f"gra = {gra}")
    chi = automol.graph.chi(gra, stereo=True)
    print(ref_chi)
    print(chi)
    assert chi == ref_chi


def test__weird_valencies():
    """atoms with more than two unpaired electrons used to cause problems"""
    ich = "InChI=1S/C"
    geo = (("C", (0.0, 0.0, 0.0)),)
    gra = ({0: ("C", 0, None)}, {})
    smi = "[C]"

    assert automol.smiles.inchi(smi) == ich
    assert automol.geom.inchi(geo) == ich
    assert automol.graph.inchi(gra) == ich

    ich = "InChI=1S/CF/c1-2"
    geo = (("C", (0.0, 0.0, 0.0)), ("F", (0.0, 0.0, 2.4)))
    gra = ({0: ("C", 0, None), 1: ("F", 0, None)}, {frozenset({0, 1}): (1, None)})
    smi = "[C]F"

    assert automol.smiles.inchi(smi) == ich
    assert automol.geom.inchi(geo) == ich
    assert automol.graph.inchi(gra) == ich


def test__symmetry_removal():
    """make sure the geometries we generate avoid planar symmetry"""
    geo = automol.smiles.geometry(r"F/C=C/F")
    gra = automol.geom.graph(geo)
    cis_dihs = automol.graph.geometry_dihedrals_near_value(
        gra, geo, 0.0, tol=4.9, degree=True
    )
    print(cis_dihs)
    for dih in cis_dihs:
        print(automol.geom.dihedral_angle(geo, *dih, degree=True))
    assert not cis_dihs

    trans_dihs = automol.graph.geometry_dihedrals_near_value(
        gra, geo, 180.0, tol=4.9, degree=True
    )
    print(trans_dihs)
    for dih in trans_dihs:
        print(automol.geom.dihedral_angle(geo, *dih, degree=True))
    assert not trans_dihs


if __name__ == "__main__":
    test__geom__with_stereo()
    test__graph__with_stereo()
    test__smiles__with_stereo()
    # test__graph__misc()
    # test__inchi_geometry()
    # test__multiple_rings()
    # test__weird_valencies()
    # test__symmetry_removal()
