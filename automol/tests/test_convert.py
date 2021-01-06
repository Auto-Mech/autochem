""" test automol.automol.convert
"""
import os
import numpy
import automol

PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
ICHS_NO_STEREO = numpy.loadtxt(
    os.path.join(DATA_PATH, 'heptane_inchis_no_stereo.txt'), dtype=str)
ICHS_WITH_STEREO = numpy.loadtxt(
    os.path.join(DATA_PATH, 'heptane_inchis_with_stereo.txt'), dtype=str)

NSAMP = 50


def test__geom__with_stereo():
    """ test geom conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print(ref_ich)
        geo = automol.inchi.geometry(ref_ich)
        ich = automol.geom.inchi(geo)
        print(ich)
        print()
        assert ich == ref_ich

        assert automol.geom.formula(geo) == automol.inchi.formula(ich)


def test__graph__with_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print(ref_ich)
        gra = automol.inchi.graph(ref_ich)
        ich = automol.graph.inchi(gra)
        print(ich)
        print()
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__graph__no_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print(ref_ich)
        gra = automol.inchi.graph(ref_ich)
        gra = automol.graph.without_stereo_parities(gra)
        ich = automol.graph.inchi(gra)
        print(ich)
        print()
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__zmatrix__with_stereo():
    """ test zmatrix conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print(ref_ich)
        ref_geo = automol.inchi.geometry(ref_ich)
        zma = automol.geom.zmatrix(ref_geo)
        geo = automol.zmat.geometry(zma)
        ich = automol.geom.inchi(geo)
        assert ich == ref_ich

        assert automol.zmat.formula(zma) == automol.inchi.formula(ich)


def test__smiles__with_stereo():
    """ test smiles conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        smi = automol.inchi.smiles(ref_ich)
        ich = automol.smiles.inchi(smi)
        assert ich == ref_ich


def test__graph__misc():
    """ test graph conversions
    """
    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra)
    print(ich, ref_ich)
    assert ich == ref_ich

    # assert automol.graph.formula(gra) == automol.inchi.formula(ich)

    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1-,4-2-;'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra)
    print(ich, ref_ich)
    assert ich == ref_ich

    # assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__inchi_geometry():
    """ test automol.inchi.geometry
    """
    ref_ich = 'InChI=1S/H2S/h1H2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich

    ref_ich = 'InChI=1S/Cl2/c1-2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich


def test__geom__zmatrix():
    """ test automol.geom.zmatrix
    """
    geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))
    zma = automol.geom.zmatrix(geo)

    ref_zma = (
        (('C', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('R1', None, None)),
         ('H', (0, 1, None), ('R2', 'A2', None)),
         ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
         ('H', (0, 1, 2), ('R4', 'A4', 'D4')),
         ('H', (1, 0, 2), ('R5', 'A5', 'D5'))),
        {'R1': 2.67535, 'R2': 2.06501, 'A2': 1.9116242,
         'R3': 2.06501, 'A3': 1.9116242, 'D3': 2.108497362,
         'R4': 2.06458, 'A4': 1.9020947, 'D4': 4.195841334,
         'R5': 1.83748, 'A5': 1.8690905, 'D5': 5.228936625})
    assert automol.zmatrix.almost_equal(zma, ref_zma)

    # CH2O + OH Abstraction
    syms = ['C', 'O', 'H', 'H', 'O', 'H']
    xyzs = [[0.00000, 0.00000, 0.00000],
            [0.00000, 0.00000, 1.20170],
            [0.93428, 0.00000, -0.59285],
            [-0.95532, 0.00002, -0.58901],
            [-2.76283, 0.00158, -0.68933],
            [-2.90953, 0.00149, 0.27760]]
    ts_bnds = [{0, 3}, {3, 4}]
    geo = automol.geom.from_data(syms, xyzs, angstrom=True)

    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)

    ref_zma = (
        (('C', (None, None, None), (None, None, None)),
         ('O', (0, None, None), ('R1', None, None)),
         ('H', (0, 1, None), ('R2', 'A2', None)),
         ('H', (0, 1, 2), ('R3', 'A3', 'D3')),
         ('O', (3, 0, 1), ('R4', 'A4', 'D4')),
         ('H', (4, 3, 0), ('R5', 'A5', 'D5'))),
        {'R1': 2.27088, 'R2': 2.09099, 'A2': 2.1362306445634998,
         'R3': 2.12085, 'A3': 2.1232977548062215, 'D3': 3.1416101068823132,
         'R4': 3.42096, 'A4': 2.6445228826218083, 'D4': 0.001763341049874911,
         'R5': 1.84814, 'A5': 1.6659202616870914, 'D5': 6.281457431220113})
    assert automol.zmatrix.almost_equal(zma, ref_zma)

    # Alkoxy Rearrangment
    syms = ['O', 'C', 'C', 'H', 'H', 'O', 'H', 'H', 'N', 'O', 'O']
    xyzs = [[0.00000, 0.00000, 0.00000],
            [0.00000, 0.00000, 1.43013],
            [1.49249, 0.00000, 1.42455],
            [-0.45693, -0.89518, 1.85040],
            [-0.44376, 0.90025, 1.85542],
            [2.04097, 1.16894, 1.92249],
            [1.27109, 0.01640, 0.11377],
            [2.03282, -0.91188, 1.65781],
            [3.41118, 1.36156, 1.58786],
            [3.93221, 0.51367, 0.92911],
            [3.84108, 2.37118, 2.03747]]
    ts_bnds = [{0, 6}, {2, 6}]
    geo = automol.geom.from_data(syms, xyzs, angstrom=True)

    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)

    ref_zma = (
        (('O', (None, None, None), (None, None, None)),
         ('C', (0, None, None), ('R1', None, None)),
         ('H', (0, 1, None), ('R2', 'A2', None)),
         ('C', (1, 0, 2), ('R3', 'A3', 'D3')),
         ('H', (1, 0, 3), ('R4', 'A4', 'D4')),
         ('H', (1, 0, 3), ('R5', 'A5', 'D5')),
         ('O', (3, 1, 0), ('R6', 'A6', 'D6')),
         ('H', (3, 1, 6), ('R7', 'A7', 'D7')),
         ('N', (6, 3, 1), ('R8', 'A8', 'D8')),
         ('O', (8, 6, 3), ('R9', 'A9', 'D9')),
         ('O', (8, 6, 9), ('R10', 'A10', 'D10'))),
        {'R1': 2.70255, 'R2': 2.41181, 'A2': 1.4815349528601507,
         'R3': 2.82042, 'A3': 1.56705783153712, 'D3': 6.270287324007349,
         'R4': 2.05864, 'A4': 1.966863893950, 'D4': 4.240434497352903,
         'R5': 2.05993, 'A5': 1.971593736223, 'D5': 2.028770722518209,
         'R6': 2.6152, 'A6': 1.9768471772714, 'D6': 1.974979674971743,
         'R7': 2.05093, 'A7': 2.090991710352, 'D7': 2.484912522526927,
         'R8': 2.69016, 'A8': 1.994544815887, 'D8': 3.437740121068191,
         'R9': 2.25531, 'A9': 2.041843238616, 'D9': 6.28194612341067,
         'R10': 2.24098, 'A10': 1.955205094547, 'D10': 3.1468286413458})
    assert automol.zmatrix.almost_equal(zma, ref_zma)

    # EGDN Beta Scission
    syms = ['O', 'C', 'C', 'H', 'H', 'O', 'H', 'H', 'N', 'O', 'O']
    xyzs = [[0.00000, 0.00000, 0.00000],
            [0.00000, 0.00000, 1.22745],
            [2.06399, 0.00000, 1.71761],
            [-0.16514, -0.93179, 1.79642],
            [-0.15635, 0.93075, 1.80036],
            [2.21856, 0.18590, 3.05580],
            [2.35175, 0.89622, 1.18872],
            [2.33417, -0.96682, 1.31894],
            [2.15304, -1.00873, 3.85747],
            [2.25822, -0.76166, 5.01008],
            [2.01188, -2.04079, 3.28024]]
    ts_bnds = [{1, 2}]
    geo = automol.geom.from_data(syms, xyzs, angstrom=True)

    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)
    ref_zma = (
        (('O', (None, None, None), (None, None, None)),
         ('C', (0, None, None), ('R1', None, None)),
         ('C', (1, 0, None), ('R2', 'A2', None)),
         ('H', (1, 0, 2), ('R3', 'A3', 'D3')),
         ('H', (1, 0, 2), ('R4', 'A4', 'D4')),
         ('O', (2, 1, 0), ('R5', 'A5', 'D5')),
         ('H', (2, 1, 5), ('R6', 'A6', 'D6')),
         ('H', (2, 1, 5), ('R7', 'A7', 'D7')),
         ('N', (5, 2, 1), ('R8', 'A8', 'D8')),
         ('O', (8, 5, 2), ('R9', 'A9', 'D9')),
         ('O', (8, 5, 9), ('R10', 'A10', 'D10'))),
        {'R1': 2.31954, 'R2': 4.00885, 'A2': 1.803954861568819,
         'R3': 2.08661, 'A3': 2.112127647593458, 'D3': 4.53698339055926,
         'R4': 2.08639, 'A4': 2.116368797675804, 'D4': 1.73722568827732,
         'R5': 2.56975, 'A5': 1.9155512139413364, 'D5': 2.9958227544632265,
         'R6': 2.04033, 'A6': 1.717445871864468, 'D6': 4.282933264638965,
         'R7': 2.04115, 'A7': 1.7295497302270488, 'D7': 2.152444753314527,
         'R8': 2.72154, 'A8': 2.014406662774295, 'D8': 1.4993233485964768,
         'R9': 2.23645, 'A9': 1.944331693306723, 'D9': 3.176097812901721,
         'R10': 2.2505, 'A10': 2.0425413703164437, 'D10': 3.145048405508742})
    assert automol.zmatrix.almost_equal(zma, ref_zma)

    # HNCO + OH Addition
    syms = ['N', 'C', 'H', 'O', 'O', 'H']
    xyzs = [[0.00000, 0.00000, 0.00000],
            [0.00000, 0.00000, 1.26664],
            [0.90194, 0.00000, -0.45759],
            [-0.36796, 0.29063, 2.33250],
            [1.12388, -1.51240, 1.20161],
            [0.50305, -2.22859, 1.00044]]
    ts_bnds = [{1, 4}]
    geo = automol.geom.from_data(syms, xyzs, angstrom=True)

    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)
    ref_zma = (
        (('N', (None, None, None), (None, None, None)),
         ('C', (0, None, None), ('R1', None, None)),
         ('H', (0, 1, None), ('R2', 'A2', None)),
         ('O', (1, 0, 2), ('R3', 'A3', 'D3')),
         ('O', (1, 0, 3), ('R4', 'A4', 'D4')),
         ('H', (4, 1, 0), ('R5', 'A5', 'D5'))),
        {'R1': 2.3936, 'R2': 1.91123, 'A2': 2.0402898955813713,
         'R3': 2.20047, 'A3': 2.7271467694112195, 'D3': 2.4730791901984053,
         'R4': 3.56287, 'A4': 1.5362981487999765, 'D4': 2.878397002389048,
         'R5': 1.83101, 'A5': 1.790707812546182, 'D5': 4.919018510528298})
    assert automol.zmatrix.almost_equal(zma, ref_zma)


def test__geom__zmatrix_torsion_coordinate_names():
    """ test automol.geom.zmatrix_torsion_coordinate_names
    """
    geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))
    zma = automol.geom.zmatrix(geo)

    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    assert set(tors_names) <= set(automol.zmatrix.dihedral_angle_names(zma))


def test__multiple_rings():
    """ test graph => inchi conversion for multiple rings
    """
    ref_ich = ('InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19'
               ')4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/'
               't10-,11-,13-,16+,17+/m1/s1')
    gra = automol.graph.from_string("""
        atoms:
          1: {symbol: C, implicit_hydrogen_valence: 3, stereo_parity: null}
          2: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
          3: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
          4: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
          5: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: null}
          6: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
          7: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
          8: {symbol: C, implicit_hydrogen_valence: 2, stereo_parity: null}
          9: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          10: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: false}
          11: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: true}
          12: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          13: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: false}
          14: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          15: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
          16: {symbol: C, implicit_hydrogen_valence: 1, stereo_parity: false}
          17: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: true}
          18: {symbol: N, implicit_hydrogen_valence: 0, stereo_parity: null}
          19: {symbol: O, implicit_hydrogen_valence: 1, stereo_parity: null}
          20: {symbol: O, implicit_hydrogen_valence: 1, stereo_parity: null}
          21: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
        bonds:
          1-18: {order: 1, stereo_parity: null}
          2-4: {order: 1, stereo_parity: null}
          2-9: {order: 1, stereo_parity: null}
          3-5: {order: 1, stereo_parity: null}
          3-10: {order: 1, stereo_parity: null}
          4-12: {order: 1, stereo_parity: null}
          5-13: {order: 1, stereo_parity: null}
          6-7: {order: 1, stereo_parity: null}
          6-17: {order: 1, stereo_parity: null}
          7-18: {order: 1, stereo_parity: null}
          8-9: {order: 1, stereo_parity: null}
          8-11: {order: 1, stereo_parity: null}
          9-14: {order: 1, stereo_parity: null}
          10-11: {order: 1, stereo_parity: null}
          10-17: {order: 1, stereo_parity: null}
          11-18: {order: 1, stereo_parity: null}
          12-15: {order: 1, stereo_parity: null}
          12-19: {order: 1, stereo_parity: null}
          13-16: {order: 1, stereo_parity: null}
          13-20: {order: 1, stereo_parity: null}
          14-15: {order: 1, stereo_parity: null}
          14-17: {order: 1, stereo_parity: null}
          15-21: {order: 1, stereo_parity: null}
          16-17: {order: 1, stereo_parity: null}
          16-21: {order: 1, stereo_parity: null}
    """)
    ich = automol.graph.inchi(gra)
    assert ich == ref_ich


if __name__ == '__main__':
    # test__geom__graph()
    # test__geom__inchi()
    # test__graph__no_stereo()
    # test__graph__with_stereo()
    test__zmatrix__with_stereo()
    # test__smiles__with_stereo()
    # test__geom__with_stereo()
    # test__graph__with_stereo()
    # test__graph__misc()
    # test__geom__with_stereo()
    # test__geom__zmatrix()
    # test__geom__zmatrix_torsion_coordinate_names()
    # test__multiple_rings()
