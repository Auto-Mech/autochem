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


def test__inchi__geometry():
    """ test automol.inchi.geometry
    """
    ref_ich = 'InChI=1S/H2S/h1H2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich

    ref_ich = 'InChI=1S/Cl2/c1-2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich


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
    # test__zmatrix__with_stereo()
    # test__smiles__with_stereo()
    # test__geom__with_stereo()
    # test__graph__with_stereo()
    # test__graph__misc()
    # test__geom__with_stereo()
    # test__geom__zmatrix()
    # test__geom__zmatrix_torsion_coordinate_names()
    # test__multiple_rings()
    test__inchi__geometry()
