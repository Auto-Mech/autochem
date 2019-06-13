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

NSAMP = 100


def test__geom__with_stereo():
    """ test geom conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        geo = automol.inchi.geometry(ref_ich)
        ich = automol.geom.inchi(geo)
        print(ich, ref_ich)
        assert ich == ref_ich

        assert automol.geom.formula(geo) == automol.inchi.formula(ich)


def test__zmatrix__with_stereo():
    """ test zmatrix conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        ref_geo = automol.inchi.geometry(ref_ich)
        zma = automol.geom.zmatrix(ref_geo)
        geo = automol.zmatrix.geometry(zma)
        ich = automol.geom.inchi(geo)
        print(ich, ref_ich)
        assert ich == ref_ich

        assert automol.zmatrix.formula(zma) == automol.inchi.formula(ich)


def test__smiles__with_stereo():
    """ test smiles conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        gra = automol.inchi.smiles(ref_ich)
        ich = automol.smiles.inchi(gra)
        print(ich, ref_ich)
        assert ich == ref_ich


def test__graph__with_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        gra = automol.inchi.graph(ref_ich)
        ich = automol.graph.inchi(gra)
        print(ich, ref_ich)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__graph__no_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        gra = automol.inchi.graph(ref_ich)
        gra = automol.graph.without_stereo_parities(gra)
        ich = automol.graph.inchi(gra)
        print(ich, ref_ich)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__graph__misc():
    """ test graph conversions
    """
    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra)
    print(ich, ref_ich)
    assert ich == ref_ich

    assert automol.graph.formula(gra) == automol.inchi.formula(ich)

    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1-,4-2-;'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra)
    print(ich, ref_ich)
    assert ich == ref_ich

    assert automol.graph.formula(gra) == automol.inchi.formula(ich)


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


if __name__ == '__main__':
    # test__geom__graph()
    # test__geom__inchi()
    # test__zmatrix__with_stereo()
    # test__geom__with_stereo()
    # test__graph__no_stereo()
    # test__graph__with_stereo()
    test__graph__misc()
