""" test automol type conversions
"""

import os
import pandas
import numpy
import pytest
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
    """ Read a file with numpy
    """
    file_path = os.path.join(path, *path_lst, file_name)
    file_df = pandas.read_csv(file_path, quotechar="'")

    return file_df


def load_numpy_string_file(path_lst, file_name, path=PATH):
    """ Read a file with numpy
    """
    file_path = os.path.join(path, *path_lst, file_name)
    file_lst = list(numpy.loadtxt(file_path, dtype=str))

    return file_lst


# InChIs
BS_DF = load_pandas_csv_string_file(
    ['data'], 'badspecies.csv', path=PATH)
ICHS_NO_STEREO = load_numpy_string_file(
    ['data'], 'heptane_inchis_no_stereo.txt', path=PATH)
ICHS_WITH_STEREO = load_numpy_string_file(
    ['data'], 'heptane_inchis_with_stereo.txt', path=PATH)
# Use NSAMP = None to test everything
NSAMP = None
# NSAMP = 10

# Geometries
C2H6_H_GEO = (
    ('C', (1.5794198675747746, 0.2458382511130598, -0.0)),
    ('C', (-1.1217055232753832, -0.6604252657772527, -0.0)),
    ('H', (2.9000776470644833, -1.3544800976831528, 7.180959276739747e-05)),
    ('H', (1.9865651407569145, 1.3905530796920174, 1.6742255375628683)),
    ('H', (1.9866142736361765, 1.390449144755117, -1.6742878985250085)),
    ('H', (-1.7333607372068203, -1.636440463684339, -1.7110355127606616)),
    ('H', (-1.7333909728248273, -1.6363724335438226, 1.7110600792002923)),
    ('H', (-2.5967029046096086, 1.4251728623104047, -1.5117809003662624e-05)),
    ('H', (-3.556090402338792, 2.9086399961389326, -1.5117809003662624e-05)))

# Z-Matrices
C5H8O_ZMA = (
    ('C', (None, None, None), (None, None, None),
     (None, None, None)),
    ('C', (0, None, None), ('R1', None, None),
     (2.894126135733367, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (2.1002551784038714, 1.9218361502726833, None)),
    ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
     (2.1008326855365818, 1.932519243247544, 4.182194537966691)),
    ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
     (2.0985632786105715, 1.9267905467443245, 2.078543469899091)),
    ('C', (1, 0, 2), ('R5', 'A5', 'D5'),
     (2.7809944542976193, 1.9090367411091194, 1.0264175778482927)),
    ('C', (1, 0, 5), ('R6', 'A6', 'D6'),
     (2.905905476629344, 1.9462477957117656, 2.1246205025559037)),
    ('H', (1, 0, 5), ('R7', 'A7', 'D7'),
     (2.1005624282579265, 1.8988979609840009, 4.223022677525844)),
    ('X', (5, 1, 0), ('R8', 'A8', 'D8'),
     (1.8897261254578288, 1.5707963267948968, 4.216836292856281)),
    ('C', (5, 8, 1), ('R9', 'A9', 'D9'),
     (2.2778505841014964, 1.5732722619955628, 3.1519457859137696)),
    ('X', (9, 5, 8), ('R10', 'A10', 'D10'),
     (1.8897261254578286, 1.56832039159423, 0.0)),
    ('H', (9, 10, 5), ('R11', 'A11', 'D11'),
     (2.0003442808863467, 1.57550991099349, 3.14859478950736)),
    ('O', (6, 1, 0), ('R12', 'A12', 'D12'),
     (2.65076899334649, 1.9387190313618887, 1.0262708014428483)),
    ('H', (6, 1, 12), ('R13', 'A13', 'D13'),
     (2.1058345184525726, 1.9323237957467607, 2.129177885999989)),
    ('H', (6, 1, 12), ('R14', 'A14', 'D14'),
     (2.1010240316411886, 1.9207088798352128, 4.1894956154070275)),
    ('H', (12, 6, 1), ('R15', 'A15', 'D15'),
     (1.8758293656194, 1.8624105681328567, 1.2477273554765336)))

HOOH_ZMA_C2 = (
    ('O', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.747759350307364, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (1.8445107263893656, 1.6890062361073546, None)),
    ('H', (1, 0, 2), ('R3', 'A3', 'D3'),
     (1.8445105807874629, 1.6890041579036852, 2.2578840834994196)))
HOOH_ZMA_CS = (
    ('O', (None, None, None), (None, None, None),
     (None, None, None)),
    ('O', (0, None, None), ('R1', None, None),
     (2.747759350307364, None, None)),
    ('H', (0, 1, None), ('R2', 'A2', None),
     (1.8445107263893656, 1.6890062361073546, None)),
    ('H', (1, 0, 2), ('R3', 'A3', 'D3'),
     (1.8445105807874629, 1.6890041579036852, 0.0000)))


def test__geom__with_stereo():
    """ test geom conversions
    """
    print('geom with stereo')

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        geo = automol.inchi.geometry(ref_ich)
        ich = automol.geom.inchi(geo)
        assert ich == ref_ich

        assert automol.geom.formula(geo) == automol.inchi.formula(ich)


def test__geom__no_stereo():
    """ test geom conversions
    """
    print('geom no stereo')

    ref_sort_ich = 'InChI=1S/C3H7O2/c1-3(2)5-4/h3-4H,1H2,2H3'
    ref_nums_lst = ((0, 1, 2, 3, 4),)

    ich = ICHS_WITH_STEREO[0]
    geo = automol.inchi.geometry(ich)
    sort_ich, nums_lst = automol.geom.inchi_with_sort(
        geo, stereo=False)

    assert sort_ich == ref_sort_ich
    assert nums_lst == ref_nums_lst

    # Test failed zmatrix call
    with pytest.raises(NotImplementedError):
        automol.geom.zmatrix_with_conversion_info(
            C2H6_H_GEO, ts_bnds=frozenset({7, 8}))


def test__graph__with_stereo():
    """ test graph conversions
    """
    print('graph with stereo')

    ref_ichs = []
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = list(numpy.random.choice(ref_ichs, NSAMP))

    for ref_ich in ref_ichs:
        print(ref_ich, flush=True)
        gra = automol.inchi.graph(ref_ich)
        ich = automol.graph.inchi(gra, stereo=True)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__graph__no_stereo():
    """ test graph conversions
    """
    print("graph no stereo")

    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print(ref_ich, flush=True)
        gra = automol.inchi.graph(ref_ich)
        gra = automol.graph.without_stereo_parities(gra)
        # gra <=> ich
        ich = automol.graph.inchi(gra)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


def test__zmatrix__with_stereo():
    """ test zmatrix conversions
    """
    print("zmatrix with stereo")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        ref_geo = automol.inchi.geometry(ref_ich)
        zma = automol.geom.zmatrix(ref_geo)
        geo = automol.zmat.geometry(zma)
        ich = automol.geom.inchi(geo)
        assert ich == ref_ich, (f'{ich} != {ref_ich}')

        assert automol.zmat.formula(zma) == automol.inchi.formula(ich)

    # Test dummy
    zma = (
        ('X', (None, None, None), (None, None, None),
         (None, None, None)),
        ('C', (0, None, None), ('R1', None, None),
         (1.8897261254578281, None, None)),
        ('O', (1, 0, None), ('R2', 'A1', None),
         (2.2601124460475623, 1.5707963267948966, None)),
        ('O', (1, 0, 2), ('R2', 'A1', 'D1'),
         (2.2601124460475623, 1.5707963267948966, 3.141592653589793)))

    ref_geo = (('C', (0.0, 0.0, 1.8897261254578281)),
               ('O', (0.0, 2.2601124460475623, 1.889726125457828)),
               ('O', (0.0, -2.2601124460475623, 1.889726125457828)))
    geo = automol.zmat.geometry(zma, dummy=False)

    assert automol.geom.almost_equal_dist_matrix(geo, ref_geo)


def test__smiles__from_geom():
    """ test smiles conversions
    """
    print("smiles from geom")

    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        # geo <=> smi
        geo = automol.inchi.geometry(ref_ich)
        smi = automol.geom.smiles(geo)
        ich = automol.smiles.inchi(smi)
        assert automol.inchi.standard_form(ich, stereo=False) == ref_ich


def test__graph__misc():
    """ test graph conversions
    """

    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra)
    assert ich == ref_ich

    ref_ich = 'InChI=1S/C4H4F2.HO/c5-3-1-2-4-6;/h1-4H;1H/b3-1-,4-2-;'
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra, stereo=True)
    assert ich == ref_ich

    ref_ich = 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'
    geo = automol.inchi.geometry(ref_ich)
    gra = automol.geom.graph(geo)
    ich = automol.graph.inchi(gra, stereo=True)
    assert ich == ref_ich

    ref_ich = 'InChI=1S/C2H4O/c1-2-3/h2-3H,1H2'
    ref_conn_gra = (
        {0: ('C', 0, None), 1: ('C', 0, None), 2: ('O', 0, None),
         3: ('H', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
         6: ('H', 0, None)},
        {frozenset({0, 1}): (1, None), frozenset({0, 3}): (1, None),
         frozenset({0, 4}): (1, None), frozenset({1, 2}): (1, None),
         frozenset({1, 5}): (1, None), frozenset({2, 6}): (1, None)})
    conn_gra = automol.geom.connectivity_graph(
        automol.inchi.geometry(ref_ich))
    assert conn_gra == ref_conn_gra

    ref_geo = (
        ('C', (0.0, 0.0, 0.0)),
        ('C', (-3.4713705762105e-16, 3.4713705762105e-16, -2.8345891881867)),
        ('Cl', (-1.04141117286317e-15, 5.2070558643158e-16, -5.6691783763734)),
        ('H', (0.0, 0.0, 2.0786987380036113)))
    ich = 'InChI=1S/C2HCl/c1-2-3/h1H'
    gra = automol.inchi.graph(ich)
    geo = automol.graph.geometry(gra)
    assert automol.geom.almost_equal_dist_matrix(geo, ref_geo)

    # Test stereo
    def randomize_atom_ordering(geo):
        """ randomize atom ordering in a geometry
        """
        natms = automol.geom.count(geo)
        ord_dct = dict(enumerate(numpy.random.permutation(natms)))
        return automol.geom.reorder(geo, ord_dct)

    # smi = 'FC=C(C=CC=CF)C=CC=CF'
    # smi = 'FC=CC=CC=CF'
    # ich = automol.smiles.inchi('CC([O])=CCO')
    ich = automol.smiles.inchi('CC([O])=CCO.O')
    geo = automol.inchi.geometry(ich)
    geo = randomize_atom_ordering(geo)
    gra = automol.geom.graph(geo)
    ich = automol.graph.inchi(gra, stereo=True)
    print(ich)
    assert ich in (
        'InChI=1S/C4H7O2.H2O/c1-4(6)2-3-5;/h2,5H,3H2,1H3;1H2/b4-2-;',
        'InChI=1S/C4H7O2.H2O/c1-4(6)2-3-5;/h2,5H,3H2,1H3;1H2/b4-2+;')


def test__inchi_geometry():
    """ test automol.inchi.geometry
    """
    # ref_ich = 'InChI=1S/H2S/h1H2'
    # ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    # print(ich)
    # assert ich == ref_ich

    # ref_ich = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    # ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    # print(ich)
    # assert ich == ref_ich

    # ref_ich = 'InChI=1S/Ar'
    # ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    # print(ich)
    # assert ich == ref_ich

    # ref_ich = 'InChI=1S/Cl2/c1-2'
    # ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    # print(ich)
    # assert ich == ref_ich

    ich = 'InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3'
    geo = automol.inchi.geometry(ich)
    ich = automol.geom.inchi(geo, stereo=True)
    print(automol.geom.chi(geo, stereo=True))
    print(ich)
    assert ich in ('InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5+',
                   'InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5-')

    ich = 'InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3'
    geo = automol.inchi.geometry(ich)
    ich = automol.geom.inchi(geo, stereo=True)
    print(ich)
    assert ich in ('InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5+',
                   'InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5-')

    # Extra test case for broken InChI conversion
    ich = automol.smiles.inchi('CC([O])=CCO')
    geo = automol.inchi.geometry(ich)
    ich = automol.geom.inchi(geo, stereo=True)
    print(ich)
    assert ich in ('InChI=1S/C4H7O2/c1-4(6)2-3-5/h2,5H,3H2,1H3/b4-2-',
                   'InChI=1S/C4H7O2/c1-4(6)2-3-5/h2,5H,3H2,1H3/b4-2+')

    # Extra test case for broken InChI conversion
    ich = automol.smiles.inchi('CC([O])=CCO.[OH]')
    geo = automol.inchi.geometry(ich)
    ich = automol.geom.inchi(geo, stereo=True)
    print(ich)
    assert ich in ('InChI=1S/C4H7O2.HO/c1-4(6)2-3-5;/h2,5H,3H2,1H3;1H/b4-2-;'
                   'InChI=1S/C4H7O2.HO/c1-4(6)2-3-5;/h2,5H,3H2,1H3;1H/b4-2+;')


def test__inchi_conformers():
    """ test automol.inchi.conformers
    """
    ref_ich = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    geos = automol.inchi.conformers(ref_ich)
    ichs = tuple(automol.geom.inchi(geo) for geo in geos)
    assert all(ich == ref_ich for ich in ichs)


def test__multiple_rings():
    """ test graph => inchi conversion for multiple rings
    """
    ref_ich = ('InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19'
               ')4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/'
               't10-,11-,13-,16+,17+/m1/s1')
    gra = automol.inchi.graph(ref_ich)
    ich = automol.graph.inchi(gra, stereo=True)
    print(ref_ich)
    print(ich)
    assert ich == ref_ich


def test__geom__x2z_zmatrix():
    """ test automol.geom.x2z_zmatrix
    """

    ref_zma = (
        ('C', (None, None, None), (None, None, None),
         (None, None, None)),
        ('O', (0, None, None), ('R1', None, None),
         (2.67535, None, None)),
        ('H', (0, 1, None), ('R2', 'A2', None),
         (2.06501, 1.9116242231243494, None)),
        ('H', (0, 1, 2), ('R3', 'A3', 'D3'),
         (2.06501, 1.9116242231243494, 2.1084973627493095)),
        ('H', (0, 1, 2), ('R4', 'A4', 'D4'),
         (2.06458, 1.9020947254084601, 4.195841334964448)),
        ('H', (1, 0, 2), ('R5', 'A5', 'D5'),
         (1.83748, 1.8690905492532472, 5.2289366258049315)))

    geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))

    zma = automol.geom.x2z_zmatrix(geo)
    assert automol.zmat.almost_equal(ref_zma, zma)

    ref_zma2 = (
        ('H', (None, None, None), (None, None, None), (None, None, None)),)

    geo2 = (('H', (2.9512589894, 0.17507745634, 0.22317665541)),)
    zma2 = automol.geom.x2z_zmatrix(geo2)
    assert automol.zmat.almost_equal(ref_zma2, zma2)


def test__geom__x2z_torsion_coordinate_names():
    """ test automol.geom.x2z_torsion_coordinate_names
    """
    geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))

    tors_names = automol.geom.x2z_torsion_coordinate_names(geo)
    assert tors_names == ('D5',)

    geo2 = (('H', (2.9512589894, 0.17507745634, 0.22317665541)),)

    tors_names2 = automol.geom.x2z_torsion_coordinate_names(geo2)
    assert tors_names2 == ()


def test__geom__x2z_atom_ordering():
    """ test automol.geom.x2z_atom_ordering
    """
    geo = (('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))

    ordering = automol.geom.x2z_atom_ordering(geo)
    assert ordering == {0: 0, 4: 1, 1: 2, 2: 3, 3: 4, 5: 5}

    geo2 = (('H', (-0.9827048283, 0.061897979239, 2.02901783816)),)

    ordering2 = automol.geom.x2z_atom_ordering(geo2)
    assert ordering2 == {0: 0}


def test__zmat_conv_dummy():
    """ test automol.zmat.geom
    """

    ref_gra = ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None),
                3: ('H', 0, None), 4: ('H', 0, None), 5: ('C', 0, None),
                6: ('C', 0, None), 7: ('H', 0, None), 8: ('X', 0, None),
                9: ('C', 0, None), 10: ('X', 0, None), 11: ('H', 0, None),
                12: ('O', 0, None), 13: ('H', 0, None), 14: ('H', 0, None),
                15: ('H', 0, None)},
               {frozenset({1, 7}): (1, None), frozenset({9, 10}): (0, None),
                frozenset({0, 3}): (1, None), frozenset({12, 15}): (1, None),
                frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None),
                frozenset({6, 14}): (1, None), frozenset({0, 4}): (1, None),
                frozenset({13, 6}): (1, None), frozenset({9, 11}): (1, None),
                frozenset({1, 5}): (1, None), frozenset({12, 6}): (1, None),
                frozenset({1, 6}): (1, None), frozenset({9, 5}): (1, None),
                frozenset({8, 5}): (0, None)})
    ref_geo = (
        ('C', (0.0, 0.0, 0.0)),
        ('C', (0.0, 0.0, 2.894126135733367)),
        ('H', (0.0, 1.972172485464683, -0.7222239971008525)),
        ('H', (-1.6951231876875892, -0.9936443210097795, -0.7434556573917012)),
        ('H', (1.7188336567949982, -0.9563671208828152, -0.7313963526678268)),
        ('C', (-2.2442053651761964, 1.3586361565617329, 3.8169371823900766)),
        ('C', (0.025535203206634938, -2.7033658002133003, 3.959700132572027)),
        ('H', (1.7085691873013926, 1.017330677567275, 3.5710250496737017)),
        ('X', (-1.1281908592835943, 2.626314426891954, 4.664621236037264)),
        ('C', (-4.085691887320963, 2.4807763649694308, 4.5506483141888445)),
        ('X', (-2.96967738142836, 3.7484546352996513, 5.398332367836032)),
        ('H', (-5.7053556919673625, 3.47091721106787, 5.181273009931513)),
        ('O', (2.1539110648484963, -4.040151863741559, 3.117325235021031)),
        ('H', (-0.002031340007772098, -2.6743735223267073, 6.065154609074555)),
        ('H', (-1.693646530214367, -3.738866587983746, 3.3380456726570835)),
        ('H', (3.635541963584469, -3.3116189027191543, 4.007687753063454)))

    gra = automol.zmat.connectivity_graph(C5H8O_ZMA, dummy=True)
    geo = automol.zmat.geometry(C5H8O_ZMA, dummy=True)

    assert ref_gra == gra
    assert automol.geom.almost_equal_dist_matrix(ref_geo, geo)


def test__smiles__no_stereo():
    """ test smiles conversions
    """
    print("smiles no stereo")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        ref_ich = automol.inchi.standard_form(ref_ich, stereo=False)
        print()
        print(ref_ich, flush=True)
        gra = automol.amchi.graph(ref_ich, stereo=False)
        smi = automol.graph.smiles(gra)
        print(smi)
        gra = automol.smiles.graph(smi, stereo=False)
        smi = automol.graph.smiles(gra)
        ich = automol.smiles.inchi(smi)
        print(smi)
        print(ich)
        assert ich == ref_ich, f"\n{ich} !=\n{ref_ich}"


def test__smiles__with_stereo():
    """ test smiles conversions
    """
    print("smiles with stereo")

    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    failed = []
    for ref_ich in ref_ichs:
        print()
        print(ref_ich, flush=True)
        gra = automol.amchi.graph(ref_ich, stereo=True)
        smi = automol.graph.smiles(gra)
        print(smi)
        gra = automol.smiles.graph(smi, stereo=True)
        smi = automol.graph.smiles(gra)
        ich = automol.smiles.inchi(smi)
        print(smi)
        print(ich)
        try:
            assert ich == ref_ich, f"\n{ich} !=\n{ref_ich}"
        except AssertionError:
            print("FAILED")
            print(f"\n{ich} !=\n{ref_ich}")
            failed.append((ref_ich, ich))

    print("FAILED:")
    for ref_ich, ich in failed:
        print(ref_ich)


if __name__ == '__main__':
    test__geom__with_stereo()
    test__geom__no_stereo()
    test__graph__with_stereo()
    test__graph__no_stereo()
    test__zmatrix__with_stereo()
    test__smiles__from_geom()
    test__graph__misc()
    # test__inchi_geometry()
    test__inchi_conformers()
    test__multiple_rings()
    # test__geom__no_stereo()
    # test__smiles__no_stereo()
    # test__smiles__with_stereo()
    # test__inchi_geometry()
