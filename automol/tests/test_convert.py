""" test automol.automol.convert
"""

import numpy
import automol
from automol.convert import _rdkit as rdkit
from automol.convert import _pybel as pybel
from _util import load_numpy_string_file
from _util import load_pandas_csv_string_file


BS_DF = load_pandas_csv_string_file(
    ['data'], 'badspecies.csv')
ICHS_NO_STEREO = load_numpy_string_file(
    ['data'], 'heptane_inchis_no_stereo.txt')
ICHS_WITH_STEREO = load_numpy_string_file(
    ['data'], 'heptane_inchis_with_stereo.txt')
NSAMP = 50


def test__geom__with_stereo():
    """ test geom conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        geo = automol.inchi.geometry(ref_ich)
        ich = automol.geom.inchi(geo)
        assert ich == ref_ich

        assert automol.geom.formula(geo) == automol.inchi.formula(ich)


def test__geom__no_stereo():
    """ test geom conversions
    """

    ref_sort_ich = 'InChI=1S/C3H7O2/c1-3(2)5-4/h3-4H,1H2,2H3'
    ref_nums = (0, 1, 2, 3, 4)

    ich = ICHS_WITH_STEREO[0]
    geo = automol.inchi.geometry(ich)
    sort_ich, nums = automol.convert.geom.inchi_with_sort(
        geo, stereo=False)

    assert sort_ich == ref_sort_ich
    assert nums == ref_nums


def test__graph__with_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_WITH_STEREO
    if NSAMP is not None:
        ref_ichs = list(numpy.random.choice(ref_ichs, NSAMP))

    # AVC note to self -- fix this case:
    ref_ichs.append(
        'InChI=1S/C5H10O3/c1-4-2-5(8-4)3-7-6/h4-6H,2-3H2,1H3/t4-,5-/m1/s1')

    for ref_ich in ref_ichs:
        print(ref_ich)
        gra = automol.inchi.graph(ref_ich)
        ich = automol.graph.inchi(gra)
        print(ich)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)
        print()


def test__graph__no_stereo():
    """ test graph conversions
    """
    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
        gra = automol.inchi.graph(ref_ich)
        gra = automol.graph.without_stereo_parities(gra)
        # gra <=> ich
        ich = automol.graph.inchi(gra)
        assert ich == ref_ich

        assert automol.graph.formula(gra) == automol.inchi.formula(ich)


# def test__zmatrix__with_stereo():
#     """ test zmatrix conversions
#     """
#     ref_ichs = ICHS_WITH_STEREO
#     if NSAMP is not None:
#         ref_ichs = numpy.random.choice(ref_ichs, NSAMP)
#
#     for ref_ich in ref_ichs:
#         ref_geo = automol.inchi.geometry(ref_ich)
#         zma = automol.geom.zmatrix(ref_geo)
#         geo = automol.zmat.geometry(zma)
#         ich = automol.geom.inchi(geo)
#         assert ich == ref_ich
#
#         assert automol.zmat.formula(zma) == automol.inchi.formula(ich)
#
#     # Test dummy
#     zma = (
#         ('X', (None, None, None), (None, None, None),
#          (None, None, None)),
#         ('C', (0, None, None), ('R1', None, None),
#          (1.8897261254578281, None, None)),
#         ('O', (1, 0, None), ('R2', 'A1', None),
#          (2.2601124460475623, 1.5707963267948966, None)),
#         ('O', (1, 0, 2), ('R2', 'A1', 'D1'),
#          (2.2601124460475623, 1.5707963267948966, 3.141592653589793)))
#
#     ref_geo = (('C', (0.0, 0.0, 1.8897261254578281)),
#                ('O', (0.0, 2.2601124460475623, 1.889726125457828)),
#                ('O', (0.0, -2.2601124460475623, 1.889726125457828)))
#     geo = automol.zmat.geometry(zma, dummy=False)
#
#     assert automol.geom.almost_equal_dist_matrix(geo, ref_geo)


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


def test__smiles__from_geom():
    """ test smiles conversions
    """
    ref_ichs = ICHS_NO_STEREO
    if NSAMP is not None:
        ref_ichs = numpy.random.choice(ref_ichs, NSAMP)

    for ref_ich in ref_ichs:
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
    ich = automol.graph.inchi(gra)
    assert ich == ref_ich

    ref_ich = 'InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+'
    geo = automol.inchi.geometry(ref_ich)
    gra = automol.geom.graph(geo)
    ich = automol.graph.inchi(gra)
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


def test__inchi_geometry():
    """ test automol.inchi.geometry
    """
    ref_ich = 'InChI=1S/H2S/h1H2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich

    ref_ich = 'InChI=1S/Cl2/c1-2'
    ich = automol.geom.inchi(automol.inchi.geometry(ref_ich))
    assert ich == ref_ich

    ref_ich = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    pbm = pybel.from_inchi(ref_ich)
    ich = pybel.to_inchi(pbm)
    assert ich == ref_ich
    rdm = rdkit.from_inchi(ref_ich)
    ich = automol.geom.inchi(rdkit.to_geometry(rdm))
    assert ich == ref_ich

    ref_ich = 'InChI=1S/Ar'
    rdm = rdkit.from_inchi(ref_ich)
    ich = automol.geom.inchi(rdkit.to_geometry(rdm))
    assert ich == ref_ich


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
    assert set(tors_names) <= set(automol.zmat.torsion_coordinate_names(zma))

    geo2 = (('H', (2.9512589894, 0.17507745634, 0.22317665541)),)

    tors_names2 = automol.geom.zmatrix_torsion_coordinate_names(geo2)
    assert tors_names2 == ()


def test__geom__zmatrix_atom_ordering():
    """ test automol.geom.zmatrix_atom_ordering
    """
    geo = (('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))

    ordering = automol.geom.zmatrix_atom_ordering(geo)
    assert ordering == {0: 0, 4: 1, 1: 2, 2: 3, 3: 4, 5: 5}

    geo2 = (('H', (-0.9827048283, 0.061897979239, 2.02901783816)),)

    ordering2 = automol.geom.zmatrix_atom_ordering(geo2)
    assert ordering2 == {0: 0}


# def test__sarah_badpsecies():
#     """ test bad species list from sarah
#     """
#     bs_ichs = list(BS_DF['inchi'])
#
#     for ref_ich in bs_ichs:
#         print(ref_ich)
#         ich = automol.inchi.add_stereo(ref_ich)
#         print(ich)
#         print()


if __name__ == '__main__':
    test__graph__with_stereo()
