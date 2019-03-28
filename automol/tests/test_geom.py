""" test the automol.geom module
"""
import numpy
import automol
from automol import geom

C2H2CLF_ICH = 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H'
C2H2CLF_STE_ICH = 'InChI=1S/C2H2ClF/c3-1-2-4/h1-2H/b2-1+'
C2H2CLF_GEO = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
C2H2CLF_CGR = ({0: ('F', 0, None), 1: ('C', 0, None), 2: ('C', 0, None),
                3: ('Cl', 0, None), 4: ('H', 0, None), 5: ('H', 0, None)},
               {frozenset({0, 1}): (1, None), frozenset({1, 2}): (1, None),
                frozenset({1, 4}): (1, None), frozenset({2, 3}): (1, None),
                frozenset({2, 5}): (1, None)})
C2H2CLF_XYZ_STR = """6
charge: 0, mult: 1
F    1.584823  -0.748487  -0.427122
C    0.619220   0.190166  -0.271639
C   -0.635731  -0.183914  -0.180364
Cl  -1.602333   0.736678  -0.026051
H    0.916321   1.229946  -0.227127
H   -0.882300  -1.224388  -0.229636
"""


def test__is_valid():
    """ test geom.is_valid
    """
    assert geom.is_valid(C2H2CLF_GEO)


def test__from_data():
    """ test geom.from_data
    """
    assert C2H2CLF_GEO == geom.from_data(
        symbols=geom.symbols(C2H2CLF_GEO),
        coordinates=geom.coordinates(C2H2CLF_GEO),
    )


def test__zmatrix():
    """ test geom.zmatrix
    """
    geo = (('C', (-0.70116587131, 0.0146227007587, -0.016166607003)),
           ('O', (1.7323365056, -0.9538524899, -0.5617192010)),
           ('H', (-0.9827048283, 0.061897979239, 2.02901783816)),
           ('H', (-0.8787925682, 1.91673409124, -0.80019507919)),
           ('H', (-2.12093033745, -1.21447973767, -0.87411360631)),
           ('H', (2.9512589894, 0.17507745634, 0.22317665541)))
    zma = geom.zmatrix(geo)

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

    tors_names = geom.zmatrix_torsion_coordinate_names(geo)
    assert set(tors_names) <= set(automol.zmatrix.dihedral_names(zma))


def test__set_coordinates():
    """ test geom.set_coordinates
    """
    ref_geo = (('F', (0., 0., 0.)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (1., 1., 1.)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    geo = geom.set_coordinates(C2H2CLF_GEO, {0: [0., 0., 0.],
                                             3: [1., 1., 1.]})
    assert geom.almost_equal(geo, ref_geo)


def test__connectivity_graph():
    """ test geom.connectivity_graph
    """
    assert geom.connectivity_graph(C2H2CLF_GEO) == C2H2CLF_CGR


def test__inchi():
    """ test geom.inchi
    """
    assert geom.inchi(C2H2CLF_GEO) == C2H2CLF_ICH


def test__stereo_inchi():
    """ test geom.inchi
    """
    assert geom.stereo_inchi(C2H2CLF_GEO) == C2H2CLF_STE_ICH


def test__from_string():
    """ test geom.from_string
    """
    assert geom.almost_equal(
        geom.from_string(geom.string(C2H2CLF_GEO)), C2H2CLF_GEO)


def test__from_xyz_string():
    """ test geom.from_xyz_string
    """
    assert geom.almost_equal(
        geom.from_xyz_string(geom.xyz_string(C2H2CLF_GEO)), C2H2CLF_GEO)


def test__formula():
    """ test geom.formula
    """
    assert geom.formula(C2H2CLF_GEO) == {'F': 1, 'C': 2, 'Cl': 1, 'H': 2}


def test__coulomb_spectrum():
    """ test geom.coulomb_spectrum
    """
    ref_coul_spec = (
        0.23373850982000086, 0.26771181015927226, 21.472418990888897,
        38.92412503488664, 104.1418603738336, 456.00384141400343)

    assert numpy.allclose(geom.coulomb_spectrum(C2H2CLF_GEO), ref_coul_spec)

    for _ in range(10):
        axis = numpy.random.rand(3)
        angle = numpy.random.rand()
        geo = geom.rotate(C2H2CLF_GEO, axis, angle)
        assert numpy.allclose(geom.coulomb_spectrum(geo), ref_coul_spec)


def test__argunique_coulomb_spectrum():
    """ test geom.argunique_coulomb_spectrum
    """
    ref_idxs = (0, 3, 5, 8)

    geo = C2H2CLF_GEO
    natms = len(geom.symbols(geo))

    geos = []
    for idx in range(10):
        axis = numpy.random.rand(3)
        angle = numpy.random.rand()
        geo = geom.rotate(geo, axis, angle)

        if idx in ref_idxs and idx != 0:
            idx_to_change = numpy.random.randint(0, natms)
            new_xyz = numpy.random.rand(3)
            geo = geom.set_coordinates(geo, {idx_to_change: new_xyz})

        geos.append(geo)

    idxs = geom.argunique_coulomb_spectrum(geos)
    assert idxs == ref_idxs


if __name__ == '__main__':
    # test__connectivity_graph()
    # test__inchi()
    # test__stereo_inchi()
    # test__from_xyz_string()
    # test__from_data()
    # test__is_valid()
    # test__formula()
    # test__coulomb_spectrum()
    # test__set_coordinates()
    # test__argunique_coulomb_spectrum()
    test__from_string()
    test__from_xyz_string()
    test__zmatrix()
