""" test the automol.geom module
"""
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


if __name__ == '__main__':
    # test__connectivity_graph()
    # test__inchi()
    # test__stereo_inchi()
    # test__from_xyz_string()
    # test__from_data()
    # test__is_valid()
    test__from_string()
