""" test automol.geom
"""

# import pytest
import numpy
import automol
from automol import geom


C2H2CLF_GEO = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
C2H6_GEO = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('H', (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)))
C2H6_GEO_2 = (
    ('C', (-1.4283035320563338, 0.013425343735546437, -0.030302158896694683)),
    ('C', (1.4283027358735494, -0.013425597530894248, 0.0303022919384165)),
    ('H', (-2.1448124562913287, -1.5396513482615042, -1.191852168914227)),
    ('H', (2.1448121742707795, 1.539654946791746, 1.1918517388178247)),
    ('H', (2.1972712765396953, 0.1922944277301287, -1.8778395029874426)),
    ('H', (2.121312248031497, -1.7927029137609576, 0.8231123911174519)),
    ('H', (-2.1972722614281355, -0.19229727219177065, 1.8778380427620682)),
    ('H', (-2.121310184939721, 1.792702413487708, -0.8231106338374065)))
C2H5_GEO = (
    ('C', (1.609011843391662, -0.13776060061651263, -0.04009265860394636)),
    ('C', (-1.1834240897842176, 0.11044868962990169, -0.032202756665620995)),
    ('H', (2.789386604699099, 1.4908722862660362, 0.3439484428110317)),
    ('H', (2.4787984007354567, -1.9869974801329033, -0.1705361344367882)),
    ('H', (-1.7411967628258649, 1.9663729732966975, -0.7497840216103133)),
    ('H', (-2.036027793494113, -1.334976675748384, -1.2381516164697437)),
    ('H', (-1.9165482027220222, -0.10795919269481989, 1.8868187449753886)))
H2O_GEO = (('H', (-1.0375307411, 1.1972456627, 0.1366217261)),
           ('O', (-0.3546124207, 0.6455743406, 0.5947534694)),
           ('H', (-0.8777308204, -0.0012671711, 1.1319188373)))
H_GEO = (('H', (-0.9827048283, 0.061897979239, 2.02901783816)),)
H2_GEO = (('H', (0.6595322581760189, 0.0, 0.0)),
          ('H', (-0.6595322581760189, 0.0, 0.0)))
HCCH_GEO = (('C', (-1.13372239064879, 0.0082038553557789, 0.268047200455629)),
            ('C', (1.133722408498827, -0.0082033046582439, 0.210166207903664)),
            ('H', (-3.14709200757641, 0.0227720059492319, 0.319442062072525)),
            ('H', (3.147091989726373, -0.022772556646775, 0.158770666805506)))
BAD_GEO = ((2.994881276150, -1.414434615111),
           (1.170155936996, 0.359360756989),
           (-1.201356763194, -0.347546894407))

GEO_STR = """4
Energy:  1.25828 kJ/mol (UFF)
O       -0.9741512421      0.1822933377     -0.1126377411
O        0.3289735103      0.1763176636     -0.3024822120
H       -1.2826135304     -0.6998243267     -0.4438019186
H        0.4659939252      0.6943122175     -1.1366628384
"""

HABS_ZRXN_STR1 = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  17: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 1, stereo_parity: null}
  2-5: {order: 1, stereo_parity: null}
  2-6: {order: 1, stereo_parity: null}
  4-7: {order: 1, stereo_parity: null}
  4-8: {order: 1, stereo_parity: null}
  4-9: {order: 1, stereo_parity: null}
  7-10: {order: 1, stereo_parity: null}
  7-11: {order: 1, stereo_parity: null}
  10-12: {order: 0.9, stereo_parity: null}
  10-13: {order: 1, stereo_parity: null}
  10-14: {order: 1, stereo_parity: null}
  12-15: {order: 0, stereo_parity: null}
  12-16: {order: 0.1, stereo_parity: null}
  16-17: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
- [16, 17]
backward TS atoms:
  1: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  9: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  10: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  11: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  12: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  13: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  14: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  15: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  16: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 0.9, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  2-4: {order: 0.1, stereo_parity: null}
  4-7: {order: 1, stereo_parity: null}
  4-10: {order: 1, stereo_parity: null}
  4-11: {order: 1, stereo_parity: null}
  5-6: {order: 1, stereo_parity: null}
  5-7: {order: 1, stereo_parity: null}
  5-12: {order: 1, stereo_parity: null}
  5-13: {order: 1, stereo_parity: null}
  6-8: {order: 1, stereo_parity: null}
  6-14: {order: 1, stereo_parity: null}
  6-15: {order: 1, stereo_parity: null}
  7-9: {order: 1, stereo_parity: null}
  8-16: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3]
- [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
"""
HABS_GEO1 = (
    ('O', (0.0, 0.0, 0.0)),
    ('C', (0.0, 0.0, 2.6914294153467657)),
    ('H', (0.0, 1.7534707783960422, -0.5517555484489822)),
    ('C', (-2.623299394112204, 0.6566385034518767, 3.6480149247207985)),
    ('H', (0.5636237182797816, -1.88605768213913, 3.3289129545424223)),
    ('H', (1.4260767249009725, 1.345507128950089, 3.35470487715432)),
    ('C', (-2.746065179408126, 0.6724245480726039, 6.5028326846462985)),
    ('H', (-3.1961339563014963, 2.518941659280532, 2.9483507693607525)),
    ('H', (-3.988069056064281, -0.7224120831906182, 2.9241998624297434)),
    ('C', (-5.340762043167491, 0.8914938898502375, 7.652631925159079)),
    ('O', (-0.845554426010156, 0.5689858109064159, 7.8373344196418655)),
    ('H', (-6.502001079557306, -0.7071271436769692, 7.052333079032717)),
    ('H', (-6.21125009870658, 2.6681404185712547, 7.061161472892827)),
    ('H', (-5.194262776062822, 0.8701632805783928, 9.713822470820245)),
    ('O', (-8.5007619646225, -2.928405157276118, 6.591025292968647)),
    ('H', (-9.48421314672397, -1.6400034467104858, 5.7038001145820685))
)


HABS_ZRXN_STR2 = """
reaction class: hydrogen abstraction
forward TS atoms:
  1: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: X, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  8: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
forward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 1, stereo_parity: null}
  1-4: {order: 1, stereo_parity: null}
  1-5: {order: 0.9, stereo_parity: null}
  5-6: {order: 0, stereo_parity: null}
  5-7: {order: 0.1, stereo_parity: null}
  7-8: {order: 1, stereo_parity: null}
reactants keys:
- [1, 2, 3, 4, 5, 6]
- [7, 8]
backward TS atoms:
  1: {symbol: O, implicit_hydrogen_valence: 0, stereo_parity: null}
  2: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  3: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  4: {symbol: C, implicit_hydrogen_valence: 0, stereo_parity: null}
  5: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  6: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
  7: {symbol: H, implicit_hydrogen_valence: 0, stereo_parity: null}
backward TS bonds:
  1-2: {order: 1, stereo_parity: null}
  1-3: {order: 0.9, stereo_parity: null}
  3-4: {order: 0.1, stereo_parity: null}
  4-5: {order: 1, stereo_parity: null}
  4-6: {order: 1, stereo_parity: null}
  4-7: {order: 1, stereo_parity: null}
products keys:
- [1, 2, 3]
- [4, 5, 6, 7]
"""
HABS_GEO2 = (
    ('C', (0.0, 0.0, 0.0)),
    ('H', (0.0, 0.0, 2.0639575202459945)),
    ('H', (0.0, 1.9459185282046763, -0.6879842041977405)),
    ('H', (-1.6852142339096723, -0.9729596887682943, -0.687986206970089)),
    ('H', (1.6852153438276314, -0.9729582966864149, -0.687984673708243)),
    ('O', (3.875541447997114, -2.4707898221289817, -2.1374322380005975)),
    ('H', (4.3062792658506215, -0.8227905834100768, -2.8535934936888445))
)


def test__from_data():
    """ test getters
    """
    assert C2H2CLF_GEO == geom.from_data(
        symbs=geom.symbols(C2H2CLF_GEO),
        xyzs=geom.coordinates(C2H2CLF_GEO),
    )


def test__is_valid():
    """ test geom.is_valid
    """

    # Check validity
    assert geom.is_valid(C2H2CLF_GEO)
    # with pytest.raises(ValueError):
    #     geom.is_valid(BAD_GEO)

    # Check empty geom returns empty information
    assert not geom.symbols(())
    assert not geom.coordinates(())


def test__struct_check():
    """ test geom.is_atom
        test geom.is_linear
    """

    assert geom.is_atom(H_GEO)
    assert not geom.is_atom(C2H2CLF_GEO)

    assert not geom.is_linear(H_GEO)
    assert geom.is_linear(H2_GEO)
    assert not geom.is_linear(C2H2CLF_GEO)
    assert geom.is_linear(HCCH_GEO)


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

    assert geom.xyz_string_comment(GEO_STR) == 'Energy:  1.25828 kJ/mol (UFF)'


def test__formula():
    """ test geom.formula
    """
    assert geom.formula(C2H2CLF_GEO) == {'F': 1, 'C': 2, 'Cl': 1, 'H': 2}


def test__atom_indices():
    """ test geom.atom_indices
        test geom.atom_indices
    """

    hidxs = geom.atom_indices(C2H2CLF_GEO, 'H', match=True)
    heavyidxs = geom.atom_indices(C2H2CLF_GEO, 'H', match=False)

    assert hidxs == (4, 5)
    assert heavyidxs == (0, 1, 2, 3)

    h_count = geom.atom_count(C2H2CLF_GEO, 'H', match=True)
    heavy_count = geom.atom_count(C2H2CLF_GEO, 'H', match=False)

    assert h_count == 2
    assert heavy_count == 4


def test__dist_analysis():
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

    ref_dist_mat = (
        (0.0000000, 2.5616993, 4.3547794, 6.6877445, 3.9644122, 4.7627773),
        (2.5616993, 0.0000000, 2.4806333, 4.3481308, 2.0452679, 3.8991098),
        (4.3547794, 2.4806333, 0.0000000, 2.5392899, 3.9684476, 2.0228115),
        (6.6877445, 4.3481308, 2.5392899, 0.0000000, 4.8648485, 3.9664779),
        (3.9644122, 2.0452679, 3.9684476, 4.8648485, 0.0000000, 5.7501111),
        (4.7627773, 3.8991098, 2.0228115, 3.9664779, 5.7501111, 0.0000000))
    assert numpy.allclose(geom.distance_matrix(geo), ref_dist_mat)


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


def test__mass():
    """ test geom.masses
        test geom.center_of_mass()
        test geom.mass_centered()
    """

    ref_masses = (1.00782503223, 15.99491461957, 1.00782503223)
    assert numpy.allclose(geom.masses(H2O_GEO, amu=True), ref_masses)

    ref_masses2 = (1837.1526464817923, 29156.94568388855, 1837.1526464817923)
    assert numpy.allclose(geom.masses(H2O_GEO, amu=False), ref_masses2)

    ref_red_mass = 0.9544182343494726
    assert numpy.isclose(geom.reduced_mass(H2O_GEO, H_GEO), ref_red_mass)

    # make sure the COM for an uncentered geometry is non-zero
    cm_xyz = automol.geom.center_of_mass(C2H2CLF_GEO)
    assert not numpy.allclose(cm_xyz, 0.)

    # now make sure centering it yields a COM at the origin
    geo = automol.geom.mass_centered(C2H2CLF_GEO)
    cm_xyz = automol.geom.center_of_mass(geo)
    assert numpy.allclose(cm_xyz, 0.)


def test__rotation_properties():
    """ test geom.principal_axes()
        test geom.rotational_constants()
    """

    ref_axes = (
        (-0.10203856485038033, 0.99478044375795, 1.0460242302053948e-07),
        (0.765298080540252, 0.07849958398718046, 0.6388713980413003),
        (-0.6355367646365401, -0.06518960063212392, 0.7693135490583429))
    axes = automol.geom.principal_axes(H2O_GEO)
    assert numpy.allclose(axes, ref_axes)

    ref_cons = (
        2.9448201404714373e-06, 1.566795638659629e-07, 1.4876452660293155e-07)
    cons = geom.rotational_constants(C2H2CLF_GEO)
    assert numpy.allclose(cons, ref_cons)


def test__swap_coordinates():
    """ test geom.swap_coordinates
    """
    ref_geo = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    swp_geo = automol.geom.swap_coordinates(C2H2CLF_GEO, 1, 4)
    assert swp_geo == ref_geo


def test__reflect_coordinates():
    """ test geom.reflect_coordinates
    """
    ref_geo1 = (('F', (2.99488127615, -1.414434615111, -0.807144415388)),
                ('C', (-1.170155936996, 0.359360756989, -0.513323178859)),
                ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
                ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
                ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
                ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    ref_geo2 = (('F', (2.99488127615, -1.414434615111, -0.807144415388)),
                ('C', (-1.170155936996, -0.359360756989, 0.513323178859)),
                ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
                ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
                ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
                ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    ref_geo3 = (('F', (2.99488127615, -1.414434615111, -0.807144415388)),
                ('C', (-1.170155936996, -0.359360756989, -0.513323178859)),
                ('C', (1.201356763194, 0.347546894407, -0.3408392500119)),
                ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
                ('H', (-1.731596406235, -2.324260256203, -0.4292070203467)),
                ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    geo1 = automol.geom.reflect_coordinates(C2H2CLF_GEO, [1], ['x'])
    geo2 = automol.geom.reflect_coordinates(C2H2CLF_GEO, [1], ['x', 'y', 'z'])
    geo3 = automol.geom.reflect_coordinates(C2H2CLF_GEO, [1, 2, 4], ['x', 'y'])

    assert automol.geom.almost_equal_dist_matrix(geo1, ref_geo1, thresh=0.001)
    assert automol.geom.almost_equal_dist_matrix(geo2, ref_geo2, thresh=0.001)
    assert automol.geom.almost_equal_dist_matrix(geo3, ref_geo3, thresh=0.001)


def test__remove():
    """ test geom.remove
    """
    ref_geo1 = (('C', (1.170155936996, 0.359360756989, -0.513323178859)),
                ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
                ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    geo1 = automol.geom.remove(C2H2CLF_GEO, [0, 2, 4])

    assert automol.geom.almost_equal_dist_matrix(geo1, ref_geo1, thresh=0.001)


def test__translate():
    """ test geom.translate
    """
    ref_geo1 = (('F', (12.994881276150, -1.414434615111, -0.807144415388)),
                ('C', (11.170155936996, 0.359360756989, -0.513323178859)),
                ('C', (8.79864323681, -0.347546894407, -0.3408392500119)),
                ('Cl', (6.97202912502, 1.39211904938, -0.0492290974807)),
                ('H', (11.731596406235, 2.324260256203, -0.4292070203467)),
                ('H', (8.33269401879, -2.31375855306, -0.433949091252)))
    geo1 = automol.geom.translate(C2H2CLF_GEO, (10.0, 0.0, 0.0))

    assert automol.geom.almost_equal_dist_matrix(geo1, ref_geo1, thresh=0.001)


def test__rotate():
    """ test geom.euler_rotate
    """
    ref_geo1 = (('F', (-3.0558758007514, -1.1509331162337, -0.9789776801786)),
                ('C', (-1.0681262247754, -0.5805345601543, 0.532908762282)),
                ('C', (1.1776424945939, -0.26506900172108, -0.4723503888912)),
                ('Cl', (3.163419853813, 0.30082019331320, 1.0056384097658)),
                ('H', (-1.3983240017954, -0.4010882123173, 2.5433538391085)),
                ('H', (1.4140728611330, -0.46140633610660, -2.471679960480)))

    theta, phi, psi = 15.0/numpy.pi, 30.0/numpy.pi, 20.0/numpy.pi
    geo1 = automol.geom.euler_rotate(C2H2CLF_GEO, theta, phi, psi)

    assert automol.geom.almost_equal_dist_matrix(geo1, ref_geo1, thresh=0.001)


def test__closest_unbonded_atoms():
    """ test geom.closest_unbonded_atoms
    """

    ref_bnd_key = frozenset({1, 5})
    ref_dist_val = 3.8991097995323956
    bnd_key, dist_val = automol.geom.closest_unbonded_atoms(C2H2CLF_GEO)

    assert bnd_key == ref_bnd_key
    assert numpy.isclose(dist_val, ref_dist_val)


def test__permutations():
    """ test geom.rot_permutated_geoms
        test.geom.permutation
    """

    ref_perm_geo = (
        ('C', (1.609011843391662, -0.13776060061651263, -0.04009265860394636)),
        ('C', (-1.183424089784217, 0.1104486896299016, -0.03220275666562099)),
        ('H', (2.789386604699099, 1.4908722862660362, 0.3439484428110317)),
        ('H', (2.4787984007354567, -1.9869974801329033, -0.1705361344367882)),
        ('H', (-1.9165482027220222, -0.10795919269481989, 1.8868187449753886)),
        ('H', (-1.7411967628258649, 1.9663729732966975, -0.7497840216103133)),
        ('H', (-2.036027793494113, -1.334976675748384, -1.2381516164697437)))
    perm_geos = automol.geom.rot_permutated_geoms(C2H5_GEO)
    assert any(automol.geom.almost_equal_dist_matrix(geo, ref_perm_geo)
               for geo in perm_geos)

    ref_perm_idxs = (0, 1, 4, 5, 6, 7, 2, 3)
    perm_idxs = automol.geom.permutation(C2H6_GEO_2, C2H6_GEO)
    assert perm_idxs == ref_perm_idxs


def test__traj():
    """ test geom.from_xyz_trajectory_string
    """

    ref_traj_str = """ 3
comment 1
C    0.000000   0.000000   0.000000
C    0.000000   0.000000   1.000000
C    0.000000   0.000000   2.000000
 3
comment 2
C    0.000000   0.000000   3.000000
C    0.000000   0.000000   4.000000
C    0.000000   0.000000   5.000000
 3
comment 3
C    0.000000   0.000000   6.000000
C    0.000000   0.000000   7.000000
C    0.000000   0.000000   8.000000"""

    traj = automol.geom.from_xyz_trajectory_string(ref_traj_str)
    geoms = tuple(geo for geo, _ in traj)
    comments = tuple(comment for _, comment in traj)
    traj_str = automol.geom.xyz_trajectory_string(geoms, comments=comments)
    assert ref_traj_str == traj_str


def test__insert_dummies():
    """ test geom.insert_dummies
    """
    # 1. Generate a z-matrix to start with
    ich = automol.smiles.inchi('CC#CC#CCCCC#CC')
    zma = automol.geom.zmatrix(automol.inchi.geometry(ich))

    # 2. Convert to cartesians and remove dummy atoms
    geo = automol.zmat.geometry(zma)

    # Assume this geometry was changed somehow, such as by optimization, and we
    # want to update the z-matrix from it. The problem is that the geometry
    # doesn't contain the dummy atoms in the z-matrix.

    # 3. Determine dummy atom keys from the original z-matrix
    dummy_key_dct = automol.zmat.dummy_key_dictionary(zma)

    # 4. Insert dummy atoms to the new geometry at the appropriate positions
    geo_wdummy = automol.geom.insert_dummies(geo, dummy_key_dct)
    assert automol.geom.symbols(geo_wdummy) == automol.zmat.symbols(zma)

    # 5. Update the z-matrix from the geometry.
    zma_new = automol.zmat.from_geometry(zma, geo_wdummy)

    # 6. Check that the geometry was correctly encoded by converting the new
    # z-matrix back and comparing.
    geo_new = automol.zmat.geometry(zma_new)
    assert automol.geom.almost_equal_dist_matrix(geo, geo_new)


def test__hydrogen_bonded_structure():
    """ test automol.geom.hydrogen_bonded_strcure
    """

    ich = automol.smiles.inchi('CCCC[O]')
    geo = automol.inchi.geometry(ich)
    print(automol.geom.hydrogen_bonded_structure(geo, grxn=None))

    ich = automol.smiles.inchi('CC(=O)CO')
    geo = automol.inchi.geometry(ich)
    print(automol.geom.hydrogen_bonded_structure(geo, grxn=None))

    # H-abstraction with OH. (1) Has internal H-Bond, (2) No internal H-Bond
    zrxn = automol.reac.from_string(HABS_ZRXN_STR1)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.hydrogen_bonded_structure(HABS_GEO1, grxn=grxn))

    zrxn = automol.reac.from_string(HABS_ZRXN_STR2)
    grxn = automol.reac.relabel_for_geometry(zrxn)
    print(automol.geom.hydrogen_bonded_structure(HABS_GEO2, grxn=grxn))


def test__change_zmatrix_row_values():
    """ test automol.geom.change_zmatrix_row_values
    """
    ref_geo = (('O', (0.0, 0.0, 0.0)),
               ('C', (0.0, 0.0, 2.6100437640951895)),
               ('H', (0.0, 1.795784435, -0.9004532847666584)),
               ('N', (2.149523799, 1.6103368939, 3.8411840715311034)),
               ('H', (3.101869570, 0.04614797866, 4.810112185573969)),
               ('H', (1.052147147, 2.2771937846, 5.430572110076959)))
    ref_zma = geom.zmatrix(ref_geo)
    ref_vals = geom.zmatrix_row_values(ref_geo, 3, idx1=1, idx2=0, idx3=2)

    geo = geom.change_zmatrix_row_values(ref_geo, 3,
                                         dist=1.0, idx1=1,
                                         ang=120., idx2=0,
                                         dih=0., idx3=2)

    # First, check that the change resulted in these values
    vals = geom.zmatrix_row_values(geo, 3, idx1=1, idx2=0, idx3=2)
    assert numpy.allclose(vals, (1., 120., 0.))

    # Now, check that the other z-matrix values weren't affected. To do so,
    # change back to the original values and compare z-matrices.
    geo = geom.change_zmatrix_row_values(geo, 3,
                                         dist=ref_vals[0], idx1=1,
                                         ang=ref_vals[1], idx2=0,
                                         dih=ref_vals[2], idx3=2)
    zma = geom.zmatrix(geo)
    assert automol.zmat.almost_equal(zma, ref_zma)


if __name__ == '__main__':
    # test__hydrogen_bonded_structure()
    test__change_zmatrix_row_values()
