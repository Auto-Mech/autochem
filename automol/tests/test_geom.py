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
    assert geom.is_valid(C2H2CLF_GEO)
    # with pytest.raises(ValueError):
    #     geom.is_valid(BAD_GEO)


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


def test__remove_coordinates():
    """ test geom.remove_coordinates
    """
    ref_geo1 = (('C', (1.170155936996, 0.359360756989, -0.513323178859)),
                ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
                ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    geo1 = automol.geom.remove_coordinates(C2H2CLF_GEO, [0, 2, 4])

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


def test__symmetry_factor():
    """ test geom.external_symmety_factor
        test geom.end_group_symmetry_factor
    """
    ref_sym_num1 = 12
    ref_sym_num2 = 1
    ref_sym_num3 = 0.5
    ref_sym_num4 = 1

    methane_geo = (('C', (1.2069668249, 1.9997649792, -0.0000004209)),
                   ('H', (3.3034303116, 1.9997688296, -0.0000006619)),
                   ('H', (0.5081497935, 0.1480118251, 0.6912478445)),
                   ('H', (0.5081445849, 2.3270011544, -1.9492882688)),
                   ('H', (0.5081426638, 3.5242779889, 1.2580393046)))
    c2h5of_geo = (('C', (-4.67963119210, -2.785693400767, -0.04102938592633)),
                  ('C', (-1.806009533535, -2.594940600449, -0.1025157659970)),
                  ('H', (-5.39544527869, -3.740953123044, -1.774996188159)),
                  ('H', (-5.501156952723, -0.854962636480, -0.01317371318990)),
                  ('O', (-5.48010155525, -4.07121874876, 2.132641999777597)),
                  ('H', (-1.208455201406, -1.52313066520, -1.804561025201)),
                  ('F', (-0.745999108314, -4.9827454242, -0.1878162481225)),
                  ('H', (-1.11738479998, -1.591763680324, 1.607521773704)),
                  ('H', (-5.30771129777, -5.90407965309, 1.771303996279)))
    methane_sym_num = automol.geom.external_symmetry_factor(methane_geo)
    h_sym_num = automol.geom.external_symmetry_factor(H_GEO)
    c2h5of_sym_num = automol.geom.external_symmetry_factor(c2h5of_geo)
    c2h2clf_sym_num = automol.geom.external_symmetry_factor(C2H2CLF_GEO)
    assert methane_sym_num == ref_sym_num1
    assert h_sym_num == ref_sym_num2
    assert c2h5of_sym_num == ref_sym_num3
    assert c2h2clf_sym_num == ref_sym_num4

    ref_end_sym_num1 = 1.0
    ref_end_sym_num2 = 1.0

    prop_geo = (('C', (-1.271060, -0.260158, 0.000015)),
                ('C', (0.000030, 0.588062, 0.000031)),
                ('H', (-1.310561, -0.906763, -0.884416)),
                ('H', (-2.170753, 0.364089, -0.000108)),
                ('H', (-1.310878, -0.906697, 0.884491)),
                ('C', (1.271039, -0.260179, 0.000001)),
                ('H', (0.000092, 1.246476, 0.877632)),
                ('H', (0.000037, 1.246109, -0.877909)),
                ('H', (1.310415, -0.906937, -0.884321)),
                ('H', (2.170786, 0.363983, -0.000215)),
                ('H', (1.310810, -0.906613, 0.884564)))

    prop_h_geo = (('H', (-2.925536, 0.892159, 0.000083)),
                  ('H', (-2.235843, 0.265746, 0.000083)),
                  ('C', (-1.087894, -0.458479, -0.000022)),
                  ('C', (0.048056, 0.535538, -0.000030)),
                  ('H', (-1.181615, -1.061889, -0.906084)),
                  ('H', (-1.181466, -1.061854, 0.906078)),
                  ('C', (1.420131, -0.148928, 0.000015)),
                  ('H', (-0.035403, 1.185804, -0.879125)),
                  ('H', (-0.035486, 1.185744, 0.879103)),
                  ('H', (1.542921, -0.783078, 0.884904)),
                  ('H', (1.543226, -0.782626, -0.885152)),
                  ('H', (2.227441, 0.591206, 0.000337)))

    _, end_sym_num1 = automol.geom.end_group_symmetry_factor(prop_geo)
    _, end_sym_num2 = automol.geom.end_group_symmetry_factor(
        prop_h_geo,
        frm_bnd_keys=frozenset({0, 1}),
        brk_bnd_keys=frozenset({1, 2}))

    assert numpy.isclose(end_sym_num1, ref_end_sym_num1)
    assert numpy.isclose(end_sym_num2, ref_end_sym_num2)


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


if __name__ == '__main__':
    test__from_data()
    # test__is_valid()
    # test__struct_check()
    # test__mass()
    # test__rotation_properties()
    # test__atom_indices()
    # test__set_coordinates()
    # test__swap_coordinates()
    test__remove_coordinates()
    # test__dist_analysis()
    # test__symmetry_factor()
    # test__closest_unbonded_atoms()
    # test__rotate()
    # test__reflect_coordinates()
    test__permutations()
    test__traj()
