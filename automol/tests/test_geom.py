""" test automol.geom
"""
import numpy
import automol
from automol import geom

C2H2CLF_GEO = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))


def test__from_data():
    """ test getters
    """
    assert C2H2CLF_GEO == geom.from_data(
        syms=geom.symbols(C2H2CLF_GEO),
        xyzs=geom.coordinates(C2H2CLF_GEO),
    )


def test__is_valid():
    """ test geom.is_valid
    """
    assert geom.is_valid(C2H2CLF_GEO)


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
        geo = geom.rotated(C2H2CLF_GEO, axis, angle)
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
        geo = geom.rotated(geo, axis, angle)

        if idx in ref_idxs and idx != 0:
            idx_to_change = numpy.random.randint(0, natms)
            new_xyz = numpy.random.rand(3)
            geo = geom.set_coordinates(geo, {idx_to_change: new_xyz})

        geos.append(geo)

    idxs = geom.argunique_coulomb_spectrum(geos)
    assert idxs == ref_idxs


def test__mass_centered():
    """ test geom.mass_centered()
    """
    # make sure the COM for an uncentered geometry is non-zero
    geo = C2H2CLF_GEO
    cm_xyz = automol.geom.center_of_mass(geo)
    assert not numpy.allclose(cm_xyz, 0.)

    # now make sure centering it yields a COM at the origin
    geo = automol.geom.mass_centered(geo)
    cm_xyz = automol.geom.center_of_mass(geo)
    assert numpy.allclose(cm_xyz, 0.)


def test__rotational_constants():
    """ test geom.rotational_constants()
    """
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


def test__move_coordinates():
    """ test geom.move_coordinates
    """
    ref_geo = (('F', (2.994881276150, -1.414434615111, -0.807144415388)),
               ('C', (-1.201356763194, -0.347546894407, -0.3408392500119)),
               ('Cl', (-3.027970874978, 1.39211904938, -0.0492290974807)),
               ('H', (1.731596406235, 2.324260256203, -0.4292070203467)),
               ('C', (1.170155936996, 0.359360756989, -0.513323178859)),
               ('H', (-1.66730598121, -2.31375855306, -0.433949091252)))
    mv_geo = automol.geom.move_coordinates(C2H2CLF_GEO, 1, 4)
    assert mv_geo == ref_geo


def test__external_symmetry_number():
    """ test geom.external_symmety_number
    """
    ref_sym_num1 = 12
    ref_sym_num2 = 2
    ref_sym_num3 = 0.5
    ref_sym_num4 = 1

    methane_geo = (('C', (1.2069668249, 1.9997649792, -0.0000004209)),
                   ('H', (3.3034303116, 1.9997688296, -0.0000006619)),
                   ('H', (0.5081497935, 0.1480118251, 0.6912478445)),
                   ('H', (0.5081445849, 2.3270011544, -1.9492882688)),
                   ('H', (0.5081426638, 3.5242779889, 1.2580393046)))
    water_geo = (('O', (-5.3344419198110174e-05, 0.7517614816502209, 0.0)),
                 ('H', (-1.4427184990730881, -0.3759830919088236, 0.0)),
                 ('H', (1.4427718434922905, -0.3757783897413961, 0.0)))
    c2h5of_geo = (('C', (-4.67963119210, -2.785693400767, -0.04102938592633)),
                  ('C', (-1.806009533535, -2.594940600449, -0.1025157659970)),
                  ('H', (-5.39544527869, -3.740953123044, -1.774996188159)),
                  ('H', (-5.501156952723, -0.854962636480, -0.01317371318990)),
                  ('O', (-5.48010155525, -4.07121874876, 2.132641999777597)),
                  ('H', (-1.208455201406, -1.52313066520, -1.804561025201)),
                  ('F', (-0.745999108314, -4.9827454242, -0.1878162481225)),
                  ('H', (-1.11738479998, -1.591763680324, 1.607521773704)),
                  ('H', (-5.30771129777, -5.90407965309, 1.771303996279)))
    methane_sym_num = automol.geom.external_symmetry_number(methane_geo)
    water_sym_num = automol.geom.external_symmetry_number(water_geo)
    c2h5of_sym_num = automol.geom.external_symmetry_number(c2h5of_geo)
    c2h2clf_sym_num = automol.geom.external_symmetry_number(C2H2CLF_GEO)
    assert methane_sym_num == ref_sym_num1
    assert water_sym_num == ref_sym_num2
    assert c2h5of_sym_num == ref_sym_num3
    assert c2h2clf_sym_num == ref_sym_num4

def test__rot_permutated_geoms():
    smi = 'C[CH2]'
    ich = automol.smiles.inchi(smi)
    geo = automol.inchi.geometry(ich)
    geos = automol.geom.rot_permutated_geoms(geo)
    for geom in geos:
         print(automol.geom.xyz_string(geom))

if __name__ == '__main__':
    # test__from_data()
    # test__is_valid()
    # test__set_coordinates()
    # test__swap_coordinates()
    # test__move_coordinates()
    #test__external_symmetry_number()
    test__rot_permutated_geoms()
