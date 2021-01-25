""" test automol.cart not utilized by other tests
"""

import pytest
import numpy
import automol.util.mat
import automol.util.vec
from _util import read_file


MAT = (
    (-2.3779010433, 5.2665623735, 0.0368733734),
    (-1.7871641824, 4.2084900234, -0.6608528628),
    (-1.9843935085, 5.3116891951, 1.0755446343),
    (-2.1223381321, 6.2190137934, -0.4713637932),
    (-2.0307815005, 3.3844774267, -0.1645507122),
    (-3.4841693472, 5.1584160162, 0.0370840829)
)


def test__vec():
    """ test automol.automol.util.cart.vec
    """

    # Test angle calculators
    ref_perp = (-0.90180687, -0.40614043, -0.14762896)
    perp = automol.util.cart.vec.unit_perpendicular(
        MAT[0], MAT[1], allow_parallel=False)
    assert numpy.allclose(perp, ref_perp)
    ref_perp = (0.000, 0.000, 0.000)
    perp = automol.util.cart.vec.unit_perpendicular(
        MAT[0], MAT[1], allow_parallel=True)
    assert numpy.allclose(perp, ref_perp)
    with pytest.raises(ValueError):
        automol.util.cart.vec.unit_perpendicular(
            MAT[0], MAT[0], allow_parallel=False)

    ref_angle = 0.28211376550390677
    angle = automol.util.cart.vec.projected_central_angle(
        MAT[0], MAT[1], MAT[2])
    assert numpy.isclose(angle, ref_angle)

    # Test the string writer
    ref_vec_str = read_file(['data'], 'vec.dat')
    vec_str = automol.util.cart.vec.string((MAT[0] + MAT[1]), num_per_row=3)

    assert vec_str == ref_vec_str


def test__mat():
    """ test automol.automol.util.cart.mat
    """

    # Various matrix builder functions
    ref_dist_mat = (
        (0.0, 1.3983236526476774, 1.1116306032993546, 1.1094057537096793,
         1.9243979282426587, 1.1115418297165136),
        (1.3983236526476774, 0.0, 2.0666456124475165, 2.0470597240408956,
         0.9923013563876134, 2.0662289082596796),
        (1.1116306032993546, 2.0666456124475165, 0.0,
         1.7986640400394194, 2.2921896768095347, 1.8306339204036375),
        (1.1094057537096793, 2.0470597240408956, 1.7986640400394194,
         0.0, 2.852562584345432, 1.7994363412282537),
        (1.9243979282426587, 0.9923013563876134, 2.2921896768095347,
         2.852562584345432, 0.0, 2.30214051327184),
        (1.1115418297165136, 2.0662289082596796, 1.8306339204036375,
         1.7994363412282537, 2.30214051327184, 0.0))
    dist_mat = automol.util.cart.mat.distance_matrix(MAT)
    assert numpy.allclose(dist_mat, ref_dist_mat)

    rand_rot_mat = automol.util.cart.mat.random_rotation_matrix()
    assert len(rand_rot_mat) == 3
    assert all(len(row) == 3 and all(isinstance(val, float) for val in row)
               for row in rand_rot_mat)

    ref_rot_mat = (
        (1.0, 0.0, 0.0),
        (0.0, -0.9922575676015892, 0.12419709955299955),
        (0.0, -0.12419709955299955, -0.9922575676015892)
    )
    rot_mat = automol.util.cart.mat.rotation_matrix(
        (1.0, 0.0, 0.0), 30.0/numpy.pi)

    assert numpy.allclose(rot_mat, ref_rot_mat)

    ref_axis_align_mat = (
        (-0.41149979255451374, 0.9113875155894119, 0.006380999557419693),
        (0.0, 0.0, 0.0),
        (-0.9018068740366958, -0.4061404341551566, -0.14762895950464514)
    )
    axis_align_mat = automol.util.cart.mat.axis_alignment_matrix(
        MAT[0], MAT[1])

    assert numpy.allclose(axis_align_mat, ref_axis_align_mat)

    ref_superimp_mat = (
        (0.9931073252359884, 0.07107456632823939, 0.0932000353212035),
        (-0.05722539832771756, 0.9879718232559653, -0.1436555959178673),
        (-0.10228926800429702, 0.13733201547176194, 0.9852293251713573)
    )

    xyz1 = (MAT[0], MAT[1])
    xyz2 = (MAT[2], MAT[3])
    superimp_mat = automol.util.cart.mat.superimposition_matrix(
        xyz1, xyz2, keep_origin=True)

    assert numpy.allclose(superimp_mat, ref_superimp_mat)

    # Test the string writer
    ref_mat_str = read_file(['data'], 'mat.dat')
    mat_str = automol.util.cart.mat.string(MAT)

    assert mat_str == ref_mat_str
