""" test automol.cart not utilized by other tests
"""

import os
import pytest
import numpy
import automol.util.mat
import automol.util.vec
from ioformat import read_text_file


PATH = os.path.dirname(os.path.realpath(__file__))
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
    perp = automol.util.vec.unit_perpendicular(
        MAT[0], MAT[1], allow_parallel=False)
    assert numpy.allclose(perp, ref_perp)

    ref_perp = (0.000, 0.000, 0.000)
    perp = automol.util.vec.unit_perpendicular(
        MAT[0], MAT[0], allow_parallel=True)
    assert numpy.allclose(perp, ref_perp)

    with pytest.raises(ValueError):
        automol.util.vec.unit_perpendicular(
            MAT[0], MAT[0], allow_parallel=False)

    ref_angle = 0.28211376550390677
    angle = automol.util.vec.projected_central_angle(
        MAT[0], MAT[1], MAT[2])
    assert numpy.isclose(angle, ref_angle)

    # Test the string writer
    ref_vec_str = read_text_file(['data'], 'vec.dat', path=PATH)
    vec_str = automol.util.vec.string((MAT[0] + MAT[1]), num_per_row=3)

    assert vec_str == ref_vec_str


def test__mat():
    """ test automol.automol.util.cart.mat
    """

    # Various matrix builder functions
    rand_rot_mat = automol.util.mat.random_rotation_matrix()
    assert len(rand_rot_mat) == 3
    assert all(len(row) == 3 and all(isinstance(val, float) for val in row)
               for row in rand_rot_mat)

    ref_rot_mat = (
        (1.0, 0.0, 0.0),
        (0.0, -0.9922575676015892, 0.12419709955299955),
        (0.0, -0.12419709955299955, -0.9922575676015892)
    )
    rot_mat = automol.util.mat.rotation_matrix(
        (1.0, 0.0, 0.0), 30.0/numpy.pi)

    assert numpy.allclose(rot_mat, ref_rot_mat)

    ref_axis_align_mat = (
        (-0.41149979255451374, 0.9113875155894119, 0.006380999557419693),
        (0.0, 0.0, 0.0),
        (-0.9018068740366958, -0.4061404341551566, -0.14762895950464514)
    )
    axis_align_mat = automol.util.mat.axis_alignment_matrix(
        MAT[0], MAT[1])

    assert numpy.allclose(axis_align_mat, ref_axis_align_mat)

    ref_superimp_mat = (
        (0.9931073252359884, 0.07107456632823939, 0.0932000353212035),
        (-0.05722539832771756, 0.9879718232559653, -0.1436555959178673),
        (-0.10228926800429702, 0.13733201547176194, 0.9852293251713573)
    )

    xyz1 = (MAT[0], MAT[1])
    xyz2 = (MAT[2], MAT[3])
    superimp_mat = automol.util.mat.superimposition_matrix(
        xyz1, xyz2, keep_origin=True)

    assert numpy.allclose(superimp_mat, ref_superimp_mat)

    # Test the string writer
    ref_mat_str = read_text_file(['data'], 'mat.dat', path=PATH)
    mat_str = automol.util.mat.string(MAT)

    assert mat_str == ref_mat_str


def test__highd_mat():
    """ test automol.util.highd_mat.from_string
        test automol.util.highd_mat.string
        test automol.util.highd_mat.string_submat_3d
        test automol.util.highd_mat.string_submat_4d
    """

    def _chk_mat_strs(str1, str2):
        """ Check if two matrix strings are similar
        """
        match = True
        for line1, line2 in zip(str1.splitlines(), str2.splitlines()):
            vals1 = tuple(float(val) for val in line1.strip().split())
            vals2 = tuple(float(val) for val in line2.strip().split())
            if not numpy.allclose(vals1, vals2):
                match = False
        return match

    ref_3d_str = read_text_file(['data'], 'ch4_h.cubic', path=PATH)
    ref_4d_str = read_text_file(['data'], 'ch4_h.quartic', path=PATH)

    # Handle reprentations with full matrices and strings
    test_3d_mat = automol.util.highd_mat.from_string(ref_3d_str)
    test_4d_mat = automol.util.highd_mat.from_string(ref_4d_str)

    test_3d_str = automol.util.highd_mat.string(test_3d_mat)
    test_4d_str = automol.util.highd_mat.string(test_4d_mat)

    assert _chk_mat_strs(test_3d_str, ref_3d_str)
    assert _chk_mat_strs(test_4d_str, ref_4d_str)

    # Handle string representations by submatrices (finish)
    test_3d_submat_str = automol.util.highd_mat.string_submat_3d(test_3d_mat)
    assert (test_3d_submat_str ==
            read_text_file(['data'], 'ch4_h.cubic_submat', path=PATH))

    test_4d_submat_str = automol.util.highd_mat.string_submat_4d(test_4d_mat)
    assert (test_4d_submat_str ==
            read_text_file(['data'], 'ch4_h.quartic_submat', path=PATH))
