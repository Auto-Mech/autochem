""" test automol.cart not utilized by other tests
"""

import os

import numpy
from automol.util import matrix, tensor, vector

PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, "data")

MAT = (
    (-2.3779010433, 5.2665623735, 0.0368733734),
    (-1.7871641824, 4.2084900234, -0.6608528628),
    (-1.9843935085, 5.3116891951, 1.0755446343),
    (-2.1223381321, 6.2190137934, -0.4713637932),
    (-2.0307815005, 3.3844774267, -0.1645507122),
    (-3.4841693472, 5.1584160162, 0.0370840829),
)

VEC_STR = """  -2.378   5.267   0.037
  -1.787   4.208  -0.661"""

MAT_STR = """  -2.378   5.267   0.037
  -1.787   4.208  -0.661
  -1.984   5.312   1.076
  -2.122   6.219  -0.471
  -2.031   3.384  -0.165
  -3.484   5.158   0.037"""


def test__vec():
    """test automol.automol.util.cart.vec"""

    # Test angle calculators
    ref_perp = (-0.90180687, -0.40614043, -0.14762896)
    perp = vector.unit_perpendicular(MAT[0], MAT[1])
    assert numpy.allclose(perp, ref_perp)

    ref_perp = (0.000, 0.000, 0.000)
    perp = vector.unit_perpendicular(MAT[0], MAT[0])
    assert numpy.allclose(perp, ref_perp)

    ref_angle = 0.28211376550390677
    angle = vector.projected_central_angle(MAT[0], MAT[1], MAT[2])
    assert numpy.isclose(angle, ref_angle)
    vec_str = vector.string((MAT[0] + MAT[1]), num_per_row=3)

    assert vec_str == VEC_STR


def test__mat():
    """test automol.automol.util.cart.mat"""

    mat_str = matrix.string(MAT)

    assert mat_str == MAT_STR


def test__highd_mat():
    """test automol.util.highd_mat.from_string
    test automol.util.highd_mat.string
    test automol.util.highd_mat.string_submat_3d
    test automol.util.highd_mat.string_submat_4d
    """

    def _chk_mat_strs(str1, str2):
        """Check if two matrix strings are similar"""
        match = True
        for line1, line2 in zip(str1.splitlines(), str2.splitlines()):
            vals1 = tuple(float(val) for val in line1.strip().split())
            vals2 = tuple(float(val) for val in line2.strip().split())
            if not numpy.allclose(vals1, vals2):
                match = False
        return match

    # Read formatted files
    ref_3d_str = None
    with open(
        os.path.join(DAT_PATH, "ch4_h.cubic"),
        mode="r",
        encoding="utf-8",
    ) as fobj:
        ref_3d_str = fobj.read()

    ref_4d_str = None
    with open(
        os.path.join(DAT_PATH, "ch4_h.quartic"),
        mode="r",
        encoding="utf-8",
    ) as fobj:
        ref_4d_str = fobj.read()

    # Handle reprentations with full matrices and strings
    test_3d_mat = tensor.from_string(ref_3d_str)
    test_4d_mat = tensor.from_string(ref_4d_str)

    test_3d_str = tensor.string(test_3d_mat)
    test_4d_str = tensor.string(test_4d_mat)

    assert _chk_mat_strs(test_3d_str, ref_3d_str)
    assert _chk_mat_strs(test_4d_str, ref_4d_str)

    # Handle string representations by submatrices (finish)
    ref_3d_submat_str = None
    with open(
        os.path.join(DAT_PATH, "ch4_h.cubic_submat"),
        mode="r",
        encoding="utf-8",
    ) as fobj:
        ref_3d_submat_str = fobj.read()

    ref_4d_submat_str = None
    with open(
        os.path.join(DAT_PATH, "ch4_h.quartic_submat"),
        mode="r",
        encoding="utf-8",
    ) as fobj:
        ref_4d_submat_str = fobj.read()

    test_3d_submat_str = tensor.string_submat_3d(test_3d_mat)
    assert test_3d_submat_str == ref_3d_submat_str

    test_4d_submat_str = tensor.string_submat_4d(test_4d_mat)
    assert test_4d_submat_str == ref_4d_submat_str


if __name__ == "__main__":
    test__vec()
    # test__mat()
    # test__highd_mat()