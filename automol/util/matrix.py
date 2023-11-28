""" 3D matrices
"""


# I/O
def string(mat, val_format="{0:>8.3f}"):
    """Write a matrix to a string.

    :param mat: matrix to form string with
    :type mat: tuple(tuple(float))
    :param precision: number of integers past decimal
    :type precision: int
    :rtype: str
    """

    mat_str = ""
    for row in mat:
        mat_str += "".join(val_format.format(val) for val in row)
        mat_str += "\n"
    mat_str = mat_str.rstrip()

    return mat_str
