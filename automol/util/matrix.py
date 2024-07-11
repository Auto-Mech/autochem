"""3D matrices."""
from collections.abc import Sequence  # type: ignore


# I/O
def string(mat: Sequence[Sequence[float]], val_format: str = "{0:>8.3f}") -> str:
    """Write a matrix to a string.

    :param mat: Matrix to form string with
    :param val_format: A number-formatting string, such as "{:.3f}"
    :return: Matrix as a string
    """
    mat_str = ""
    for row in mat:
        mat_str += "".join(val_format.format(val) for val in row)
        mat_str += "\n"
    mat_str = mat_str.rstrip()

    return mat_str
