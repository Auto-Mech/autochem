"""I/O operations for higher dimension matrices of 3 or higher."""

import itertools
from collections.abc import Sequence

import numpy

from . import vector


# I/O
def string(
    arr: numpy.ndarray,
    include_zeros: bool = False,
    include_perms: bool = False,
    val_format: str = "{0:>16.8f}",
) -> str:
    """Write a higher dimensional array (3 or more) to a string.

    :param arr: higher dimensions matrix
    :return: Numbers in a string

    Format of output string:
        idx1 idx2 .. idxn  val1
        idx1 idx2 .. idxn  valn
    """

    def _chk_zero(val: str, include_zeros: bool = False):
        """Decide to write value.
        :param val: A number
        :param include_zeroes:
        :return: True if it is zero and False if not.
        """
        iszero = numpy.isclose(val, 0.0)
        return bool(not iszero or (iszero and include_zeros))

    def _chk_idxs(arr: numpy.ndarray) -> str:
        """Decide if idxs.
        check if any permutation of idxs is already written.
        :param arr: Array of high dimension matrix
        :return: Index as a string.
        """
        idxs_lst = tuple(idxs for idxs, _ in numpy.ndenumerate(arr))
        if not include_perms:
            idxs_lst = {tuple(sorted(x)) for x in idxs_lst}
        vals_lst = tuple(arr[idxs] for idxs in idxs_lst)

        re_idxs = list(zip(idxs_lst, vals_lst, strict=False))
        fin_idxs = tuple(sorted(re_idxs, key=lambda x: x[0]))

        return fin_idxs

    arr_str = ""
    for idxs, val in _chk_idxs(arr):
        if _chk_zero(val, include_zeros):
            val_str = "".join(f"{idx+1:<6d}" for idx in idxs)
            val_str += val_format.format(val)
            val_str += "\n"

            arr_str += val_str

    arr_str.rstrip()

    return arr_str


def from_string(arr_str: str, fill_perms: bool = False) -> numpy.ndarray:
    """Write a higher dimensional array (3 or more) to a string.

    :param arr_str: string containing higher dimensions matrix
    :return: Array of high dimension matrix

    Format of input string:
        idx1 idx2 .. idxn  val1
        idx1 idx2 .. idxn  valn
    """
    lines = arr_str.splitlines()

    # Get the number of values in each array dimension; initialize array
    nvals = [int(val) for val in lines[-1].strip().split()[:-1]]
    arr = numpy.zeros(nvals)

    mat_idxs, mat_vals = [], []
    for line in lines:
        # Read the values from the line
        tmp = line.strip().split()
        idxs = tuple(int(val) - 1 for val in tmp[:-1])
        val = float(tmp[-1])
        # Store values
        mat_idxs.append(idxs)
        mat_vals.append(val)

    arr = build_full_array(mat_idxs, mat_vals, fill_perms=fill_perms)

    return arr


def string_submat_4d(arr: numpy.ndarray) -> str:
    """Writes a 4-dimensional matrix to a unique string format.
    :param arr:higher dimensions matrix
    :return: Numbers in a string
    idx1
      2d submat1
    idx2
      2d submat2
    idxn
      2d submatn.
    """  # noqa: D401
    nd1, nd2, nd3, nd4 = arr.shape
    idxs = tuple(val for val in range(nd1))

    arr_str = ""
    for idx in idxs:
        # Get the string for the sub 3n matrix
        sub_arr = arr[idx, :nd2, :nd3, :nd4]
        bmat_str = string_submat_3d(sub_arr)

        # Write index
        arr_str += f"{idx+1:>6d}"
        arr_str += bmat_str
        arr_str += "\n"

    return arr_str


def string_submat_3d(arr: numpy.ndarray) -> str:
    """Writes a 3-dimensional matrix to a unique string format.
    :param arr:higher dimensions matrix
    :return: Numbers in a string
    idx1
      2d submat1
    idx2
      2d submat2
    idxn
      2d submatn.

    """  # noqa: D401
    nd1, nd2, nd3 = arr.shape
    idxs = tuple(val for val in range(nd1))

    arr_str = ""
    for idx in idxs:
        sub_arr = arr[idx, :nd2, :nd3]
        # just get a single vector of elements and print to a line
        flat_sub_arr = sub_arr.flatten()
        sub_arr_str = vector.string(flat_sub_arr)
        arr_str += sub_arr_str
        arr_str += "\n"
    arr_str.rstrip()

    return arr_str


def build_full_array(
    mat_idxs: Sequence[tuple[int, ...], object],
    mat_vals: Sequence[float],
    fill_perms: bool = False,
) -> numpy.ndarray:
    """Function to fill out the array with avail
    caps: (
        ((idx1, idx2, ..., idxn), val1),
        ((idx1, idx2, ..., idxn), val2),
        ...,
        ((idx1, idx2, ..., idxn), valn),.
    :param mat_idxs: Matrix of the indez
    :param mat_vals: Matrix
    :param fill_perms: False if ...
    :return: Array.
    """  # noqa: D401

    def _gen_idxs(mat_idxs, mat_vals):
        """Permute idxs."""
        full_idxs, full_vals = [], []
        for idxs, val in zip(mat_idxs, mat_vals, strict=False):
            idx_perms = tuple(itertools.permutations(idxs))
            for perm in idx_perms:
                full_idxs.append(perm)
                full_vals.append(val)

        return full_idxs, full_vals

    # Build out the idxs and vals to include permutations, if needed
    if fill_perms:
        mat_idxs, mat_vals = _gen_idxs(mat_idxs, mat_vals)

    # Get dimensionality of matrix (assumes 0 idx of mat)
    ncoords = max(max(idxs) for idxs in mat_idxs) + 1
    ndim = len(mat_idxs[0])

    dims = tuple(ncoords for _ in range(ndim))

    # Build the force constant matrix
    mat = numpy.zeros(dims)
    if fill_perms:
        mat_idxs, mat_vals = _gen_idxs(mat_idxs, mat_vals)
    for idxs, val in zip(mat_idxs, mat_vals, strict=False):
        mat[idxs] = val

    return mat
