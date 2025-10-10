"""Math utilities for working with scalar functions of multiple variables.

Example: The rate constant k(T, P)
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, DTypeLike, NDArray


def normalize_arguments(
    args: Sequence[ArrayLike], *, dtype: DTypeLike = np.float64
) -> tuple[NDArray, ...]:
    """Normalize arguments for ready formula evaluation.

    :param args: Arguments
    :return: Argument meshgrid
    """
    return np.meshgrid(*(np.array(arg, dtype=dtype) for arg in args))


def normalize_values(
    vals: ArrayLike, args: Sequence[ArrayLike], *, dtype: DTypeLike = np.float64
) -> NDArray:
    """Normalize values to match arguments.

    Ensures that shape of the values matches the arguments:

        Argument 1      | Argument 2        | Shape
        100.            | 1.                | ()
        [100., 200.]    | 1.                | (2,)
        100.            | [0.1, 1., 10.]    | (3,)
        [100., 200.]    | [0.1, 1., 10.]    | (2, 3)

    :param args: Arguments
    :param vals: Values
    :return: Values
    """
    shape = sum(map(np.shape, args), ())
    return np.reshape(vals, shape).astype(dtype)


def translation_indices(
    args1: Sequence[ArrayLike], args2: Sequence[ArrayLike]
) -> tuple[tuple[NDArray, ...], tuple[NDArray, ...]]:
    """Get indices for translating values between functions.

    Given function data `vals1` and `vals2` defined over arguments `args1` and
    `args2`, respectively, this provides the translation between matching
    arguments.

    For example, one could use the data from `vals2` to set the data in `vals1`:

        ix1, ix2 = translation_indices(args1, args2)
        vals1[ix1] = vals2[ix2]

    :param args1: Arguments
    :param args2: Arguments
    :return: Indices for args1, indices for args2
    """
    ix_pairs = []
    for arg1, arg2 in zip(args1, args2, strict=True):
        arg2_, arg1_ = np.meshgrid(arg2, arg1)
        ix_pairs.append(np.nonzero(np.isclose(arg1_, arg2_)))

    ix1s, ix2s = zip(*ix_pairs, strict=True)
    return np.ix_(*ix1s), np.ix_(*ix2s)
