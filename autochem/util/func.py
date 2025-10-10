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
