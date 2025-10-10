"""Utilities for working with Arrhenius function data."""

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import RegularGridInterpolator

from . import func


def interpolator(
    T: ArrayLike,  # noqa: N803
    P: ArrayLike,  # noqa: N803
    rate_data: NDArray,
    *,
    method: str = "linear",
    bounds_error: bool = False,
) -> Callable[[ArrayLike, ArrayLike], NDArray[np.float64]]:
    """Define an Arrhenius interpolator.

    :param T: Temperatures
    :param P: Pressures
    :param rate_data: Rate data, as an array of k(T, P) values
    :param method: Interpolation method
        (see scip.interpolate.RegularGridInterpolator)
    :param bounds_error: Whether to raise bounds errors
        (see scip.interpolate.RegularGridInterpolator)
    Whether to raise an error when out of bonds
    :raises ValueError: _description_
    :return: Interpolator
    """
    points = (np.divide(1000.0, T), P)
    values = np.log10(rate_data)
    transformed_interp_ = RegularGridInterpolator(
        points, values, method=method, bounds_error=bounds_error
    )

    def interp_(t: ArrayLike, p: ArrayLike) -> NDArray[np.float64]:
        # 1. Perform interpolation
        t_, p_ = func.normalize_arguments((t, p))
        x_ = np.divide(1000.0, t_)
        rate_interp = np.power(10, transformed_interp_((x_, p_)).T)
        # 2. Place original data where arguments match
        # (Floating point precision may cause some of these to be missing)
        ix_data, ix_interp = func.translation_indices((T, P), (t, p))
        rate_interp[ix_interp] = rate_data[ix_data]
        return func.normalize_values(rate_interp, (t, p))

    return interp_
