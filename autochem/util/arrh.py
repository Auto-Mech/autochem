from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import RegularGridInterpolator

from . import func


def interpolator(
    T: ArrayLike, P: ArrayLike, rate_data: ArrayLike, method: str = "linear"
) -> Callable[[ArrayLike, ArrayLike], NDArray[np.float64]]:
    """Define an Arrhenius interpolator.

    :param T: Temperatures
    :param P: Pressures
    :param rate_data: Rate data, as an array of k(T, P) values
    :raises ValueError: _description_
    :return: Interpolator
    """
    if np.less_equal(rate_data, 0).any():
        msg = "All values must be positive for logarithmic interpolation."
        raise ValueError(msg)

    points = (np.divide(1000.0, T), P)
    values = np.log10(rate_data)
    transformed_interp_ = RegularGridInterpolator(points, values, method=method)

    def interp_(T: ArrayLike, P: ArrayLike) -> NDArray[np.float64]:
        T_, P_ = func.normalize_arguments((T, P))
        x_ = np.divide(1000.0, T_)
        k = np.power(10, transformed_interp_((x_, P_)))
        return func.normalize_values(k, (T, P))

    return interp_
