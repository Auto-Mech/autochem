"""Plotting functions."""

from collections.abc import Mapping

import altair
import numpy
import pandas
import pint
from numpy.typing import ArrayLike


def arrhenius_plot(
    T: ArrayLike,  # noqa: N803
    kT_: ArrayLike | Mapping[str, ArrayLike],  # noqa: N803
    x_label: str = "1000/T",
    x_unit: str | pint.Unit = "1/K",
    y_label: str = "k",
    y_unit: str | pint.Unit | None = None,
) -> altair.Chart:
    """Generate an Arrhenius plot.

    :param T: Temperature array
    :param kT_: Rate array(s), as a dictionary by legend label if multiple
    :param x_label: X-axis label
    :param x_unit: X-axis unit
    :param y_label: Y-axis label
    :param y_unit: Y-axis unit
    :return: Chart
    """
    x_unit = format(pint.Unit(x_unit), "~P")
    y_unit = None if y_unit is None else format(pint.Unit(y_unit), "~P")
    x_label = f"{x_label} ({x_unit})"
    y_label = y_label if y_unit is None else f"{y_label} ({y_unit})"

    # Add data to DataFrame, dropping negative values for log plot
    has_series_labels = isinstance(kT_, Mapping)
    kTs = kT_ if isinstance(kT_, Mapping) else {"y": kT_}
    kTs = {lb: numpy.where(numpy.less(a, 0), numpy.nan, a) for lb, a in kTs.items()}
    data = pandas.DataFrame({"x": numpy.divide(1000, T), **kTs})

    # Determine the exponent range
    k_array = numpy.array(list(kTs.values()))
    is_nan = numpy.isnan(k_array)
    e_array = numpy.log10(k_array, where=~is_nan)
    e_array[is_nan] = 0.0
    e_array = numpy.rint(e_array).astype(int)
    e_max = numpy.max(e_array).item()
    e_min = numpy.min(e_array).item()
    values = [10**x for x in range(e_min, e_max + 2)]

    # Prepare encoding parameters
    x = altair.X("x", title=x_label)
    y = (
        altair.Y("value:Q", title=y_label)
        .scale(type="log")
        .axis(format=".1e", values=values)
    )
    color = "key:N" if has_series_labels else altair.Undefined

    # Create chart
    return (
        altair.Chart(data)
        .mark_line()
        .transform_fold(fold=list(kTs.keys()))
        .encode(x=x, y=y, color=color)
    )
