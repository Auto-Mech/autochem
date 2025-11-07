"""Plotting helpers."""

import itertools
from collections.abc import Callable, Sequence
from typing import Any

import altair as alt
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.interpolate import CubicSpline

from .. import unit_
from ..unit_ import UNITS, Units, UnitsData


class Color:
    """Color hex values."""

    # Core colors:
    blue = "#0066ff"
    red = "#ff0000"
    green = "#1ab73a"
    orange = "#ef7810"
    purple = "#8533ff"
    pink = "#d0009a"
    yellow = "#ffcd00"
    # Extra colors:
    teal = "#008080"
    cyan = "#00ffff"
    magenta = "#ff00ff"
    lime = "#00ff00"
    navy = "#000080"
    maroon = "#800000"
    olive = "#808000"
    coral = "#ff7f50"
    gold = "#ffd700"
    sky_blue = "#87ceeb"
    violet = "#ee82ee"
    indigo = "#4b0082"
    salmon = "#fa8072"
    mint = "#98ff98"
    peach = "#ffdab9"
    forest_green = "#228b22"
    mustard = "#ffdb58"
    steel_blue = "#4682b4"
    plum = "#dda0dd"
    ochre = "#cc7722"
    # Point colors:
    black = "#000000"
    gray = "#808080ff"
    light_gray = "#bfbfbfff"
    brown = "#916e6e"
    sienna = "#a0522d"

    # ChatGPT generated colors:
    # # --- Core colors (your originals) ---
    # blue = "#0066ff"
    # red = "#ff0000"
    # green = "#1ab73a"
    # orange = "#ef7810"
    # purple = "#8533ff"
    # pink = "#d0009a"
    # yellow = "#ffcd00"
    # black = "#000000"
    # gray = "#808080"
    # light_gray = "#bfbfbf"
    # brown = "#916e6e"

    # # --- Strong blues ---
    # navy = "#003f5c"
    # royal_blue = "#4169e1"
    # sky_blue = "#1ca8dd"
    # teal = "#008080"
    # turquoise = "#00a0b0"

    # # --- Strong reds / magentas ---
    # crimson = "#dc143c"
    # brick = "#b22222"
    # carmine = "#960018"
    # magenta = "#ff00ff"
    # raspberry = "#b03060"

    # # --- Oranges & yellows ---
    # amber = "#ffbf00"
    # burnt_orange = "#cc5500"
    # gold = "#d4a017"
    # ochre = "#c07a28"
    # mustard = "#e1ad01"

    # # --- Greens ---
    # forest_green = "#228b22"
    # emerald = "#50c878"
    # lime_green = "#32cd32"
    # olive = "#6b8e23"
    # jade = "#00a86b"

    # # --- Purples & violets ---
    # indigo = "#4b0082"
    # violet = "#7b68ee"
    # plum = "#8e4585"
    # mauve = "#7d5ba6"
    # orchid = "#da70d6"

    # # --- Cyans & aquas ---
    # cyan = "#17becf"
    # sea_green = "#20b2aa"
    # steel_blue = "#4682b4"
    # cerulean = "#007ba7"
    # azure = "#007fff"

    # # --- Earth tones & neutrals ---
    # sienna = "#a0522d"
    # copper = "#b87333"
    # chocolate = "#7b3f00"
    # slate_gray = "#708090"
    # charcoal = "#36454f"

    # # --- Extras for balance ---
    # coral = "#ff7f50"
    # maroon = "#800000"
    # olive_drab = "#556b2f"
    # midnight_blue = "#191970"
    # royal_purple = "#7851a9"


LINE_COLOR_CYCLE = [
    Color.blue,
    Color.red,
    Color.green,
    Color.purple,
    Color.pink,
    Color.yellow,
    Color.orange,
    Color.teal,
    Color.cyan,
    Color.magenta,
    Color.lime,
    Color.navy,
    Color.maroon,
    Color.olive,
    Color.coral,
    Color.gold,
    Color.sky_blue,
]


POINT_COLOR_CYCLE = [
    Color.black,
    Color.gray,
    Color.light_gray,
    Color.brown,
]


class Mark:
    """Altair mark types."""

    point = "point"
    line = "line"


def regular_scale(val_range: tuple[float, float]) -> alt.Scale:
    """Generate a regular scale specification.

    :param val_range: Range
    :return: Scale
    """
    return alt.Scale(domain=val_range)


def regular_scale_axis(val_range: tuple[float, float]) -> alt.Axis:
    """Generate a nice regular scale axis.

    :param val_range: Range
    :return: Axis
    """
    val_min, val_max = val_range
    val_scale = val_max - val_min
    if val_scale < 1:
        fmt = ".2f"
    elif val_scale < 3:
        fmt = ".1f"
    else:
        fmt = ".0f"
    return alt.Axis(format=fmt)


def log_scale(val_range: tuple[float, float]) -> alt.Scale:
    """Generate a log scale specification.

    :param val_range: Range
    :return: Scale
    """
    return alt.Scale(type="log", domain=log_scale_domain(val_range))


def log_scale_axis(val_range: tuple[float, float]) -> alt.Axis:
    """Generate a nice log scale axis.

    :param val_range: Range
    :return: Axis
    """
    max_exp = np.max(np.abs(np.log10(val_range)))
    fmt = ".0e" if max_exp > 3 else alt.Undefined
    vals = log_scale_values(val_range)
    label_expr_condition = " ||\n".join(
        f"(abs(datum.value - {v}) / abs({v}) < 1e-5)" for v in vals
    )
    label_expr = f"({label_expr_condition}) ? datum.label : ''"
    return alt.Axis(format=fmt, values=log_scale_ticks(val_range), labelExpr=label_expr)


def log_scale_domain(val_range: tuple[float, float]) -> tuple[float, float]:
    """Determine log scale ticks for a given range.

    :param val_range: Range
    :return: Ticks
    """
    # Determine tick min and max
    val_min, val_max = val_range
    mant_min, exp_min = decompose_base10(val_min)
    mant_max, exp_max = decompose_base10(val_max)
    start = recompose_base10(mant=np.floor(mant_min), exp=exp_min)
    stop = recompose_base10(mant=np.ceil(mant_max), exp=exp_max)
    return start, stop


def log_scale_values(val_range: tuple[float, float]) -> list[float]:
    """Determine log scale ticks for a given range.

    :param val_range: Range
    :return: Ticks
    """
    val_min, val_max = log_scale_domain(val_range)
    _, exp_min = decompose_base10(val_min)
    _, exp_max = decompose_base10(val_max)

    # Add power of 10 steps in between
    exp_start = exp_min + 1
    exp_stop = exp_max
    exp_count = exp_stop - exp_start + 1
    powers_of_10 = []
    if exp_count > 0:
        powers_of_10 = np.logspace(
            exp_start, exp_stop, num=exp_count, endpoint=True
        ).tolist()
    return [val_min, *powers_of_10, val_max]


def log_scale_ticks(val_range: tuple[float, float]) -> list[float]:
    """Determine log scale ticks for a given range.

    :param val_range: Range
    :return: Ticks
    """
    vals = log_scale_values(val_range)
    if len(vals) > 5:
        return vals

    # If the scale is not too large, add intervening ticks
    bounds = [*vals, None]
    ticks = []
    for start, stop in itertools.pairwise(bounds):
        ticks.append(start)
        if stop is not None:
            _, exp = decompose_base10(start)
            vals = [recompose_base10(m, exp) for m in range(2, 10)]
            vals = [v for v in vals if start < v and v < stop]
            ticks.extend(vals)
    return ticks


def decompose_base10(val: float) -> tuple[float, int]:
    """Decompose value into base-10 mantissa and exponent.

    :param val: Value
    :return: Mantissa and exponent
    """
    exp = np.floor(np.log10(val)).astype(int)
    mant = val / (10.0**exp)
    return float(mant), int(exp)


def recompose_base10(mant: float, exp: int) -> float:
    """Recompose value from base-10 mantissa and exponent

    :param mant: Mantissa
    :param exp: Exponent
    :return: Value
    """
    return float(mant * 10.0**exp)


MARKS = (Mark.point, Mark.line)


def transformed_spline_interpolator(
    x_data: ArrayLike,
    y_data: ArrayLike,
    x_trans: Callable[[ArrayLike], ArrayLike] = lambda x: x,
    y_trans: Callable[[ArrayLike], ArrayLike] = lambda y: y,
    y_trans_inv: Callable[[ArrayLike], ArrayLike] = lambda y: y,
) -> Callable[[Any], np.ndarray]:
    """Generate an inerpolator from data.

    :param y_data: Y data
    :param x_data: X data
    :return: Y interpolator
    """
    valid = np.isfinite(y_data)
    x_data = np.compress(valid, x_data)
    y_data = np.compress(valid, y_data)
    interp_trans_ = CubicSpline(x_trans(x_data), y_trans(y_data))

    def interp_(x: Any) -> np.ndarray:
        return np.asarray(y_trans_inv(interp_trans_(x_trans(x))))

    return interp_


def general(
    y_data: Sequence[Sequence[float]],
    x_data: Sequence[float],  # noqa: N803
    labels: Sequence[str],
    *,
    colors: Sequence[str] | None = None,
    x_label: str | None = None,  # noqa: RUF001
    y_label: str | None = None,  # noqa: RUF001
    x_scale: alt.Scale | None = None,
    y_scale: alt.Scale | None = None,
    x_axis: alt.Axis | None = None,
    y_axis: alt.Axis | None = None,
    mark: str = Mark.line,
    mark_kwargs: dict | None = None,
    legend: bool = True,
) -> alt.Chart:
    """Display as simple plot.

    :return: Chart
    """
    x_label = "" if x_label is None else x_label
    y_label = "" if y_label is None else y_label
    x_scale_ = alt.Undefined if x_scale is None else x_scale
    y_scale_ = alt.Undefined if y_scale is None else y_scale
    x_axis_ = alt.Undefined if x_axis is None else x_axis
    y_axis_ = alt.Undefined if y_axis is None else y_axis

    assert mark in MARKS, f"{mark} not in {MARKS}"
    color_cycle = (
        LINE_COLOR_CYCLE
        if mark == Mark.line
        else [*POINT_COLOR_CYCLE, *LINE_COLOR_CYCLE]
    )

    ny, nx = np.shape(y_data)  # noqa: N806
    colors = colors or list(itertools.islice(itertools.cycle(color_cycle), ny))
    assert len(x_data) == nx, f"{x_data} !~ {y_data}"
    assert len(labels) == ny, f"{labels} !~ {y_data}"

    # Gather data from functons
    data_dct = dict(zip(labels, y_data, strict=True))
    data = pd.DataFrame({"x": x_data, **data_dct})

    # Prepare encoding parameters
    x = alt.X("x", title=x_label, scale=x_scale_, axis=x_axis_)
    y = alt.Y("value:Q", title=y_label, scale=y_scale_, axis=y_axis_)
    color = alt.Color(
        "key:N",
        scale=alt.Scale(domain=labels, range=colors),
        legend=alt.Undefined if legend else None,
    )

    chart = alt.Chart(data)
    kwargs = {} if mark_kwargs is None else mark_kwargs
    if mark == Mark.point:
        chart = chart.mark_point(**kwargs)
    else:
        chart = chart.mark_line(**kwargs)

    # Create chart
    return chart.transform_fold(fold=list(data_dct.keys())).encode(
        x=x, y=y, color=color
    )


def simple(
    ks: Sequence[Sequence[float]],
    T: Sequence[float],  # noqa: N803
    *,
    order: int = 1,
    units: UnitsData | None = None,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    x_label: str | None = "𝑇",  # noqa: RUF001
    y_label: str | None = "𝑘",  # noqa: RUF001
    x_unit: str | None = None,
    y_unit: str | None = None,
    mark: str = Mark.line,
) -> alt.Chart:
    """Display as simple plot.

    :param others: Other rate constants
    :param others_labels: Labels for other rate constants
    :param T_range: Temperature range
    :param P: Pressure
    :param units: Units
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :param point: Whether to mark with points instead of a line
    :return: Chart
    """
    x_label = "𝑇" if x_label is None else x_label
    y_label = "𝑘" if y_label is None else y_label

    assert mark in MARKS, f"{mark} not in {MARKS}"
    color_cycle = (
        LINE_COLOR_CYCLE
        if mark == Mark.line
        else [*POINT_COLOR_CYCLE, *LINE_COLOR_CYCLE]
    )

    nk, nT = np.shape(ks)  # noqa: N806
    colors = colors or list(itertools.islice(itertools.cycle(color_cycle), nk))
    keep_legend = labels is not None
    labels = labels or [f"k{i + 1}" for i in range(nk)]
    assert len(T) == nT, f"{T} !~ {ks}"
    assert len(labels) == nk, f"{labels} !~ {ks}"

    # Process units
    units = UNITS if units is None else Units.model_validate(units)
    x_unit = unit_.pretty_string(units.temperature) if x_unit is None else x_unit
    y_unit = (
        unit_.pretty_string(units.rate_constant(order)) if y_unit is None else y_unit
    )

    # Add units to labels
    if x_unit:
        x_label = f"{x_label} ({x_unit})"

    if y_unit:
        y_label = f"{y_label} ({y_unit})"

    # Gather data from functons
    data_dct = dict(zip(labels, ks, strict=True))
    data = pd.DataFrame({"x": T, **data_dct})

    # Prepare encoding parameters
    x = alt.X("x", title=x_label, scale=alt.Scale(zero=False))
    y = alt.Y("value:Q", title=y_label)
    color = (
        alt.Color("key:N", scale=alt.Scale(domain=labels, range=colors))
        if keep_legend
        else alt.value(colors[0])
    )

    chart = alt.Chart(data)
    chart = (
        chart.mark_point(filled=True, opacity=1)
        if mark == Mark.point
        else chart.mark_line()
    )

    # Create chart
    return chart.transform_fold(fold=list(data_dct.keys())).encode(
        x=x,
        y=y,
        color=color,
    )


def arrhenius(  # noqa: PLR0913
    ks: ArrayLike,
    T: Sequence[float],  # noqa: N803
    *,
    order: int = 1,
    units: UnitsData | None = None,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    x_label: str | None = "1000/𝑇",  # noqa: RUF001
    y_label: str | None = "𝑘",  # noqa: RUF001
    x_unit: str | None = None,
    y_unit: str | None = None,
    mark: str = Mark.line,
    mark_kwargs: dict | None = None,
) -> alt.Chart:
    """Display as Arrhenius plot.

    :param others: Other rate constants
    :param others_labels: Labels for other rate constants
    :param T_range: Temperature range
    :param P: Pressure
    :param units: Units
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :param point: Whether to mark with points instead of a line
    :return: Chart
    """
    x_label = "1000/𝑇" if x_label is None else x_label
    y_label = "𝑘" if y_label is None else y_label

    assert mark in MARKS, f"{mark} not in {MARKS}"
    color_cycle = LINE_COLOR_CYCLE if mark == Mark.line else POINT_COLOR_CYCLE

    ks = np.where(np.less_equal(ks, 0), np.nan, ks)
    nk, nT = np.shape(ks)  # noqa: N806
    colors = colors or list(itertools.islice(itertools.cycle(color_cycle), nk))
    keep_legend = labels is not None
    labels = labels or [f"k{i + 1}" for i in range(nk)]
    assert len(T) == nT, f"{T} !~ {ks}"
    assert len(labels) == nk, f"{labels} !~ {ks}"

    # Process units
    units = UNITS if units is None else Units.model_validate(units)
    x_unit = unit_.pretty_string(units.temperature**-1) if x_unit is None else x_unit
    y_unit = (
        unit_.pretty_string(units.rate_constant(order)) if y_unit is None else y_unit
    )

    # Add units to labels
    if x_unit:
        x_label = f"{x_label} ({x_unit})"

    if y_unit:
        y_label = f"{y_label} ({y_unit})"

    # Gather data from functons
    data_dct = dict(zip(labels, ks, strict=True))
    data = pd.DataFrame({"x": np.divide(1000, T), **data_dct})

    # Determine exponent range
    vals_arr = np.array(list(data_dct.values()))
    y_range = (np.nanmin(vals_arr), np.nanmax(vals_arr))

    # Prepare encoding parameters
    x = alt.X("x", title=x_label, scale=alt.Scale(zero=False))
    y = alt.Y(
        "value:Q",
        title=y_label,
        scale=log_scale(y_range),
        axis=log_scale_axis(y_range),
    )
    color = alt.Color(
        "key:N",
        scale=alt.Scale(domain=labels, range=colors),
        legend=alt.Undefined if keep_legend else None,
    )

    chart = alt.Chart(data)
    if mark == Mark.point:
        kwargs = {"filled": True, "opacity": 1} if mark_kwargs is None else mark_kwargs
        chart = chart.mark_point(**kwargs)
    else:
        kwargs = {} if mark_kwargs is None else mark_kwargs
        chart = chart.mark_line(**kwargs)

    # Create chart
    return chart.transform_fold(fold=list(data_dct.keys())).encode(
        x=x,
        y=y,
        color=color,
    )
