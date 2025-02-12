"""Unit utilities."""

import pint


def string(unit: pint.Unit) -> str:
    return format(unit, "~")
