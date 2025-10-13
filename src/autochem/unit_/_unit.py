"""Unit utilities."""

import pint


def string(unit: pint.Unit) -> str:
    return format(unit, "~C")


def pretty_string(unit: pint.Unit) -> str:
    return format(unit, "~P")
