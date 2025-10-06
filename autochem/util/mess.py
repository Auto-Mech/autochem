"""Utility functions for reading and writing MESS input/output data."""

import numpy as np
import pydantic
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from .type_ import NDArray_


class MessOutputChannelParseResults(pydantic.BaseModel):
    """Mess output channel parse results."""

    T: list[float]
    P: list[float]
    k_data: NDArray_
    k_high: NDArray_ | None = None
    id1: str | None = None
    id2: str | None = None


class Key:
    """Pyparsing token keys."""

    eq = "eq"
    T = "T"
    data = "data"
    high = "high"


NAN = pp.Keyword("***").set_parse_action(pp.replace_with(np.nan))
NUMBERS = pp.OneOrMore(ppc.number, stop_on=pp.LineEnd())
NUMBERS_ = pp.OneOrMore(ppc.number | NAN, stop_on=pp.LineEnd())
ID = pp.Word(pp.alphanums)
CHAN = pp.Combine(ID + pp.Literal("->") + pp.Opt(ID))
TEMP_LINE = pp.Suppress(pp.Literal(r"P\T")) + NUMBERS
RATE_LINE = ppc.number + NUMBERS_
RATE_LINES = pp.OneOrMore(pp.Group(RATE_LINE))
HIGH_LINE = pp.Suppress(pp.Keyword("O-O")) + NUMBERS_


def parse_output_channel(mess_chan_out: str) -> MessOutputChannelParseResults:
    """Parse MESS output for one channel.

    :param mess_chan_out: MESS output for one channel
    :return: Parse results
    """
    expr = (
        pp.Opt(CHAN)(Key.eq)
        + TEMP_LINE(Key.T)
        + RATE_LINES(Key.data)
        + HIGH_LINE(Key.high)
    )
    res = expr.parse_string(mess_chan_out).as_dict()
    data = res[Key.data]
    eq: str = res[Key.eq]
    id1, id2 = eq.split("->")
    return MessOutputChannelParseResults(
        T=res[Key.T],
        P=[row[0] for row in data],
        k_data=[row[1:] for row in data],
        k_high=res[Key.high],
        id1=id1,
        id2=id2,
    )
