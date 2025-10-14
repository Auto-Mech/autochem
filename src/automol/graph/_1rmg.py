"""Conversion to and from RMG adjacency lists
"""

import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from .base import from_data

BOND_ORDER_DCT = {"S": 1, "D": 2, "T": 3, "B": 1.5}

KEY = ppc.integer("key")
LABEL = pp.Suppress("*") + ppc.integer("label")
SYMBOL = pp.Word(pp.alphas, max=2)("symbol")
UNPAIRED = pp.Suppress("u") + ppc.integer("u")
PAIRED = pp.Suppress("p") + ppc.integer("p")
CHARGE = pp.Suppress("c") + ppc.signed_integer("c")
BOND_ORDER = pp.Char(BOND_ORDER_DCT.keys())("bond_order")
BOND = pp.Suppress("{") + KEY + pp.Suppress(",") + BOND_ORDER + pp.Suppress("}")
BONDS = pp.OneOrMore(pp.Group(BOND))("bonds")

ADJ_LINE = (
    KEY
    + pp.Opt(LABEL)
    + SYMBOL
    + UNPAIRED
    + pp.Opt(PAIRED)
    + pp.Opt(CHARGE)
    + pp.Opt(BONDS)
)
RMG_ADJACENCY_LIST = pp.OneOrMore(pp.Group(ADJ_LINE))("adj_list")


def from_rmg_adjacency_list(adj_str):
    """Get a graph from an RMG adjacency list

    :param adj_str: An adjacency list string
    :return: The graph
    """
    adj_par_ret = RMG_ADJACENCY_LIST.parseString(adj_str).asDict()["adj_list"]
    return from_parsed_rmg_adjacency_list(adj_par_ret)


def from_parsed_rmg_adjacency_list(adj_par_ret: dict):
    """Get a graph from the parse result of an RMG adjacency list

    :param adj_par_ret: A parse result from the RMG_ADJACENCY_LIST parser
    :return: The graph
    """
    symb_dct = {}
    ord_dct = {}

    for adj_par_row in adj_par_ret:
        key = adj_par_row["key"]
        symb = adj_par_row["symbol"]
        bonds = adj_par_row.get("bonds", ())

        symb_dct[key] = symb
        for bond in bonds:
            nkey = bond["key"]

            order = BOND_ORDER_DCT[bond["bond_order"]]
            ord_dct[frozenset({key, nkey})] = order

    gra = from_data(symb_dct, ord_dct.keys(), bnd_ord_dct=ord_dct)
    return gra
