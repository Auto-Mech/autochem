"""Utility functions for parsing Chemkin data."""

import re
from collections import defaultdict
from collections.abc import Sequence

import numpy
import pydantic
import pyparsing as pp
from pyparsing import common as ppc

COMMENT_REGEX = re.compile(r"# .*$|!.*$", flags=re.M)


# Readers
def read_comments(chem_str: str) -> list[str]:
    """Get comments from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Comments
    """
    return [c[1:].strip() for c in re.findall(COMMENT_REGEX, chem_str)]


def read_without_comments(chem_str: str) -> str:
    """Remove comments from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Chemkin string
    """
    return re.sub(COMMENT_REGEX, "", chem_str).strip()


def read_extract_comments(chem_str: str) -> tuple[list[str], str]:
    """Extract comments from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Comments and chemkin string with comments removed.
    """
    return read_comments(chem_str), read_without_comments(chem_str)


def read_equation_reagents(chem_str: str) -> tuple[list[str], list[str]]:
    """Extact equation reagents from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Reactants and products
    """
    res = parse_equation(chem_str)
    return (res.reactants, res.products)


class ChemkinRateParseResults(pydantic.BaseModel):
    reactants: list[str]
    products: list[str]
    reversible: bool
    pressure_dependent: bool
    arrhenius: list[float] = pydantic.Field(default_factory=list)
    efficiencies: dict[str, float] = pydantic.Field(default_factory=dict)
    aux_numbers: dict[str, list[float]] = pydantic.Field(default_factory=dict)
    aux_misc: dict[str, str] = pydantic.Field(default_factory=dict)
    comments: list[str] = pydantic.Field(default_factory=list)


def parse_equation(chem_str: str) -> ChemkinRateParseResults:
    """Extract equation information from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Parse results including reactants, products, reversible,
        pressure_dependent, and efficiencies (only contains third-body, with value 1.0)
    """
    # Split the reaction line off from the auxiliary lines
    chem_str = read_without_comments(chem_str).strip()
    eq, *_ = re.split(r"\n|$|\s[+-]?\d", chem_str, maxsplit=1)
    rct, arrow, prd = re.split(r"(<=>|=>|=)", eq)

    # Assess reversibility
    reversible = arrow != "=>"

    # Extract (+M) or (+X) third body, if present
    efficiencies = {}
    pressure_dependent = False
    rct_res = REAC_SIDE_FALLOFF.parse_string(rct)
    prd_res = REAC_SIDE_FALLOFF.parse_string(prd)
    if Key.falloff in rct_res:
        (third_body,) = rct_res.get(Key.falloff)
        rct = rct_res.get(Key.reagents)
        prd = prd_res.get(Key.reagents)
        pressure_dependent = True
        efficiencies[third_body] = 1.0

    # Determine reactants and products
    rcts, prds = [[s.strip() for s in re.split(r"\+(?!\+)", r)] for r in (rct, prd)]

    # Extract +M third body, if present
    if "M" in rcts or "M" in prds:
        rcts.remove("M")
        prds.remove("M")
        efficiencies["M"] = 1.0

    return ChemkinRateParseResults(
        reactants=rcts,
        products=prds,
        reversible=reversible,
        pressure_dependent=pressure_dependent,
        efficiencies=efficiencies,
    )


def parse_rate(chem_str: str) -> ChemkinRateParseResults:
    """Extract full rate information from a Chemkin string.

    :param chem_str: Chemkin string
    :return: Parse results including reactants, products, reversible,
        pressure_dependent, and efficiencies (only contains third-body, with value 1.0)
    """
    # Extract comments
    comments, chem_str = read_extract_comments(chem_str)

    # Split the reaction line off from the auxiliary lines
    rxn_line, aux_lines = re.split(r"\n|$", chem_str, maxsplit=1)

    # Parse the reaction line
    rxn_res = REAC_LINE.parse_string(rxn_line)

    # Parse the equation
    eq = rxn_res.get(Key.reaction)
    res = parse_equation(eq)

    # Parse the arrhenius parameters
    res.arrhenius = rxn_res.get(Key.arrhenius).as_list()

    # Add comments
    res.comments = comments

    # Parse the auxiliary lines
    res.aux_numbers = defaultdict(list)
    for r, *_ in AUX_ENTRY.scan_string(aux_lines):
        (key, val) = r.as_list()
        if Key.collider in r:
            res.efficiencies[key] = val
        elif Key.numbers in r:
            res.aux_numbers[key].extend(val)
        else:
            res.aux_misc[key] = val

    return res


# Writers
def write_equation(
    reactants: Sequence[str],
    products: Sequence[str],
    reversible: bool = True,
    third_body: str | None = None,
    pressure_dependent: bool = False,
) -> str:
    """Write Chemkin equation to a string.

    :param reactants: Reactants
    :param products: Products
    :param reversible: Whether the reaction is reversible
    :param third_body: Third body
    :param pressure_dependent: Whether the third body influence is pressure-dependent
    :return: Chemkin equation
    """
    reac = " + ".join(reactants)
    prod = " + ".join(products)

    if third_body is not None:
        reac += f" (+ {third_body})" if pressure_dependent else f" + {third_body}"
        prod += f" (+ {third_body})" if pressure_dependent else f" + {third_body}"

    arrow = " = " if reversible else " => "
    return arrow.join([reac, prod])


def write_with_dup(reac_str: str, dup: bool = True) -> str:
    """Format Chemkin reaction string with DUP keyword.

    :param reac_str: Chemkin reaction string
    :param dup: Whether to append the DUP keyword
    :return: Chemkin rate string
    """
    return "\n".join([reac_str, write_aux("DUP")]) if dup else reac_str


def write_efficiencies(
    efficiencies: dict[str, float], third_body: str | None = "M", indent: int = 4
) -> str | None:
    """Format efficiencies as Chemkin string.

    :param efficiencies: Collider efficiencies
    :param third_body: Third body
    :param indent: Indentation, defaults to 4
    :return: Chemkin string
    """
    eff = efficiencies.copy()
    if third_body in eff:
        eff.pop(third_body)

    if not eff:
        return None

    return " " * indent + "  ".join(f"{k} /{v:.3f}/" for k, v in eff.items())


def write_aux(
    key: str,
    val: str | float | Sequence[float] | None = None,
    head_width: int = 55,
    key_width: int = 5,
    indent: int = 4,
    digits: int = 4,
    always_sci: bool = False,
    as_int: bool = False,
) -> str:
    """Format auxiliary Chemkin reaction data.

    :param key: Key, e.g. 'PLOG'
    :param val: Value(s), e.g. '5.000  8500' or ['5.000', '8500']
    :param head_width: Width of the header line, for alignment purposes
    :param key_width: Key column width, defaults to 5
    :param indent: Indentation, defaults to 4
    :return: Chemkin string
    """
    if val is None:
        return " " * indent + key

    val_width = max(head_width - indent - key_width - 2, 0)
    val_str = val
    if isinstance(val, Sequence) and not isinstance(val, str):
        val_str = write_numbers(
            val, digits=digits, always_sci=always_sci, as_int=as_int
        )
    if isinstance(val, float):
        val_str = write_number(val, digits=digits, always_sci=always_sci, as_int=as_int)
    assert isinstance(val_str, str), f"val_str = {val_str}"
    return " " * indent + f"{key:<{key_width}} /{val_str:>{val_width}}/"


# Helpers
def write_numbers(
    nums: Sequence[float],
    digits: int = 4,
    always_sci: bool = False,
    as_int: bool = False,
) -> str:
    """Write a sequence of numbers to a formatted string.

    :param nums: The numbers
    :param digits: How many digits to include, defaults to 4
    :param always_sci: Whether to always use scientific notation; if given as a list,
        this can be used to set scientific notation for individual numbers
    :param as_int: Whether to write intgeger values
    :return: The formatted number sequence string
    """
    return " ".join(
        write_number(n, always_sci=always_sci, digits=digits, as_int=as_int)
        for n in nums
    )


def write_number(
    num: float | int, digits: int = 4, always_sci: bool = False, as_int: bool = False
) -> str:
    """Write a number to a formatted string.

    :param num: The number
    :param digits: How many digitst to include, defaults to 4
    :param always_sci: Whether to always use scientific notation
    :param as_int: Whether to write integer values
    :return: The formatted number string
    """
    # Exact width of scientific notation with 2-digit exponent:
    max_width = digits + 6  # from general formula: digits + 4 + |_log(log(num))_|

    if as_int:
        num = int(num)
        return f"{num:>{max_width}d}"

    exp = int(numpy.floor(numpy.log10(numpy.abs(num)))) if num else 0
    float_width = max(exp + 1, digits + 1) if exp > 0 else numpy.abs(exp) + 1 + digits

    if always_sci or float_width > max_width:
        decimals = digits - 1
        return f"{num:>{max_width}.{decimals}E}"

    decimals = max(0, digits - exp - 1)
    return f"{num:>{max_width}.{decimals}f}"


# Parse helpers
def parenthetical(
    expr: pp.ParserElement, open: str = "(", close: str = ")"
) -> pp.ParserElement:
    """Add slashes on either side of a parse expression.

    :param expr: Expression
    :param open: Open parenthesis
    :param close: Close parenthesis
    :return: Expression
    """
    return pp.Suppress(pp.Literal(open)) + expr + pp.Suppress(pp.Literal(close))


# Helpers
#  - Pyparsing token keys
class Key:
    reaction = "reaction"
    arrhenius = "arrhenius"
    collider = "collider"
    numbers = "numbers"
    misc = "misc"
    reagents = "reagents"
    falloff = "falloff"


#  - Pyparsing expressions
#   - Generic
NUMBER_LIST = pp.Group(pp.OneOrMore(ppc.number))
#   - Reaction line
ARRH_EXPR = pp.WordStart() + ppc.number[3]
REAC_EXPR = pp.SkipTo(ARRH_EXPR).set_parse_action(lambda t: str.rstrip(t[0]))
REAC_LINE = REAC_EXPR(Key.reaction) + ARRH_EXPR(Key.arrhenius)
#   - Auxiliary entry
COLL_KEY = pp.Word(pp.printables, exclude_chars="/")
COLL_VAL = parenthetical(ppc.number, "/", "/")
NUMS_KEY = pp.Word(pp.alphanums)
NUMS_VAL = parenthetical(NUMBER_LIST, "/", "/")
MISC_KEY = pp.Word(pp.alphanums)
MISC_VAL = parenthetical(pp.SkipTo("/"), "/", "/")
COLL_ENTRY = COLL_KEY + COLL_VAL
NUMS_ENTRY = NUMS_KEY + NUMS_VAL
MISC_ENTRY = MISC_KEY + MISC_VAL
AUX_ENTRY = COLL_ENTRY(Key.collider) ^ NUMS_ENTRY(Key.numbers) ^ MISC_ENTRY(Key.misc)
#   - Reaction side with third body
PLUS = pp.Suppress(pp.Literal("+"))
FALLOFF = parenthetical(PLUS + pp.Word(pp.alphanums))
REAC_SIDE_FALLOFF = pp.SkipTo(FALLOFF ^ pp.StringEnd())(Key.reagents) + pp.Opt(FALLOFF)(
    Key.falloff
)
