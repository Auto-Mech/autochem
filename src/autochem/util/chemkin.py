"""Utility functions for reading and writing Chemkin data."""

import re
from collections import defaultdict
from collections.abc import Callable, Sequence

import more_itertools as mit
import numpy as np
import pydantic
import pyparsing as pp
from pyparsing import common as ppc

COMMENT_REGEX = re.compile(r"# .*$|!.*$", flags=re.MULTILINE)


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


# Thermo
class ChemkinThermoParseResults(pydantic.BaseModel):
    """Chemkin thermo parse results."""

    name: str
    formula: dict[str, int]
    T_min: float
    T_max: float
    coeffs: list[float]
    phase: str = "G"
    T_mid: float | None = None
    date: str | None = None
    comments: list[str] = pydantic.Field(default_factory=list)


def parse_thermo(therm_str: str) -> ChemkinThermoParseResults:
    """Extract all thermo information from a Chemkin thermo string.

    :param therm_str: Chemkin string
    :return: Parse results including name, formula, coefficients, etc.
    """
    comments, therm_str = read_extract_comments(therm_str)
    line1, _, coeff_str = therm_str.strip().partition("\n")

    name = line1[:18].strip()
    date = line1[18:24].strip()
    form_str = line1[24:44].strip() + line1[73:78].strip()
    form_dct = dict(FORM_ENTRIES.parse_string(form_str).as_list())
    phase = line1[44]
    temp_expr = THERM_TEMP(Key.min) + THERM_TEMP(Key.max) + pp.Opt(THERM_TEMP)(Key.mid)
    temp_dct = temp_expr.parse_string(line1[45:73]).as_dict()
    coeffs = COEFFS.parse_string(coeff_str).as_list()

    return ChemkinThermoParseResults(
        name=name,
        formula=form_dct,
        T_min=temp_dct.get(Key.min),
        T_max=temp_dct.get(Key.max),
        coeffs=coeffs,
        phase=phase,
        T_mid=temp_dct.get(Key.mid),
        date=date,
        comments=comments,
    )


# Rates
class ChemkinRateParseResults(pydantic.BaseModel):
    """Chemkin rate parse results."""

    reactants: list[str]
    products: list[str]
    reversible: bool
    pressure_dependent: bool
    arrhenius: list[float] = pydantic.Field(default_factory=list)
    efficiencies: dict[str, float] = pydantic.Field(default_factory=dict)
    aux_numbers: dict[str, list[float]] = pydantic.Field(default_factory=dict)
    aux_misc: dict[str, str] = pydantic.Field(default_factory=dict)
    comments: list[str] = pydantic.Field(default_factory=list)


def parse_equation(reac_str: str) -> ChemkinRateParseResults:
    """Extract equation information from a Chemkin reaction string.

    :param chem_str: Chemkin string
    :return: Parse results including reactants, products, reversible,
        pressure_dependent, and efficiencies (only contains third-body, with value 1.0)
    """
    # Split the reaction line off from the auxiliary lines
    reac_str = read_without_comments(reac_str).strip()
    eq, *_ = re.split(r"\n|$|\s[+-]?\d", reac_str, maxsplit=1)
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


def parse_rate(rate_str: str) -> ChemkinRateParseResults:
    """Extract all rate information from a Chemkin rate string.

    :param chem_str: Chemkin string
    :return: Parse results including reactants, products, reversible,
        pressure_dependent, and efficiencies (only contains third-body, with value 1.0)
    """
    # Extract comments
    comments, rate_str = read_extract_comments(rate_str)

    # Split the reaction line off from the auxiliary lines
    rxn_line, aux_lines = re.split(r"\n|$", rate_str, maxsplit=1)

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
    *,
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


def write_with_dup(reac_str: str, *, dup: bool = True) -> str:
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


def write_aux(  # noqa: PLR0913
    key: str,
    val: str | float | Sequence[float] | None = None,
    *,
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


def write_therm_entry_header(  # noqa: PLR0913
    name: str,
    form_dct: dict[str, int],
    T_min: float,  # noqa: N803
    T_max: float,  # noqa: N803
    *,
    T_mid: float | None = None,  # noqa: N803
    charge: int = 0,
    phase: str = "G",
    date: str = "",
) -> str:
    """Write the header line for a Chemkin thermo entry.

    :param name: Name of the species
    :param form_dct: Dictionary of the species formula
    :param T_min: Minimum temperature
    :param T_max: Maximum temperature
    :param T_mid: Mid temperature, optional
    :param charge: Charge, defaults to 0
    :param date: Date, defaults to empty string
    :return: Chemkin thermo entry header line
    """
    # Build formula strings
    assert len(form_dct) <= 5, f"Too many elements for Chemkin: |{form_dct}| > 5"
    form_items = [(k, form_dct.get(k)) for k in sorted(form_dct, key=symbol_sort_key())]
    if charge:
        form_items.append(("E", charge))
    form_items1 = form_items[:4]
    form_items2 = form_items[4:]
    form_str1 = "".join(f"{k: <2}{v: >3}" for k, v in form_items1)
    form_str2 = "".join(f"{k: <2}{v: >3}" for k, v in form_items2)

    # Build temperature string
    temp_str = f"{T_min:>10.1f}{T_max:>10.1f}"
    temp_str += f"{T_mid:>8.1f}" if T_mid else ""

    return (
        f"{name: <18}{date: <6}{form_str1: <20}{phase:<1}{temp_str: <28}{form_str2: <5}"
    )


def write_therm_entry_coefficient_lines(coeffs: Sequence[float]) -> list[str]:
    """Write the coefficient lines for a Chemkin thermo entry.

    :param coeffs: Coefficients
    :return: Chemkin thermo entry coefficient lines
    """
    num_fmt = "{:>15.8E}"
    return ["".join(map(num_fmt.format, cs)) for cs in mit.chunked(coeffs, 5)]


def symbol_sort_key(
    symbs_first: Sequence[str] = ("C", "H"), symbs_last: Sequence[str] = ()
) -> Callable[[str], tuple[int, str]]:
    """Sort key for atomic symbols.

    :param symbs_first: Symbols to sort first
    :return: Sort ke
    """
    symbs_first = list(map(str.title, symbs_first))
    symbs_last = list(map(str.title, symbs_last))

    def _key(symb: str) -> tuple[int, str]:
        """Symbol sort key."""
        symb = symb.title()
        if symb in symbs_first:
            val = symbs_first.index(symb)
        elif symb in symbs_last:
            val = len(symbs_first) + symbs_last.index(symb) + 1
        else:
            val = len(symbs_first)
        return val, symb

    return _key


# Helpers
def write_numbers(
    nums: Sequence[float],
    *,
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
    num: float, *, digits: int = 4, always_sci: bool = False, as_int: bool = False
) -> str:
    """Write a number to a formatted string.

    :param num: The number
    :param digits: How many digits to include, defaults to 4
    :param always_sci: Whether to always use scientific notation
    :param as_int: Whether to write integer values
    :return: The formatted number string
    """
    # Exact width of scientific notation with 2-digit exponent:
    max_width = digits + 6  # from general formula: digits + 4 + |_log(log(num))_|

    if as_int:
        num = int(num)
        return f"{num:>{max_width}d}"

    exp = int(np.floor(np.log10(np.abs(num)))) if num else 0
    float_width = max(exp + 1, digits + 1) if exp > 0 else np.abs(exp) + 1 + digits

    if always_sci or float_width > max_width:
        decimals = digits - 1
        return f"{num:>{max_width}.{decimals}E}"

    decimals = max(0, digits - exp - 1)
    return f"{num:>{max_width}.{decimals}f}"


# Parse helpers
def parenthetical(
    expr: pp.ParserElement, open_: str = "(", close: str = ")"
) -> pp.ParserElement:
    """Add slashes on either side of a parse expression.

    :param expr: Expression
    :param open: Open parenthesis
    :param close: Close parenthesis
    :return: Expression
    """
    return pp.Suppress(pp.Literal(open_)) + expr + pp.Suppress(pp.Literal(close))


# Helpers
class Key:
    """Pyparsing token keys."""

    reaction = "reaction"
    arrhenius = "arrhenius"
    collider = "collider"
    numbers = "numbers"
    misc = "misc"
    reagents = "reagents"
    falloff = "falloff"
    mid = "mid"
    min = "min"
    max = "max"


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
    Key.falloff,
)
#   - Thermo entry formula
FORM_KEY = pp.Word(pp.alphas, max=2)
FORM_VAL = ppc.signed_integer
FORM_ENTRY = pp.Group(FORM_KEY + FORM_VAL)
FORM_ENTRIES = pp.OneOrMore(FORM_ENTRY)
#   - Thermo entry temperatures
THERM_TEMP = ppc.number
#   - Coefficient numbers
LINE_NUM = ppc.integer
COEFF = ppc.sci_real
COEFFS = pp.OneOrMore(pp.OneOrMore(COEFF) + pp.Suppress(LINE_NUM))
