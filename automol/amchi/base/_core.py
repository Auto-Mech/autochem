""" Level 2 AMChI functions (depend on L1)

The parsing functions apply equally well to InChI or AMChI strings, so the
documentation simply refers to "ChI" strings.
"""

import functools
import itertools
import warnings
from collections import abc

import numpy
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from ... import form, util
from ...util import dict_


# Build parser
def layer_parser(key=pp.Empty(), chars=pp.printables, content_only=False):
    """Create a parser for an InChI layer"""
    layer_prefix = pp.Combine(pp.Suppress("/") + key)
    layer_content = pp.Word(chars, excludeChars="/")
    if content_only:
        layer = pp.Combine(pp.Suppress(layer_prefix) + layer_content)
    else:
        layer = pp.Group(layer_prefix + layer_content)
    return layer


PREFIX = pp.Literal("InChI") ^ pp.Literal("AMChI")
VERSION = pp.Combine(pp.Char(pp.nums) + pp.Opt("S"))
FORMULA = layer_parser(content_only=True)
CONNECTIVITY = pp.Opt(layer_parser(key="c"))
HYDROGENS = pp.Opt(layer_parser(key="h"))
CHARGE = pp.Opt(layer_parser(key="q"))
PROTON = pp.Opt(layer_parser(key="p"))
BOND_STEREO = pp.Opt(layer_parser(key="b"))
ATOM_STEREO = pp.Opt(layer_parser(key="t"))
MIRROR_IMAGE = pp.Opt(layer_parser(key="m"))
STEREO_CONFIGURATION = pp.Opt(layer_parser(key="s"))
STEREO_LAYERS = BOND_STEREO + ATOM_STEREO + MIRROR_IMAGE + STEREO_CONFIGURATION
ISOTOPE = layer_parser(key="i")
BOND_BREAKING = pp.Opt(layer_parser(key="k"))
BOND_FORMING = pp.Opt(layer_parser(key="f"))
REVERSE_DIRECTION = pp.Opt(layer_parser(key="r"))

CHI_PARSER = (
    PREFIX("prefix")
    + pp.Suppress("=")
    + VERSION("version")
    + FORMULA("formula")
    + pp.Group(CONNECTIVITY + HYDROGENS)("main")
    + pp.Group(CHARGE + PROTON)("charge")
    + pp.Group(STEREO_LAYERS)("stereo")
    + pp.Group(pp.Opt(ISOTOPE + HYDROGENS + STEREO_LAYERS))("isotope")
    + pp.Group(BOND_BREAKING + BOND_FORMING + REVERSE_DIRECTION)("ts")
)

MAIN_PFXS = ("c", "h")
CHAR_PFXS = ("q", "p")
STE_PFXS = ("b", "t", "m", "s")
ISO_NONSTE_PFXS = ("i", "h")
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS
TS_PFXS = ("k", "f", "r")


# # constructor
def from_data(
    fml_lyr,
    main_lyr_dct=None,
    char_lyr_dct=None,
    ste_lyr_dct=None,
    iso_lyr_dct=None,
    ts_lyr_dct=None,
    pfx="AMChI",
):
    """Build a ChI string from layers

    :param fml_lyr: The formula layer
    :type fml_lyr: str
    :param main_lyr_dct: main layers, specifying connectivity and implicit
        hydrogens, by key ('c' and 'h')
    :type main_lyr_dct: dict[str: str]
    :param char_lyr_dct: charge layers, by key ('q' and 'p')
    :type char_lyr_dct: dict[str: str]
    :param ste_lyr_dct: stero layers, by key ('b', 't', 'm', and 's')
    :type ste_lyr_dct: dict[str: str]
    :param iso_lyr_dct: isotope layers, by key ('i', 'h', 'b', 't', 'm',
        and 's')
    :type iso_lyr_dct: dict[str: str]
    :param ts_lyr_dct: TS layers, by key ('k', 'f', and 'r')
    :type ts_lyr_dct: dict[str: str]
    :param pfx: The prefix to use ("AMChI" or "InChI")
    :type pfx: str, optional
    :rtype: str
    """
    assert pfx in ("AMChI", "InChI")

    main_dct = dict_.empty_if_none(main_lyr_dct)
    char_dct = dict_.empty_if_none(char_lyr_dct)
    ste_dct = dict_.empty_if_none(ste_lyr_dct)
    iso_dct = dict_.empty_if_none(iso_lyr_dct)
    ts_dct = dict_.empty_if_none(ts_lyr_dct)

    main_lyrs = [
        pfx + lyr
        for pfx, lyr in zip(MAIN_PFXS, dict_.values_by_key(main_dct, MAIN_PFXS))
        if lyr
    ]
    char_lyrs = [
        pfx + lyr
        for pfx, lyr in zip(CHAR_PFXS, dict_.values_by_key(char_dct, CHAR_PFXS))
        if lyr
    ]
    ste_lyrs = [
        pfx + lyr
        for pfx, lyr in zip(STE_PFXS, dict_.values_by_key(ste_dct, STE_PFXS))
        if lyr
    ]
    iso_lyrs = [
        pfx + lyr
        for pfx, lyr in zip(ISO_PFXS, dict_.values_by_key(iso_dct, ISO_PFXS))
        if lyr
    ]
    ts_lyrs = [
        pfx + lyr
        for pfx, lyr in zip(TS_PFXS, dict_.values_by_key(ts_dct, TS_PFXS))
        if lyr
    ]

    start = f"{pfx}=1"
    chi = "/".join(
        [start, fml_lyr] + main_lyrs + char_lyrs + ste_lyrs + iso_lyrs + ts_lyrs
    )

    return chi


# # recalculate/standardize
def standard_form(
    chi, stereo=True, racem=False, racem_m_layer=True, ste_dct=None, iso_dct=None
):
    """Return a ChI string in standard form.

    Includes only the standard layers.

    :param chi: ChI string
    :type chi: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param racem: parameter to designate a racemic mixture, if chiral
    :type racem: bool
    :param racem_m: Include the m layer for a racemic ChI?, default True
    :type racem_m: bool, optional
    :param ste_dct: a dictionary to overwrite stereo information; layers
        not overwritten will be left in tact; if the attempted overwrite
        fails, the function will return None
    :type ste_dct: dict
    :param iso_dct: a dictionary to overwrite isotope stereo information;
        layers not overwritten will be left intact; if the attempted
        overwrite fails, the function will return None
    :type iso_dct: dict
    :rtype: str
    """
    pfx = prefix(chi)
    fml_slyr = formula_layer(chi)
    main_dct = main_layers(chi)
    char_dct = charge_layers(chi)
    ts_dct = ts_layers(chi)

    extra_ste_dct = ste_dct
    extra_iso_dct = iso_dct

    if stereo:
        ste_dct = stereo_layers(chi)
        iso_dct = isotope_layers(chi)
    else:
        ste_dct = {}
        iso_dct = dict_.by_key(isotope_layers(chi), ISO_NONSTE_PFXS)

    if extra_ste_dct is not None:
        ste_dct.update(extra_ste_dct)

    if extra_iso_dct is not None:
        iso_dct.update(extra_iso_dct)

    if racem:
        if "m" in ste_dct:
            if racem_m_layer:
                ste_dct["m"] = ste_dct["m"].translate(str.maketrans("01", "22"))
            else:
                ste_dct.pop("m")

            ste_dct["s"] = "3"

        if "m" in iso_dct:
            if racem_m_layer:
                ste_dct["m"] = ste_dct["m"].translate(str.maketrans("01", "22"))
            else:
                ste_dct.pop("m")

            iso_dct["s"] = "3"

    chi = from_data(
        fml_slyr,
        main_lyr_dct=main_dct,
        char_lyr_dct=char_dct,
        ste_lyr_dct=ste_dct,
        iso_lyr_dct=iso_dct,
        ts_lyr_dct=ts_dct,
        pfx=pfx,
    )

    # Make sure multi-component ChIs are properly sorted
    if has_multiple_components(chi):
        chi = join(split(chi))

    return chi


# # getters
def prefix(chi):
    """Determine the chemical identifier prefix (InChI or AMChI).

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return chi_dct["prefix"]


def version(chi):
    """Determine version number of the ChI string.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return chi_dct["version"]


def formula_layer(chi):
    """Parse the ChI string for the formula layer

    :param chi: ChI string
    :type chi: str
    :returns: the formula string
    :rtype: str
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return chi_dct["formula"]


def main_layers(chi):
    """Parse the ChI string for the main layers ('c' and 'h').

    :param chi: ChI string
    :type chi: str
    :returns: the main layers, as a dictionary keyed by layer prefixes
    :rtype: dict[str: str]
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return dict(chi_dct["main"])


def charge_layers(chi):
    """Parse the ChI string for the charge layers ('q' and 'p').

    :param chi: ChI string
    :type chi: str
    :returns: the charge layers, as a dictionary keyed by layer prefixes
    :rtype: dict[str: str]
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return dict(chi_dct["charge"])


def stereo_layers(chi):
    """Parse the ChI string for the stereo layers ('b', 't', 'm', and 's')

    :param chi: ChI string
    :type chi: str
    :returns: the stereo layers, as a dictionary keyed by layer prefixes
    :rtype: dict[str: str]
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return dict(chi_dct["stereo"])


def isotope_layers(chi):
    """Parse the ChI string for the isotope layers ('i', 'h', 'b', 't', 'm',
    and 's')

    :param chi: ChI string
    :type chi: str
    :returns: the isotope layers, as a dictionary keyed by layer prefixes
    :rtype: dict[str: str]
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return dict(chi_dct["isotope"])


def ts_layers(chi):
    """Parse the ChI string for the TS layers ('k', 'f', and 'r')

    :param chi: ChI string
    :type chi: str
    :returns: the TS layers, as a dictionary keyed by layer prefixes
    :rtype: dict[str: str]
    """
    chi_dct = CHI_PARSER.parseString(chi).asDict()
    return dict(chi_dct["ts"])


# # setters
def with_inchi_prefix(chi):
    """Return a ChI with InChI prefix, whether AMChI or InChI.

    :param chi: ChI string
    :type chi: str
    :returns: InChI string
    :rtype: str
    """
    pfx = prefix(chi)
    if pfx == "AMChI":
        chi = "InChI" + chi[5:]
    else:
        assert pfx == "InChI", f"ChI string '{chi}' has unknown prefix '{pfx}'."
    return chi


def reflect(chi):
    """If this is an enantiomer, flip to the other enantiomer by changing the
    m-layer

    :param chi: InChI string
    :type chi: str
    :returns: the other enantiomer
    :rtype: bool
    """
    if is_enantiomer(chi):
        ste_upd_dct = stereo_layers(chi)
        iso_upd_dct = isotope_layers(chi)

        refl_trans = str.maketrans("01", "10")
        if "m" in ste_upd_dct:
            ste_upd_dct["m"] = ste_upd_dct["m"].translate(refl_trans)

        if "m" in iso_upd_dct:
            iso_upd_dct["m"] = iso_upd_dct["m"].translate(refl_trans)

        chi = standard_form(chi, ste_dct=ste_upd_dct, iso_dct=iso_upd_dct)

    return chi


def canonical_enantiomer(chi):
    """Convert this ChI string to a canonical enantiomer, if it isn't already

    Works for multi-component InChIs or lists of InChIs

    :param chi: ChI string
    :type chi: str
    :returns: The reflected ChI string, if it was a non-canonical
        enantiomer; otherwise, it returns the original ChI string
    """
    if not is_canonical_enantiomer(chi):
        chi = reflect(chi) if isinstance(chi, str) else tuple(map(reflect, chi))
    return chi


def reflect_reaction(rct_chis, prd_chis):
    """Apply a reflection operation to a reaction.

    :param rct_chis: A list of ChIs for the reactants
    :type rct_chis: list[str]
    :param prd_chis: A list of ChIs for the products
    :type prd_chis: list[str]
    :returns: The reactant and product ChIs, all reflected
    :rtype: (list[str], list[str])
    """
    rct_chis = tuple(map(reflect, rct_chis))
    prd_chis = tuple(map(reflect, prd_chis))
    return (rct_chis, prd_chis)


def canonical_enantiomer_reaction(rct_chi, prd_chi):
    """Convert this reaction into a canonical combination of enantiomers

    :param rct_chi: A multi-component ChI or list of ChIs for the reactants
    :type rct_chi: str or list[str]
    :param prd_chi: A multi-component ChI or list of ChIs for the products
    :type prd_chi: str or list[str]

    :returns: The reflected reaction, if it was a non-canonical combination
        of enantiomers; otherwise, it returns the original reaction.
    """
    if not is_canonical_enantiomer_reaction(rct_chi, prd_chi):
        rct_chi = (
            reflect(rct_chi)
            if isinstance(rct_chi, str)
            else tuple(map(reflect, rct_chi))
        )
        prd_chi = (
            reflect(prd_chi)
            if isinstance(prd_chi, str)
            else tuple(map(reflect, prd_chi))
        )
    return rct_chi, prd_chi


# # conversions
def formula(chi):
    """Generate a formula dictionary from a ChI string.

    :param chi: ChI string
    :type chi: str
    :rtype: dict[str: int]
    """
    # split it up to handle hard-coded molecules in multi-component chis
    chis = split(chi)
    fml_strs = list(map(formula_layer, chis))
    fmls = list(map(form.from_string, fml_strs))
    fml = functools.reduce(form.join, fmls)
    return fml


def formula_string(chi):
    """Generate a formula string from a ChI string.

    For multi-component ChIs, this will be the *combined* formula.

    Use `formula_layer(chi)` if you want the multi-component layer.

    :param chi: ChI string
    :type chi: str
    :rtype: str
    """
    fml = formula(chi)
    fml_str = form.string(fml)
    return fml_str


def without_stereo(chi):
    """Remove all stereo layers"""

    return standard_form(chi, stereo=False)


def racemic(chi):
    """If chiral, convert the ChI into a racemic mixture

    This drops the /m layer and replaces /s1 with /s3, indicating a racemic
    mixture. The chirality of the species is still implied by the presence
    of the /s layer.

    :param chi: ChI string
    :type chi: str
    """
    return standard_form(chi, racem=True)


def connectivity(chi, parse_connection_layer=True, parse_h_layer=True):
    """Return the 'c' and 'h' layers of the connectivity string

    The user may also specify what combination of the two layers
    that they wish to return
    """

    # Read the two sublayers that are requested to be parsed
    conn_lyrs = main_layers(chi)

    if parse_connection_layer:
        cslyr = conn_lyrs.get("c", "")
    else:
        cslyr = ""

    if parse_h_layer:
        hslyr = conn_lyrs.get("h", "")
    else:
        hslyr = ""

    # Write the parts of the connectivity string based on what was parsed
    if cslyr and hslyr:
        _str = f"c{cslyr}/h{hslyr}"
    elif cslyr:
        _str = f"c{cslyr}"
    elif hslyr:
        _str = f"h{hslyr}"
    else:
        _str = None

    return _str


def are_enantiomers(chi_a, chi_b):
    """Are these InChI enantiomers of eachother?

    :param chi: InChI string
    :type chi: str
    """
    ste_dct_a = stereo_layers(chi_a)
    ste_dct_b = stereo_layers(chi_b)
    enant = False
    if main_layers(chi_a) == main_layers(chi_b):
        if len(ste_dct_b.keys()) == len(ste_dct_a.keys()) and "m" in ste_dct_a.keys():
            if ste_dct_a["m"] != ste_dct_b["m"]:
                if "t" in ste_dct_a.keys():
                    if ste_dct_a["t"] == ste_dct_b["t"]:
                        enant = True
                else:
                    enant = True
    return enant


def are_diastereomers(chi_a, chi_b):
    """Are these InChI diastereomers of each other?

    Checks if main layer is the same, if so then checks
    if the stereo layers differ in any way.

    :returns: whether or not the InChI is enantiomeric
    :rtype: bool
    """

    diast = False
    if chi_a != chi_b:  # chk not same InChIs
        if main_layers(chi_a) == main_layers(chi_b):
            ste_dct_a = stereo_layers(chi_a)
            ste_dct_b = stereo_layers(chi_b)
            # b-lyr are diastereomers; t-lyr may be, need check
            if len(ste_dct_a.keys()) == len(ste_dct_b.keys()):
                if "b" in ste_dct_a.keys():
                    if ste_dct_a["b"] != ste_dct_b["b"]:
                        diast = True
                elif "t" in ste_dct_a.keys():
                    if ste_dct_a["t"] != ste_dct_b["t"]:
                        diast = True

    return diast


# # properties
# # # prefix
def is_amchi(chi: str) -> bool:
    """Is this an AMChI string?

    :param chi: ChI string
    :type chi: str
    :return: `True` if it is, `False`, if it isn't
    :rtype: bool
    """
    return prefix(chi) == "AMChI"


def is_inchi(chi: str) -> bool:
    """Is this an InChI string?

    :param chi: ChI string
    :type chi: str
    :return: `True` if it is, `False`, if it isn't
    :rtype: bool
    """
    return prefix(chi) == "InChI"


# # # formula layer
def symbols(chi, one_indexed=False):
    """Determine the atomic symbols of backbone atoms in a ChI string

    :param chi: ChI string
    :type chi: str
    :param one_indexed: use one-indexing?
    :type one_indexed: bool
    :returns: a dictionary of atomic symbols, keyed by canonical index
    :rtype: dict[int: str]
    """
    fml = formula(chi)
    pool = list(form.sorted_symbols(fml.keys(), symbs_first=["C"]))

    # If there are only hydrogens, then one of them must be a backbone atom
    if set(pool) == {"H"}:
        symbs = ["H"]
    # Otherwise, remove all hydrogens from the list of backbone atom symbols
    else:
        symbs = [
            s for symb in pool for s in itertools.repeat(symb, fml[symb]) if symb != "H"
        ]

    shift = 1 if one_indexed else 0
    symb_dct = {k + shift: s for k, s in enumerate(symbs)}
    return symb_dct


def canonical_indices(chi, one_indexed=False):
    """Determine the list of canonical indices for a ChI string"""
    idxs = sorted(symbols(chi, one_indexed=one_indexed).keys())
    return idxs


# # # main layers
def bonds(chi, one_indexed=False):
    """Determine bonds between backbone atoms in a ChI string

    :param chi: ChI string
    :type chi: str
    :param one_indexed: use one-indexing?
    :type one_indexed: bool
    """
    # Set up the pyparsing parser
    chain = pp.delimitedList(ppc.integer, "-")
    chains = chain + pp.ZeroOrMore("," + chain)
    side_chain = pp.nestedExpr("(", ")", chains)
    parser = pp.Opt(chain + pp.ZeroOrMore(side_chain + chain))

    # Do the parsing. This produces a nested list of numbers and commas
    # mirroring the connection layer
    main_lyr_dct = main_layers(chi)
    conn_lyr = main_lyr_dct["c"] if "c" in main_lyr_dct else ""
    conn_lst = parser.parseString(conn_lyr).asList()

    shift = 0 if one_indexed else -1

    def _recurse_find_bonds(bnds, conn_lst):
        # Pop the current idx
        if conn_lst:
            idx = conn_lst.pop(0) + shift

        # If there are elements left, continue
        if conn_lst:
            # Look at the next element
            obj = conn_lst[0]

            # Deal with the case where obj is a sequence
            if isinstance(obj, abc.Sequence):
                # In this case, we have multiple branches

                # Pop the sequence
                obj = conn_lst.pop(0)

                # Split the sequence at commas
                lsts = util.breakby(obj, ",")

                # Add bonds to the first element and continue the recursion for
                # each sub list from the split
                for lst in map(list, lsts):
                    nei = lst[0] + shift
                    bnds.add(frozenset({idx, nei}))

                    _recurse_find_bonds(bnds, lst)

                # Now that the list has been dealt with, continue with the
                # element following it, which is also bonded to `idx`
                nei = conn_lst[0] + shift

                # Check that this is an integer (it should always be)
                assert isinstance(
                    nei, int
                ), f"Something is wrong. {nei} should be an integer."

                # Add the bond
                bnds.add(frozenset({idx, nei}))

                # Continue the recursion
                bnds = _recurse_find_bonds(bnds, conn_lst)
            # Deal with the case where obj is a number
            else:
                # In this case, we are continuing along a chain

                # Add the bond
                nei = obj + shift
                bnds.add(frozenset({idx, nei}))

                # Continue the recursion
                bnds = _recurse_find_bonds(bnds, conn_lst)

        return bnds

    bnds = _recurse_find_bonds(set(), conn_lst)

    return bnds


def adjacency_list(chi):
    """Generate an adjacency list for backbone atoms in a ChI string

    :param chi: ChI string
    :type chi: str
    :param one_indexed: use one-indexing?
    :type one_indexed: bool
    """
    idxs = canonical_indices(chi, one_indexed=False)
    bnds = bonds(chi, one_indexed=False)
    adjs = [sorted(next(iter(b - {i})) for b in bnds if i in b) for i in idxs]
    return adjs


def hydrogen_valences(chi, one_indexed=False):
    """Determine the hydrogen valences of backbone atoms in a ChI string

    :param chi: ChI string
    :type chi: str
    :param one_indexed: use one-indexing?
    :type one_indexed: bool
    :returns: a dictionary of hydrogen valences, keyed by canonical index
    :rtype: dict[int: int]
    """
    # Set up the parser
    integer = ppc.integer
    sep = "-" | pp.Suppress(",")
    fixedh = integer + pp.ZeroOrMore(sep + integer) + "H" + pp.Opt(integer)
    mobileh = pp.Combine(
        pp.Suppress("(")
        + "H"
        + pp.Opt(integer)
        + pp.OneOrMore("," + integer)
        + pp.Suppress(")")
    )
    parser = pp.Opt(
        pp.Opt(pp.Group(fixedh) + pp.ZeroOrMore(sep + pp.Group(fixedh)))
        + pp.Opt(sep)
        + pp.ZeroOrMore(mobileh)
    )

    mobileh_parser = (
        pp.Combine("H" + pp.Opt(integer))
        + sep
        + pp.Group(pp.delimitedList(integer, ","))
    )

    # Do the parsing
    main_lyr_dct = main_layers(chi)
    nhyd_lyr = main_lyr_dct["h"] if "h" in main_lyr_dct else ""
    nhyd_lsts = tuple(parser.parseString(nhyd_lyr).asList())

    # Interpret the list
    shift = 0 if one_indexed else -1
    all_idxs = canonical_indices(chi, one_indexed=one_indexed)
    nhyd_dct = dict_.by_key({}, all_idxs, fill_val=0)
    mob_lsts = []
    for nhyd_lst in nhyd_lsts:
        # Mobile hydrogen blocks will be strings instead of lists
        if isinstance(nhyd_lst, str) and nhyd_lst.startswith("H"):
            mob_lsts.append(nhyd_lst)
        # Otherwise, this is not a mobile hydrogen block -- proceed as usual
        else:
            if isinstance(nhyd_lst[-1], int):
                nhyd = nhyd_lst[-1]
                nhyd_lst = nhyd_lst[:-2]
            else:
                nhyd = 1
                nhyd_lst = nhyd_lst[:-1]

            lsts = list(map(list, util.breakby(nhyd_lst, "-")))
            idxs = lsts.pop(0)
            for lst in lsts:
                idxs.extend(range(idxs[-1] + 1, lst[0]))
                idxs.extend(lst)
            idxs = [k + shift for k in idxs]
            nhyd_dct.update({k: nhyd for k in idxs})

    # Add in mobile hydrogens after we get the others
    for mob_lst in mob_lsts:
        mob_lst = tuple(mobileh_parser.parseString(nhyd_lst).asList())
        nmob_str = mob_lst[0][1:]
        nmob = int(nmob_str) if nmob_str else 1
        # Add available mobile hydrogens to the first nmob atoms
        idxs = [k + shift for k in mob_lst[1]][:nmob]
        nhyd_dct.update({k: nhyd_dct[k] + 1 for k in idxs})

    return nhyd_dct


# # # charge layers
def charge(chi):
    """Determine charge from the ChI string

    :param chi: ChI string
    :type chi: str
    :rtype: int
    """
    char_lyr_dct = charge_layers(chi)
    char = int(char_lyr_dct["q"]) if "q" in char_lyr_dct else 0
    return char


# # # stereo layers
def bond_stereo_parities(chi, one_indexed=False, ordered_key=False):
    """Parse the bond stereo parities from the stereochemistry layers.

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :param ordered_key: Return ordered or unordered bond keys?
    :returns: A dictionary mapping bond keys onto parities
    :rtype: dict[frozenset[int]: bool]
    """
    lyr_dct = stereo_layers(chi)
    if "b" not in lyr_dct:
        return {}

    ste_dct = _bonds(lyr_dct["b"], one_indexed=one_indexed, ordered_key=ordered_key)
    return ste_dct


def atom_stereo_parities(chi, one_indexed=False):
    """Parse the atom stereo parities from the stereochemistry layers.

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :returns: A dictionary mapping atom keys onto parities
    :rtype: dict[int: bool]
    """
    ste_lyr_dct = stereo_layers(chi)
    atm_ste_dct = _atom_stereo_parities(ste_lyr_dct, one_indexed=one_indexed)
    return atm_ste_dct


def is_inverted_enantiomer(chi):
    """Determine enantiomer inversion from the stereo layers.

    :param chi: ChI string
    :type chi: str
    :returns: whether or not the ChI is inverted; returns None if not an
        enantiomer
    :rtype: bool
    """
    ste_lyr_dct = stereo_layers(chi)
    is_inv = _is_inverted_enantiomer(ste_lyr_dct)
    return is_inv


def is_canonical_enantiomer(chi):
    """Is this a canonical enantiomer? Also returns true for non-enantiomers.

    For multi-component InChIs, it checks that the first enantiomeric
    component (if any), in sorted order, is canonical.

    :param chi: ChI string
    :type chi: str
    :returns: False if the ChI is the inverted form of its canonical
        enantiomer; True in all other cases (including non-enantiomers).
    :rtype: bool
    """
    return is_canonical_enantiomer_list(split(chi))


def is_canonical_enantiomer_list(chis):
    """Is this list of ChIs a canonical combination of enantiomers?

    Sorts them as they would appear in a multi-component AMChI and checks
    that the first enantiomer (if any) is canonical.

    :param chis: A list of ChIs
    :type chis: list[str]
    :returns: Whether or not the list is canonical
    :rtype: bool
    """
    chis = sorted_(chis)
    invs = list(map(is_inverted_enantiomer, chis))
    inv = next((i for i in invs if i is not None), None)
    can = inv in (False, None)
    return can


def is_canonical_enantiomer_reaction(rct_chi, prd_chi):
    """Does this reaction have a canonical combination of enantiomers?

    :param rct_chi: A multi-component ChI or list of ChIs for the reactants
    :type rct_chi: str or list[str]
    :param prd_chi: A multi-component ChI or list of ChIs for the products
    :type prd_chi: str or list[str]

    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    rct_chis = split(rct_chi) if isinstance(rct_chi, str) else rct_chi
    prd_chis = split(prd_chi) if isinstance(prd_chi, str) else prd_chi

    # Switch to the canonical reaction direction
    if not is_canonical_reaction_direction(rct_chis, prd_chis):
        rct_chis, prd_chis = prd_chis, rct_chis

    chis = sorted_(rct_chis) + sorted_(prd_chis)
    can = is_canonical_enantiomer_list(chis)
    return can


def is_canonical_reaction_direction(rct_chis, prd_chis):
    """Is this the canonical reaction direction, or should it be reversed?

    :param rct_chis: A list of ChIs for the reactants
    :type rct_chis: list[str]
    :param prd_chis: A list of ChIs for the products
    :type prd_chis: list[str]
    :returns: Whether or not the reaction is canonical
    :rtype: bool
    """
    nrcts = len(rct_chis)
    nprds = len(prd_chis)

    idxs = argsort(rct_chis + prd_chis)
    rct_idxs = sorted(idxs[:nrcts])
    prd_idxs = sorted(idxs[nrcts:])

    rct_rep = (nrcts, rct_idxs)
    prd_rep = (nprds, prd_idxs)
    return rct_rep < prd_rep


def is_enantiomer_list(chis):
    """Does this list of species form an enantiomer?

    Often true if there are enantiomers in the list, but not always.  If,
    for every enantiomer in the list, its mirror image is also in the list,
    then this combination of species is achiral.

    :param chis: A list of ChIs
    :type chis: list[str]
    :returns: Whether or not the list is chiral
    :rtype: bool
    """
    orig_chis = sorted(chis)
    refl_chis = sorted(map(reflect, chis))
    return orig_chis != refl_chis


def is_enantiomer_reaction(rct_chis, prd_chis):
    """Is this an enantiomer reaction? I.e., is it chiral?

    :param rct_chis: A list of ChIs for the reactants
    :type rct_chis: list[str]
    :param prd_chis: A list of ChIs for the products
    :type prd_chis: list[str]
    :returns: Whether or not the reaction is chiral
    :rtype: bool
    """
    return is_enantiomer_list(rct_chis) or is_enantiomer_list(prd_chis)


# # # TS layers
def breaking_bond_keys(chi, one_indexed: bool = False):
    """Get the breaking bond keys of a TS AMChI

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :return: A set of bond keys
    :rtype: frozenset[frozenset[int]]
    """
    lyr_dct = ts_layers(chi)

    if "k" not in lyr_dct:
        return frozenset()

    bkeys = frozenset(_bonds(lyr_dct["k"], one_indexed=one_indexed).keys())
    return bkeys


def forming_bond_keys(chi, one_indexed: bool = False):
    """Get the forming bond keys of a TS AMChI

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :return: A set of bond keys
    :rtype: frozenset[frozenset[int]]
    """
    lyr_dct = ts_layers(chi)

    if "f" not in lyr_dct:
        return frozenset()

    bkeys = frozenset(_bonds(lyr_dct["f"], one_indexed=one_indexed).keys())
    return bkeys


def is_reversed_ts(chi):
    """Is this a reversed TS AMChI?

    :param chi: ChI string
    :type chi: str
    :return: `True` if it is a reversed TS, `False` if it is a non-reversed TS, `None`
        if it isn't a TS
    :rtype: Optional[bool]
    """
    lyr_dct = ts_layers(chi)

    if "r" not in lyr_dct:
        return None

    assert lyr_dct["r"] in "01"
    return lyr_dct["r"] == "1"


# # # isotope layers
def bond_isotope_stereo_parities(chi, one_indexed=False):
    """Parse the bond stereo parities from the isotope layers.

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :returns: A dictionary mapping bond keys onto parities
    :rtype: dict[frozenset[int]: bool]
    """
    lyr_dct = isotope_layers(chi)
    if "b" not in lyr_dct:
        return {}

    ste_dct = _bonds(lyr_dct["b"], one_indexed=one_indexed)
    return ste_dct


def atom_isotope_stereo_parities(chi, one_indexed=False):
    """Parse the atom stereo parities from the isotope layers.

    :param chi: ChI string
    :type chi: str
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    :returns: A dictionary mapping atom keys onto parities
    :rtype: dict[int: bool]
    """
    iso_lyr_dct = isotope_layers(chi)
    atm_ste_dct = _atom_stereo_parities(iso_lyr_dct, one_indexed=one_indexed)
    return atm_ste_dct


def is_inverted_isotope_enantiomer(chi):
    """Determine enantiomer inversion from the isotope layers.

    :param chi: ChI string
    :type chi: str
    :returns: whether or not the ChI is inverted; returns None if not an
        enantiomer
    :rtype: bool
    """
    iso_lyr_dct = isotope_layers(chi)
    is_inv = _is_inverted_enantiomer(iso_lyr_dct)
    return is_inv


# # other properties
def has_multiple_components(chi):
    """Determine if the ChI string has multiple components.

    :param chi: ChI string
    :type chi: str
    :rtype: bool
    """
    fstr = formula_layer(chi)
    return ";" in chi or "*" in chi or "." in fstr or str.isdigit(fstr[0])


def has_stereo(chi):
    """Determine if the ChI string has stereochemistry information.

    :param chi: ChI string
    :type chi: str
    :rtype: bool
    """
    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)
    return bool(ste_dct or any(pfx in iso_dct for pfx in STE_PFXS))


def has_mobile_hydrogens(chi):
    """Determine if the ChI string has mobile hydrogens.

    :param chi: ChI string
    :type chi: str
    :rtype: bool
    """
    main_lyr_dct = main_layers(chi)
    nhyd_lyr = main_lyr_dct["h"] if "h" in main_lyr_dct else ""
    return "(" in nhyd_lyr


def low_spin_multiplicity(chi):
    """Guess spin multiplicity based on the number of electrons.

    :param chi: ChI string
    :type chi: str
    :rtype: int
    """

    fml = formula(chi)
    nelec = form.electron_count(fml)

    if (nelec % 2) == 0:
        mult = 1
    else:
        mult = 2

    return mult


def is_enantiomer(chi, iso=True):
    """Is this ChI an enantiomer? (I.e., is it chiral?).

    Determined based on whether or not the ChI has an s-layer.

    :param ich: ChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :returns: whether or not the ChI is an enantiomer
    :rtype: bool
    """
    ret = "s" in stereo_layers(chi)
    if iso:
        ret |= "s" in isotope_layers(chi)
    return ret


def is_racemic(chi: str, iso: bool = True) -> bool:
    """Determine whether this is a racemic enantiomer.

    :param ich: ChI string
    :param iso: Include isotope stereochemistry?
    :return: `True` if it is, `False` if it isn't
    """
    ret = stereo_layers(chi).get("s", None) == "3"
    if iso:
        ret |= isotope_layers(chi).get("s", None) == "3"
    return ret


# # comparisons
def same_connectivity(chi1, chi2):
    """Determine if two ChI strings have the same connectivity.

    :param chi1: ChI string 1
    :type chi1: str
    :param chi2: ChI string 2
    :type chi2: str
    :rtype: bool
    """
    fml_str1 = formula_layer(chi1)
    fml_str2 = formula_layer(chi2)
    conn_dct1 = main_layers(chi1)
    conn_dct2 = main_layers(chi2)
    chg_dct1 = charge_layers(chi1)
    chg_dct2 = charge_layers(chi2)
    return fml_str1 == fml_str2 and conn_dct1 == conn_dct2 and chg_dct1 == chg_dct2


def equivalent(chi1, chi2):
    """Determine if two ChI strings are equivalent. Currently
    the srings are only checked up to the isotope sublayer.

    :param chi1: ChI string 1
    :type chi1: str
    :param chi2: ChI string 2
    :type chi2: str
    :rtype: bool
    """
    fml_str1 = formula_layer(chi1)
    fml_str2 = formula_layer(chi2)
    conn_dct1 = main_layers(chi1)
    conn_dct2 = main_layers(chi2)
    chg_dct1 = charge_layers(chi1)
    chg_dct2 = charge_layers(chi2)
    ste_dct1 = stereo_layers(chi1)
    ste_dct2 = stereo_layers(chi2)
    iso_dct1 = isotope_layers(chi1)
    iso_dct2 = isotope_layers(chi2)
    # Stereo layers get dropped upon split/joins, so remove these from the
    # equivalence test
    for dct in (ste_dct1, ste_dct2, iso_dct1, iso_dct2):
        if "s" in dct:
            dct.pop("s")
    return (
        fml_str1 == fml_str2
        and conn_dct1 == conn_dct2
        and chg_dct1 == chg_dct2
        and ste_dct1 == ste_dct2
        and iso_dct1 == iso_dct2
    )


# # split/join
def split(chi):
    """Split a multi-component ChI into ChIs for each of its components.

    :param chi: ChI string
    :type chi: str
    :returns: the split ChI strings
    :rtype: tuple[str]
    """
    pfx = prefix(chi)
    fml_str = formula_layer(chi)
    main_dct = main_layers(chi)
    char_dct = charge_layers(chi)
    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)
    ts_dct = ts_layers(chi)
    fml_strs = split_layer(fml_str)
    count = len(fml_strs)

    main_dcts = split_layers(main_dct, count)
    char_dcts = split_layers(char_dct, count)
    ste_dcts = split_layers(ste_dct, count)
    iso_dcts = split_layers(iso_dct, count)
    ts_dcts = split_layers(ts_dct, count)

    chis = tuple(
        from_data(
            fml_lyr=fml_str,
            main_lyr_dct=main_dct,
            char_lyr_dct=char_dct,
            ste_lyr_dct=ste_dct,
            iso_lyr_dct=iso_dct,
            ts_lyr_dct=ts_dct,
            pfx=pfx,
        )
        for fml_str, main_dct, char_dct, ste_dct, iso_dct, ts_dct in zip(
            fml_strs, main_dcts, char_dcts, ste_dcts, iso_dcts, ts_dcts
        )
    )
    return chis


def join(chis, sort=True):
    """Join separate ChI strings into one multi-component ChI string.

    :param chis: sequence of ChI strings
    :type chis: tuple[str]
    :param sort: sort the ChI strings in the standard sort order?
    :type sort: bool
    :returns: the joined ChI string
    :rtype: str
    """
    if not sort:
        return _join(chis)

    chis = list(itertools.chain(*map(split, chis)))
    chi, _ = sorted_join(chis)
    return chi


def sorted_join(chis):
    """Sort and join ChI strings, returning the sort order

    :param chis: sequence of ChI strings
    :type chis: tuple[str]
    """
    assert not any(
        len(split(c)) > 1 for c in chis
    ), f"List contains multi-component ChIs: {chis}"

    srt = argsort(chis)
    chi = _join([chis[i] for i in srt])
    return chi, srt


def _join(chis):
    """Join separate ChI strings into one multi-component ChI string.

    (Unsorted)

    :param chis: sequence of ChI strings
    :type chis: tuple[str]
    :param sort: sort the ChI strings in the standard sort order?
    :type sort: bool
    :returns: the joined ChI string
    :rtype: str
    """

    pfxs = list(map(prefix, chis))
    assert len(set(pfxs)) == 1, f"ChIs in join must share same prefix: {chis}"
    pfx, *_ = pfxs

    fml_strs = list(map(formula_layer, chis))
    fml_str = join_layer_strings(fml_strs, count_sep="", sep=".")
    main_dct, _ = join_layers(list(map(main_layers, chis)))
    char_dct, _ = join_layers(list(map(charge_layers, chis)))
    ste_dct, warn_msg1 = join_layers(list(map(stereo_layers, chis)))
    iso_dct, warn_msg2 = join_layers(list(map(isotope_layers, chis)))
    ts_dct, _ = join_layers(list(map(ts_layers, chis)))

    warn_msg = warn_msg1 if warn_msg1 else warn_msg2
    if warn_msg:
        warnings.warn(f"\n{chis}\n{warn_msg}")

    return from_data(
        fml_lyr=fml_str,
        main_lyr_dct=main_dct,
        char_lyr_dct=char_dct,
        ste_lyr_dct=ste_dct,
        iso_lyr_dct=iso_dct,
        ts_lyr_dct=ts_dct,
        pfx=pfx,
    )


# # sort
def sorted_(chis):
    """Sort a sequence of ChI strings in the standard sort order (see argsort)

    :param chis: sequence of ChI strings
    :type chis: tuple(str)
    :rtype: tuple(str)
    """
    return tuple(chis[i] for i in argsort(chis))


def argsort(chis):
    """Determine the standard sort order for multiple ChIs.

    Follows the sort order for multicomponent InChIs as much as possible.

    :param chis: sequence of ChI strings
    :type chis: tuple(str)
    """
    # 1. Formula sort values
    fmls = list(map(formula, chis))
    symbs = set(itertools.chain(*[f.keys() for f in fmls]))
    symbs = form.sorted_symbols(symbs, symbs_first=("C"))
    fml_vecs = [form.sort_vector(fml, symbs=symbs) for fml in fmls]

    # 2. Connectivity sort values
    conn_vecs = list(map(adjacency_list, chis))

    # 3. Hydrogen sort values
    idxs_lst = list(map(canonical_indices, chis))
    nhyd_dcts = list(map(hydrogen_valences, chis))
    nhyd_vecs = [list(map(h.__getitem__, i)) for i, h in zip(idxs_lst, nhyd_dcts)]

    # 4. Charge sort values
    char_vecs = list(map(charge, chis))

    # 5. Bond stereo sort values
    def _sorted_bond_keys(bnd_keys):
        srt_bnd_keys = sorted(sorted(k, reverse=True) for k in bnd_keys)
        return map(frozenset, srt_bnd_keys)

    bste_dcts = list(map(bond_stereo_parities, chis))
    bkeys_lst = list(map(_sorted_bond_keys, bste_dcts))
    bste_vecs = [list(map(s.__getitem__, k)) for k, s in zip(bkeys_lst, bste_dcts)]

    # 6. Atom stereo sort values
    aste_dcts = list(map(atom_stereo_parities, chis))
    akeys_lst = list(map(sorted, aste_dcts))
    aste_vecs = [list(map(s.__getitem__, k)) for k, s in zip(akeys_lst, aste_dcts)]

    # 7. Enantiomeric inversions
    inv_vecs = list(map(is_inverted_enantiomer, chis))

    # *. Stop here for now. Add isotope layers when we actually need them.
    if any(map(isotope_layers, chis)):
        raise NotImplementedError(
            "Isotope layers not yet implemented, but could easily be added."
        )

    # Do the actual sorting
    arr = numpy.empty(len(chis), dtype=object)
    arr[:] = list(
        zip(fml_vecs, conn_vecs, nhyd_vecs, char_vecs, bste_vecs, aste_vecs, inv_vecs)
    )
    idxs = tuple(reversed(numpy.argsort(arr)))

    return idxs


# # helpers
def _split_layer_string(lyr, count_delim=pp.Empty(), delim="."):
    """Split a layer string into components"""
    count = ppc.integer + pp.Suppress(count_delim)
    count = pp.Opt(count).setParseAction(lambda x: x[0] if x else [1])
    comp = pp.Opt(pp.Word(pp.printables, excludeChars=delim)).setParseAction(
        lambda x: x[0] if x else [""]
    )
    group = pp.Group(count + comp)
    groups = pp.delimitedList(group, delim=delim)
    grouped_comps = groups.parseString(lyr).asList()
    comps = tuple(itertools.chain(*([c] * n for n, c in grouped_comps)))
    return comps


def split_layer(lyr, key=""):
    """Split a layer into components

    :param lyr: The layer contents
    :type lyr: str
    :param key: They layer key (e.g. "t"), defaults to ""
    :type key: str, optional
    :return: The component layer contents
    :rtype: List[str]
    """
    if not key:
        lyrs = _split_layer_string(lyr)
    elif key == "m":
        lyrs = tuple(lyr)
    else:
        lyrs = _split_layer_string(lyr, count_delim="*", delim=";")
    return lyrs


def split_layers(lyr_dct, ncomps):
    """Split a multi-component layer dictionary into separate ones

    :param dct: A (potentially) multi-component layer dictionary
    :type dct: dict[str: str]
    :returns: the split layer dictionaries, along with a warning message to
        catch potential loss/corruption of information.
    :rtype: tuple[dict[str: str]], str
    """
    if not lyr_dct:
        lyr_dcts = ({},) * ncomps
        return lyr_dcts

    sval = lyr_dct.pop("s") if "s" in lyr_dct else ""

    keys = lyr_dct.keys()
    lyr_lsts = zip(*[split_layer(lyr, key) for key, lyr in lyr_dct.items()])
    lyr_dcts = [dict(zip(keys, lst)) for lst in lyr_lsts]
    lyr_dcts = [
        dict_.filter_by_value(d, lambda x: x not in ("", ".")) for d in lyr_dcts
    ]

    # Handle the /s layer (not possible to handle it perfectly for racemics)
    if sval:
        lyr_dcts = [{**d, "s": "1"} if "m" in d else d for d in lyr_dcts]

    return lyr_dcts


def join_layers(dcts):
    """Join all of the components of a ChI layer

    For the s layer, the highest value is chosen.

    Warning: For racemic mixtures with /s3 layers, there is a potential for
    corruption/loss of information, because it is no longer possible to
    tell if they are chiral or not from the /m layer.

    :param dcts: layer components, grouped by prefix
    :type dcts: tuple[dict[str: str]]
    :returns: the joined layer dictionary, along with a warning message
    :rtype: dict[str: str], str
    """

    pfxs = sorted(functools.reduce(set.union, map(set, dcts)))
    dcts = [dict_.by_key(dct, pfxs, fill_val="") for dct in dcts]
    lyrs_lst = [[dct[pfx] for dct in dcts] for pfx in pfxs]
    dct = {
        pfx: (
            join_layer_strings(lyrs)
            if pfx not in "ms"
            else (
                _join_m_layer_strings(lyrs)
                if pfx != "s"
                else _join_s_layer_strings(lyrs)
            )
        )
        for pfx, lyrs in zip(pfxs, lyrs_lst)
    }

    warn_msg = ""
    if "s" in dct and dct["s"] == "3":
        if any(("t" in d and "," in d["t"] and "m" not in d) for d in dcts):
            warn_msg = (
                "Potential loss/corruption of information!\n"
                "Potentially achiral species are being joined under "
                "an /s3 layer!"
            )

    return dct, warn_msg


def join_layer_strings(lyrs, count_sep="*", sep=";"):
    """Join layer strings into one multi-component layer string.

    :param lyrs: ChI layer contents to join
    :type lyrs: tuple[str]
    :param count_sep: delimiter for putting numbers before a repeated
        component, e.g. 2*H2O
    :type count_sep: str
    :param sep: delimiter between components, e.g. H2O;CH4
    :type sep: str
    """

    def _s(count, lyr):
        if count > 1 and lyr:
            ret = ("{:d}" + count_sep + "{:s}").format(count, lyr)
        elif lyr:
            ret = lyr
        else:
            ret = sep * (count - 1)
        return ret

    counts, lyrs = zip(*[(len(list(g)), lyr) for lyr, g in itertools.groupby(lyrs)])

    lyr = sep.join([_s(count, lyr) for count, lyr in zip(counts, lyrs)])
    return lyr


# # common multilayer properties
def _bonds(lyr, one_indexed=False, ordered_key=False):
    """Parse bond stereo parities from a given layer dictionary"""
    # Set up the parser
    integer = ppc.integer
    dash = pp.Suppress("-")
    bond = pp.Group(integer + dash + integer)
    parity = pp.Or(["+", "-", "?"])
    term = pp.Group(bond + pp.Opt(parity, default=None))
    parser = pp.Opt(pp.delimitedList(term, ","))
    key_type_ = tuple if ordered_key else frozenset

    # Do the parsing
    shift = 0 if one_indexed else -1
    lst = parser.parseString(lyr).asList()
    bnd_ste_dct = {
        key_type_((k1 + shift, k2 + shift)): (
            True if p == "+" else False if p == "-" else None
        )
        for (k1, k2), p in lst
    }
    return bnd_ste_dct


def _atom_stereo_parities(lyr_dct, one_indexed=False):
    """Parse atom stereo parities from a given layer dictionary"""
    if "t" not in lyr_dct:
        return {}

    lyr = lyr_dct["t"]

    # Set up the parser
    shift = 0 if one_indexed else -1
    integer = ppc.integer
    parity = pp.Or(["+", "-", "?"])
    term = pp.Group(integer + parity)
    parser = pp.Opt(pp.delimitedList(term, ","))

    # Do the parsing
    lst = parser.parseString(lyr).asList()
    atm_ste_dct = {k + shift: None if p == "?" else (p == "+") for k, p in lst}
    return atm_ste_dct


def _is_inverted_enantiomer(lyr_dct):
    """Determine enantiomer inversion from a given layer dictionary."""
    is_inv = None
    if "m" in lyr_dct:
        if lyr_dct["m"] == "1":
            is_inv = True
        else:
            assert lyr_dct["m"] == "0"
            is_inv = False
    return is_inv


# # split/join helpers
def _join_s_layer_strings(s_lyrs):
    s_lyrs = [s_lyr for s_lyr in s_lyrs if s_lyr]
    return max(s_lyrs, key=int)


def _join_m_layer_strings(m_lyrs):
    m_lyrs = [m_lyr if m_lyr else "." for m_lyr in m_lyrs]
    return "".join(m_lyrs)
