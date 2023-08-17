""" Level 2 AMChI functions (depend on L1)

The parsing functions apply equally well to InChI or AMChI strings, so the
documentation simply refers to "ChI" strings.

Future task: Rewrite all of this to use the pyparsing module, rather than
autoparse. It will be much cleaner.
"""

import itertools
import functools
import warnings
from collections import abc
import pyparsing as pp
import numpy
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
import automol.util
from automol.util import dict_
import automol.formula


MAIN_PFXS = ('c', 'h')
CHAR_PFXS = ('q', 'p')
STE_PFXS = ('b', 't', 'm', 's')
ISO_NONSTE_PFXS = ('i', 'h')
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS
TS_PFXS = ('k', 'f', 'r')
NONSLASH = '[^/]'
NONSLASHES = app.one_or_more(NONSLASH)
SLASH = app.escape('/')
SLASH_OR_START = app.one_of_these([SLASH, app.STRING_START])
SLASH_OR_END = app.one_of_these([SLASH, app.STRING_END])


# # constructor
def from_data(fml_str, main_lyr_dct=None,
              char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None, ts_lyr_dct=None):
    """ Build a ChI string from each of the various layers.

        :param fml_str: formula string, in Hill-sort order
        :type fml_str: str
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
        :rtype: str
    """

    main_dct = dict_.empty_if_none(main_lyr_dct)
    char_dct = dict_.empty_if_none(char_lyr_dct)
    ste_dct = dict_.empty_if_none(ste_lyr_dct)
    iso_dct = dict_.empty_if_none(iso_lyr_dct)
    ts_dct = dict_.empty_if_none(ts_lyr_dct)

    main_lyrs = [
        pfx + lyr for pfx, lyr
        in zip(MAIN_PFXS, dict_.values_by_key(main_dct, MAIN_PFXS)) if lyr]
    char_lyrs = [
        pfx + lyr for pfx, lyr
        in zip(CHAR_PFXS, dict_.values_by_key(char_dct, CHAR_PFXS)) if lyr]
    ste_lyrs = [
        pfx + lyr for pfx, lyr
        in zip(STE_PFXS, dict_.values_by_key(ste_dct, STE_PFXS)) if lyr]
    iso_lyrs = [
        pfx + lyr for pfx, lyr
        in zip(ISO_PFXS, dict_.values_by_key(iso_dct, ISO_PFXS)) if lyr]
    ts_lyrs = [
        pfx + lyr for pfx, lyr
        in zip(TS_PFXS, dict_.values_by_key(ts_dct, TS_PFXS)) if lyr]

    chi = '/'.join(['AMChI=1', fml_str] + main_lyrs + char_lyrs +
                   ste_lyrs + iso_lyrs + ts_lyrs)

    return chi


# # recalculate/standardize
def standard_form(chi, stereo=True, racem=False, ste_dct=None, iso_dct=None):
    """ Return a ChI string in standard form.

        Includes only the standard layers.

        :param chi: ChI string
        :type chi: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param racem: parameter to designate a racemic mixture, if chiral
        :type racem: bool
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
    fml_slyr = formula_string(chi)
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
        if 'm' in ste_dct:
            ste_dct.pop('m')
            ste_dct['s'] = '3'

        if 'm' in iso_dct:
            iso_dct.pop('m')
            iso_dct['s'] = '3'

    chi = from_data(fml_slyr,
                    main_lyr_dct=main_dct,
                    char_lyr_dct=char_dct,
                    ste_lyr_dct=ste_dct,
                    iso_lyr_dct=iso_dct,
                    ts_lyr_dct=ts_dct)

    # Make sure multi-component ChIs are properly sorted
    if has_multiple_components(chi):
        chi = join(split(chi))

    return chi


# # getters
def prefix(chi):
    """ Determine the chemical identifier prefix (InChI or AMChI).

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    ptt = app.capturing(NONSLASHES) + app.followed_by('=')
    std = apf.first_capture(ptt, chi)
    return std


def version(chi):
    """ Determine version number of the ChI string.

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    ptt = app.capturing(_version_pattern())
    ver = apf.first_capture(ptt, chi)
    return ver


def formula_string(chi):
    """ Parse the ChI string for the formula string.

        :param chi: ChI string
        :type chi: str
        :returns: the formula string
        :rtype: str
    """
    ptt = (_version_pattern() +
           SLASH + app.capturing(_formula_pattern()))
    lyr_str = apf.first_capture(ptt, chi)
    return lyr_str


def main_layers(chi):
    """ Parse the ChI string for the main layers ('c' and 'h').

        :param chi: ChI string
        :type chi: str
        :returns: the main layers, as a dictionary keyed by layer prefixes
        :rtype: dict[str: str]
    """
    ptt = (_version_pattern() +
           SLASH + _formula_pattern() +
           SLASH + app.capturing(_main_layers_pattern()))
    lyrs_str = apf.first_capture(ptt, chi)
    return _layers(lyrs_str)


def charge_layers(chi):
    """ Parse the ChI string for the charge layers ('q' and 'p').

        :param chi: ChI string
        :type chi: str
        :returns: the charge layers, as a dictionary keyed by layer prefixes
        :rtype: dict[str: str]
    """
    ptt = (_version_pattern() +
           SLASH + _formula_pattern() +
           app.maybe(SLASH + _main_layers_pattern()) +
           SLASH + app.capturing(_charge_layers_pattern()))
    lyrs_str = apf.first_capture(ptt, chi)
    return _layers(lyrs_str)


def stereo_layers(chi):
    """ Parse the ChI string for the stereo layers ('b', 't', 'm', and 's')

        :param chi: ChI string
        :type chi: str
        :returns: the stereo layers, as a dictionary keyed by layer prefixes
        :rtype: dict[str: str]
    """
    ptt = (_version_pattern() +
           SLASH + _formula_pattern() +
           app.maybe(SLASH + _main_layers_pattern()) +
           app.maybe(SLASH + _charge_layers_pattern()) +
           SLASH + app.capturing(_stereo_layers_pattern()))
    lyrs_str = apf.first_capture(ptt, chi)
    return _layers(lyrs_str)


def isotope_layers(chi):
    """ Parse the ChI string for the isotope layers ('i', 'h', 'b', 't', 'm',
        and 's')

        :param chi: ChI string
        :type chi: str
        :returns: the isotope layers, as a dictionary keyed by layer prefixes
        :rtype: dict[str: str]
    """
    ptt = (_version_pattern() +
           SLASH + _formula_pattern() +
           app.maybe(SLASH + _main_layers_pattern()) +
           app.maybe(SLASH + _charge_layers_pattern()) +
           app.maybe(SLASH + _stereo_layers_pattern()) +
           SLASH + app.capturing(_isotope_layers_pattern()))
    lyrs_str = apf.first_capture(ptt, chi)
    return _layers(lyrs_str)


def ts_layers(chi):
    """ Parse the ChI string for the TS layers ('k', 'f', and 'r')

        :param chi: ChI string
        :type chi: str
        :returns: the TS layers, as a dictionary keyed by layer prefixes
        :rtype: dict[str: str]
    """
    lyr_dct = {}
    for key in ('k', 'f', 'r'):
        lyr = _parse_layer(chi, key)
        if lyr is not None:
            lyr_dct[key] = lyr

    if lyr_dct:
        assert 'r' in lyr_dct and ('k' in lyr_dct or 'f' in lyr_dct), (
            f"Invalid TS AMChI: {chi}"
        )
    return lyr_dct


def _parse_layer(chi, key):
    layer_start = f"/{key}"
    lyr = None
    if layer_start in chi:
        before = pp.Suppress(pp.Literal("AMChI") + ... + pp.Literal(layer_start))
        after = pp.Suppress(pp.stringEnd ^ pp.Char('/'))

        layer_parser = before + ... + after
        lyr = layer_parser.parseString(chi).asList()
        lyr = lyr[0] if lyr else None

    return lyr


# # setters
def with_inchi_prefix(chi):
    """ Return a ChI with InChI prefix, whether AMChI or InChI.

        :param chi: ChI string
        :type chi: str
        :returns: InChI string
        :rtype: str
    """
    pfx = prefix(chi)
    if pfx == 'AMChI':
        chi = 'InChI' + chi[5:]
    else:
        assert pfx == 'InChI', (
            f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return chi


def reflect(chi):
    """ If this is an enantiomer, flip to the other enantiomer by changing the
        m-layer

        :param chi: InChI string
        :type chi: str
        :returns: the other enantiomer
        :rtype: bool
    """
    if is_enantiomer(chi):
        ste_upd_dct = stereo_layers(chi)
        iso_upd_dct = isotope_layers(chi)

        refl_trans = str.maketrans('01', '10')
        if 'm' in ste_upd_dct:
            ste_upd_dct['m'] = ste_upd_dct['m'].translate(refl_trans)

        if 'm' in iso_upd_dct:
            iso_upd_dct['m'] = iso_upd_dct['m'].translate(refl_trans)

        chi = standard_form(chi, ste_dct=ste_upd_dct, iso_dct=iso_upd_dct)

    return chi


def canonical_enantiomer(chi):
    """ Convert this ChI string to a canonical enantiomer, if it isn't already

        Works for multi-component InChIs or lists of InChIs

        :param chi: ChI string
        :type chi: str
        :returns: The reflected ChI string, if it was a non-canonical
            enantiomer; otherwise, it returns the original ChI string
    """
    if not is_canonical_enantiomer(chi):
        chi = (reflect(chi) if isinstance(chi, str)
               else tuple(map(reflect, chi)))
    return chi


def reflect_reaction(rct_chis, prd_chis):
    """ Apply a reflection operation to a reaction.

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
    """ Convert this reaction into a canonical combination of enantiomers

        :param rct_chi: A multi-component ChI or list of ChIs for the reactants
        :type rct_chi: str or list[str]
        :param prd_chi: A multi-component ChI or list of ChIs for the products
        :type prd_chi: str or list[str]

        :returns: The reflected reaction, if it was a non-canonical combination
            of enantiomers; otherwise, it returns the original reaction.
    """
    if not is_canonical_enantiomer_reaction(rct_chi, prd_chi):
        rct_chi = (reflect(rct_chi) if isinstance(rct_chi, str)
                   else tuple(map(reflect, rct_chi)))
        prd_chi = (reflect(prd_chi) if isinstance(prd_chi, str)
                   else tuple(map(reflect, prd_chi)))
    return rct_chi, prd_chi


# # conversions
def formula(chi):
    """ Generate a formula dictionary from a ChI string.

        :param chi: ChI string
        :type chi: str
        :rtype: dict[str: int]
    """
    sym_ptt = app.UPPERCASE_LETTER + app.zero_or_more(app.LOWERCASE_LETTER)
    num_ptt = app.maybe(app.UNSIGNED_INTEGER)
    ptt = app.capturing(sym_ptt) + app.capturing(num_ptt)

    def _connected_formula(chi):
        fml_str = formula_string(chi)
        fml = {s: int(n) if n else 1
               for s, n in apf.all_captures(ptt, fml_str)}
        return fml

    # split it up to handle hard-coded molecules in multi-component chis
    chis = split(chi)
    fmls = list(map(_connected_formula, chis))
    fml = functools.reduce(automol.formula.join, fmls)

    return fml


def without_stereo(chi):
    """ Remove all stereo layers
    """

    return standard_form(chi, stereo=False)


def racemic(chi):
    """ If chiral, convert the ChI into a racemic mixture

        This drops the /m layer and replaces /s1 with /s3, indicating a racemic
        mixture. The chirality of the species is still implied by the presence
        of the /s layer.

        :param chi: ChI string
        :type chi: str
    """
    return standard_form(chi, racem=True)


def connectivity(chi, parse_connection_layer=True, parse_h_layer=True):
    """ Return the 'c' and 'h' layers of the connectivity string

        The user may also specify what combination of the two layers
        that they wish to return
    """

    # Read the two sublayers that are requested to be parsed
    conn_slyrs = main_layers(chi)

    if parse_connection_layer:
        cslyr = conn_slyrs.get('c', '')
    else:
        cslyr = ''

    if parse_h_layer:
        hslyr = conn_slyrs.get('h', '')
    else:
        hslyr = ''

    # Write the parts of the connectivity string based on what was parsed
    if cslyr and hslyr:
        _str = f'c{cslyr}/h{hslyr}'
    elif cslyr:
        _str = f'c{cslyr}'
    elif hslyr:
        _str = f'h{hslyr}'
    else:
        _str = None

    return _str


def are_enantiomers(chi_a, chi_b):
    """ Are these InChI enantiomers of eachother?

        :param chi: InChI string
        :type chi: str
    """
    ste_dct_a = stereo_layers(chi_a)
    ste_dct_b = stereo_layers(chi_b)
    enant = False
    if main_layers(chi_a) == main_layers(chi_b):
        if (len(ste_dct_b.keys()) == len(ste_dct_a.keys())
                and 'm' in ste_dct_a.keys()):
            if ste_dct_a['m'] != ste_dct_b['m']:
                if 't' in ste_dct_a.keys():
                    if ste_dct_a['t'] == ste_dct_b['t']:
                        enant = True
                else:
                    enant = True
    return enant


def are_diastereomers(chi_a, chi_b):
    """ Are these InChI diastereomers of each other?

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
                if 'b' in ste_dct_a.keys():
                    if ste_dct_a['b'] != ste_dct_b['b']:
                        diast = True
                elif 't' in ste_dct_a.keys():
                    if ste_dct_a['t'] != ste_dct_b['t']:
                        diast = True

    return diast


# # properties
# # # formula layer
def symbols(chi, one_indexed=False):
    """ Determine the atomic symbols of backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
        :returns: a dictionary of atomic symbols, keyed by canonical index
        :rtype: dict[int: str]
    """
    fml = formula(chi)
    pool = list(automol.formula.sorted_symbols(fml.keys(), symbs_first=['C']))

    # If there are only hydrogens, then one of them must be a backbone atom
    if set(pool) == {'H'}:
        symbs = ['H']
    # Otherwise, remove all hydrogens from the list of backbone atom symbols
    else:
        symbs = [s for symb in pool for s in itertools.repeat(symb, fml[symb])
                 if symb != 'H']

    shift = 1 if one_indexed else 0
    symb_dct = {k+shift: s for k, s in enumerate(symbs)}
    return symb_dct


def canonical_indices(chi, one_indexed=False):
    """ Determine the list of canonical indices for a ChI string
    """
    idxs = sorted(symbols(chi, one_indexed=one_indexed).keys())
    return idxs


# # # main layers
def bonds(chi, one_indexed=False):
    """ Determine bonds between backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
    """
    # Set up the pyparsing parser
    integer = pp.Word(pp.nums)
    chain = pp.delimitedList(integer, '-')
    chains = chain + pp.ZeroOrMore(',' + chain)
    side_chain = pp.nestedExpr('(', ')', chains)
    parser = pp.Opt(chain + pp.ZeroOrMore(side_chain + chain))

    # Do the parsing. This produces a nested list of numbers and commas
    # mirroring the connection layer
    main_lyr_dct = main_layers(chi)
    conn_lyr = main_lyr_dct['c'] if 'c' in main_lyr_dct else ''
    conn_lst = list(ap_cast(parser.parseString(conn_lyr).asList()))

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
                lsts = automol.util.breakby(obj, ',')

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
                assert isinstance(nei, int), (
                    f"Something is wrong. {nei} should be an integer.")

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
    """ Generate an adjacency list for backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
    """
    idxs = canonical_indices(chi, one_indexed=False)
    bnds = bonds(chi, one_indexed=False)
    adjs = [sorted(next(iter(b-{i})) for b in bnds if i in b) for i in idxs]
    return adjs


def hydrogen_valences(chi, one_indexed=False):
    """ Determine the hydrogen valences of backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
        :returns: a dictionary of hydrogen valences, keyed by canonical index
        :rtype: dict[int: int]
    """
    # Set up the parser
    integer = pp.Word(pp.nums)
    sep = '-' | pp.Suppress(',')
    fixedh = integer + pp.ZeroOrMore(sep + integer) + 'H' + pp.Opt(integer)
    mobileh = pp.Combine(pp.Suppress('(') + 'H' + pp.Opt(integer) +
                         pp.OneOrMore(',' + integer) + pp.Suppress(')'))
    parser = pp.Opt(
        pp.Opt(pp.Group(fixedh) + pp.ZeroOrMore(sep + pp.Group(fixedh))) +
        pp.Opt(sep) + pp.ZeroOrMore(mobileh))

    mobileh_parser = (pp.Combine('H' + pp.Opt(integer)) + sep +
                      pp.Group(pp.delimitedList(integer, ',')))

    # Do the parsing
    main_lyr_dct = main_layers(chi)
    nhyd_lyr = main_lyr_dct['h'] if 'h' in main_lyr_dct else ''
    nhyd_lsts = ap_cast(parser.parseString(nhyd_lyr).asList())

    # Interpret the list
    shift = 0 if one_indexed else -1
    all_idxs = canonical_indices(chi, one_indexed=one_indexed)
    nhyd_dct = dict_.by_key({}, all_idxs, fill_val=0)
    mob_lsts = []
    for nhyd_lst in nhyd_lsts:
        # Mobile hydrogen blocks will be strings instead of lists
        if isinstance(nhyd_lst, str) and nhyd_lst.startswith('H'):
            mob_lsts.append(nhyd_lst)
        # Otherwise, this is not a mobile hydrogen block -- proceed as usual
        else:
            if isinstance(nhyd_lst[-1], int):
                nhyd = nhyd_lst[-1]
                nhyd_lst = nhyd_lst[:-2]
            else:
                nhyd = 1
                nhyd_lst = nhyd_lst[:-1]

            lsts = list(map(list, automol.util.breakby(nhyd_lst, '-')))
            idxs = lsts.pop(0)
            for lst in lsts:
                idxs.extend(range(idxs[-1]+1, lst[0]))
                idxs.extend(lst)
            idxs = [k+shift for k in idxs]
            nhyd_dct.update({k: nhyd for k in idxs})

    # Add in mobile hydrogens after we get the others
    for mob_lst in mob_lsts:
        mob_lst = ap_cast(mobileh_parser.parseString(nhyd_lst).asList())
        nmob_str = mob_lst[0][1:]
        nmob = int(nmob_str) if nmob_str else 1
        # Add available mobile hydrogens to the first nmob atoms
        idxs = [k+shift for k in mob_lst[1]][:nmob]
        nhyd_dct.update({k: nhyd_dct[k] + 1 for k in idxs})

    return nhyd_dct


# # # charge layers
def charge(chi):
    """ Determine charge from the ChI string

        :param chi: ChI string
        :type chi: str
        :rtype: int
    """
    char_lyr_dct = charge_layers(chi)
    char = int(char_lyr_dct['q']) if 'q' in char_lyr_dct else 0
    return char


# # # stereo layers
def bond_stereo_parities(chi, one_indexed=False):
    """ Parse the bond stereo parities from the stereochemistry layers.

        :param chi: ChI string
        :type chi: str
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
        :returns: A dictionary mapping bond keys onto parities
        :rtype: dict[frozenset[int]: bool]
    """
    ste_lyr_dct = stereo_layers(chi)
    bnd_ste_dct = _bond_stereo_parities(ste_lyr_dct, one_indexed=one_indexed)
    return bnd_ste_dct


def atom_stereo_parities(chi, one_indexed=False):
    """ Parse the atom stereo parities from the stereochemistry layers.

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
    """ Determine enantiomer inversion from the stereo layers.

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
    """ Is this a canonical enantiomer? Also returns true for non-enantiomers.

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
    """ Is this list of ChIs a canonical combination of enantiomers?

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
    """ Does this reaction have a canonical combination of enantiomers?

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
    """ Is this the canonical reaction direction, or should it be reversed?

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
    """ Does this list of species form an enantiomer?

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
    """ Is this an enantiomer reaction? I.e., is it chiral?

        :param rct_chis: A list of ChIs for the reactants
        :type rct_chis: list[str]
        :param prd_chis: A list of ChIs for the products
        :type prd_chis: list[str]
        :returns: Whether or not the reaction is chiral
        :rtype: bool
    """
    return is_enantiomer_list(rct_chis) or is_enantiomer_list(prd_chis)


# # # isotope layers
def bond_isotope_stereo_parities(chi, one_indexed=False):
    """ Parse the bond stereo parities from the isotope layers.

        :param chi: ChI string
        :type chi: str
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
        :returns: A dictionary mapping bond keys onto parities
        :rtype: dict[frozenset[int]: bool]
    """
    iso_lyr_dct = isotope_layers(chi)
    bnd_ste_dct = _bond_stereo_parities(iso_lyr_dct, one_indexed=one_indexed)
    return bnd_ste_dct


def atom_isotope_stereo_parities(chi, one_indexed=False):
    """ Parse the atom stereo parities from the isotope layers.

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
    """ Determine enantiomer inversion from the isotope layers.

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
    """ Determine if the ChI string has multiple components.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    fstr = formula_string(chi)
    return ';' in chi or '*' in chi or '.' in fstr or str.isdigit(fstr[0])


def has_stereo(chi):
    """ Determine if the ChI string has stereochemistry information.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)
    return bool(ste_dct or
                any(pfx in iso_dct for pfx in STE_PFXS))


def has_mobile_hydrogens(chi):
    """ Determine if the ChI string has mobile hydrogens.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    main_lyr_dct = main_layers(chi)
    nhyd_lyr = main_lyr_dct['h'] if 'h' in main_lyr_dct else ''
    return '(' in nhyd_lyr


def low_spin_multiplicity(chi):
    """ Guess spin multiplicity based on the number of electrons.

        :param chi: ChI string
        :type chi: str
        :rtype: int
    """

    fml = formula(chi)
    nelec = automol.formula.electron_count(fml)

    if (nelec % 2) == 0:
        mult = 1
    else:
        mult = 2

    return mult


def is_enantiomer(chi, iso=True):
    """ Is this ChI an enantiomer? (I.e., is it chiral?)

        Determined based on whether or not the ChI has an s-layer.

        :param ich: ChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the ChI is an enantiomer
        :rtype: bool
    """
    ret = 's' in stereo_layers(chi)
    if iso:
        ret |= 's' in isotope_layers(chi)
    return ret


# # comparisons
def same_connectivity(chi1, chi2):
    """ Determine if two ChI strings have the same connectivity.

        :param chi1: ChI string 1
        :type chi1: str
        :param chi2: ChI string 2
        :type chi2: str
        :rtype: bool
    """
    return (standard_form(chi1, stereo=False) ==
            standard_form(chi2, stereo=False))


def equivalent(chi1, chi2):
    """ Determine if two ChI strings are equivalent. Currently
        the srings are only checked up to the isotope sublayer.

        :param chi1: ChI string 1
        :type chi1: str
        :param chi2: ChI string 2
        :type chi2: str
        :rtype: bool
    """
    fml_str1 = formula_string(chi1)
    fml_str2 = formula_string(chi2)
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
        if 's' in dct:
            dct.pop('s')
    return (fml_str1 == fml_str2 and conn_dct1 == conn_dct2 and
            chg_dct1 == chg_dct2 and ste_dct1 == ste_dct2 and
            iso_dct1 == iso_dct2)


# # split/join
def split(chi):
    """ Split a multi-component ChI into ChIs for each of its components.

        :param chi: ChI string
        :type chi: str
        :returns: the split ChI strings
        :rtype: tuple[str]
    """
    fml_str = formula_string(chi)
    main_dct = main_layers(chi)
    char_dct = charge_layers(chi)
    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)
    ts_dct = ts_layers(chi)
    fml_strs = split_layer_string(
        fml_str, count_sep_ptt='', sep_ptt=app.escape('.'))
    count = len(fml_strs)

    main_dcts, _ = split_layers(main_dct, count)
    char_dcts, _ = split_layers(char_dct, count)
    ste_dcts, warn_msg1 = split_layers(ste_dct, count)
    iso_dcts, warn_msg2 = split_layers(iso_dct, count)
    ts_dcts, _ = split_layers(ts_dct, count)

    warn_msg = warn_msg1 if warn_msg1 else warn_msg2
    if warn_msg:
        warnings.warn(f"\n{chi}\n{warn_msg}")

    chis = tuple(from_data(fml_str=fml_str,
                           main_lyr_dct=main_dct,
                           char_lyr_dct=char_dct,
                           ste_lyr_dct=ste_dct,
                           iso_lyr_dct=iso_dct,
                           ts_lyr_dct=ts_dct)
                 for fml_str, main_dct, char_dct, ste_dct, iso_dct, ts_dct
                 in zip(fml_strs, main_dcts, char_dcts, ste_dcts, iso_dcts, ts_dcts))
    return chis


def join(chis, sort=True):
    """ Join separate ChI strings into one multi-component ChI string.

        :param chis: sequence of ChI strings
        :type chis: tuple[str]
        :param sort: sort the ChI strings in the standard sort order?
        :type sort: bool
        :returns: the joined ChI string
        :rtype: str
    """
    # first, make sure they are completely split up
    chis = list(itertools.chain(*map(split, chis)))
    if sort:
        chis = sorted_(chis)

    fml_strs = list(map(formula_string, chis))
    fml_str = join_layer_strings(fml_strs, count_sep='', sep='.')
    main_dct, _ = join_layers(list(map(main_layers, chis)))
    char_dct, _ = join_layers(list(map(charge_layers, chis)))
    ste_dct, warn_msg1 = join_layers(list(map(stereo_layers, chis)))
    iso_dct, warn_msg2 = join_layers(list(map(isotope_layers, chis)))
    ts_dct, _ = join_layers(list(map(ts_layers, chis)))

    warn_msg = warn_msg1 if warn_msg1 else warn_msg2
    if warn_msg:
        warnings.warn(f"\n{chis}\n{warn_msg}")

    return from_data(fml_str=fml_str,
                     main_lyr_dct=main_dct,
                     char_lyr_dct=char_dct,
                     ste_lyr_dct=ste_dct,
                     iso_lyr_dct=iso_dct,
                     ts_lyr_dct=ts_dct)


# # sort
def sorted_(chis):
    """ Sort a sequence of ChI strings in the standard sort order (see argsort)

        :param chis: sequence of ChI strings
        :type chis: tuple(str)
        :rtype: tuple(str)
    """
    return tuple(chis[i] for i in argsort(chis))


def argsort(chis):
    """ Determine the standard sort order for multiple ChIs.

        Follows the sort order for multicomponent InChIs as much as possible.

        :param chis: sequence of ChI strings
        :type chis: tuple(str)
    """
    # 1. Formula sort values
    fmls = list(map(formula, chis))
    symbs = set(itertools.chain(*[f.keys() for f in fmls]))
    symbs = automol.formula.sorted_symbols(symbs, symbs_first=('C'))
    fml_vecs = [automol.formula.sort_vector(fml, symbs=symbs) for fml in fmls]

    # 2. Connectivity sort values
    conn_vecs = list(map(adjacency_list, chis))

    # 3. Hydrogen sort values
    idxs_lst = list(map(canonical_indices, chis))
    nhyd_dcts = list(map(hydrogen_valences, chis))
    nhyd_vecs = [
        list(map(h.__getitem__, i)) for i, h in zip(idxs_lst, nhyd_dcts)]

    # 4. Charge sort values
    char_vecs = list(map(charge, chis))

    # 5. Bond stereo sort values
    def _sorted_bond_keys(bnd_keys):
        srt_bnd_keys = sorted(sorted(k, reverse=True) for k in bnd_keys)
        return map(frozenset, srt_bnd_keys)

    bste_dcts = list(map(bond_stereo_parities, chis))
    bkeys_lst = list(map(_sorted_bond_keys, bste_dcts))
    bste_vecs = [
        list(map(s.__getitem__, k)) for k, s in zip(bkeys_lst, bste_dcts)]

    # 6. Atom stereo sort values
    aste_dcts = list(map(atom_stereo_parities, chis))
    akeys_lst = list(map(sorted, aste_dcts))
    aste_vecs = [
        list(map(s.__getitem__, k)) for k, s in zip(akeys_lst, aste_dcts)]

    # 7. Enantiomeric inversions
    inv_vecs = list(map(is_inverted_enantiomer, chis))

    # *. Stop here for now. Add isotope layers when we actually need them.
    if any(map(isotope_layers, chis)):
        raise NotImplementedError("Isotope layers not yet implemented,"
                                  "but could easily be added.")

    # Do the actual sorting
    arr = numpy.empty(len(chis), dtype=object)
    arr[:] = list(zip(fml_vecs, conn_vecs, nhyd_vecs, char_vecs,
                      bste_vecs, aste_vecs, inv_vecs))
    idxs = tuple(reversed(numpy.argsort(arr)))

    return idxs


# # helpers
def version_pattern():
    """ Build the autoparse regex pattern for the ChI string version.

        :rtype: str
    """
    return _version_pattern()


def join_layers(dcts):
    """ Join all of the components of a ChI layer

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
    dcts = [dict_.by_key(dct, pfxs, fill_val='') for dct in dcts]
    lyrs_lst = [[dct[pfx] for dct in dcts] for pfx in pfxs]
    dct = {pfx: (join_layer_strings(lyrs) if pfx not in 'ms' else
                 _join_m_layer_strings(lyrs) if pfx != 's' else
                 _join_s_layer_strings(lyrs))
           for pfx, lyrs in zip(pfxs, lyrs_lst)}

    warn_msg = ''
    if 's' in dct and dct['s'] == '3':
        if any(('t' in d and ',' in d['t'] and 'm' not in d) for d in dcts):
            warn_msg = ("Potential loss/corruption of information!\n"
                        "Potentially achiral species are being joined under "
                        "an /s3 layer!")

    return dct, warn_msg


def split_layers(dct, count):
    """ Split a multi-component layer dictionary into separate ones

        See warning in join_layers doc string about s layers. If splitting from
        a string with an /s3 layer, there is no way to determine which species
        with /t layers are chiral or not (without doing a full recalculation).
        All species with /t stereochemistry will inherit the /s3 layer.

        :param dct: A (potentially) multi-component layer dictionary
        :type dct: dict[str: str]
        :returns: the split layer dictionaries, along with a warning message to
            catch potential loss/corruption of information.
        :rtype: tuple[dict[str: str]], str
    """
    warn_msg = ''

    if not dct:
        dcts = ({},) * count
    else:
        sval = dct.pop('s') if 's' in dct else ''
        pfxs = sorted(dct.keys())
        lyrs_lst = [
            split_layer_string(dct[pfx]) if pfx != 'm'
            else _split_m_layer_string(dct[pfx]) for pfx in pfxs]
        assert all(len(lyrs) == count for lyrs in lyrs_lst)
        dcts = list({pfx: lyr for pfx, lyr in zip(pfxs, lyrs) if lyr}
                    for lyrs in zip(*lyrs_lst))
        for dct_ in dcts:
            if 'm' in dct_ and dct_['m'] == '.':
                dct_.pop('m')
            if 'm' in dct:
                dct_['s'] = sval
            if sval == '3' and 't' in dct_:
                dct_['s'] = sval
                if ',' in dct_['t']:
                    warn_msg = ("Potential loss/corruption of information!\n"
                                "Potentially achiral species are being split "
                                "from under an /s3 layer.\n"
                                "Any species with a /t layer will be marked "
                                "as chiral!")
    return dcts, warn_msg


def join_layer_strings(lyrs, count_sep='*', sep=';'):
    """ Join layer strings into one multi-component layer string.

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
            ret = ('{:d}' + count_sep + '{:s}').format(count, lyr)
        elif lyr:
            ret = lyr
        else:
            ret = sep * (count - 1)
        return ret

    counts, lyrs = zip(*[
        (len(list(g)), lyr) for lyr, g in itertools.groupby(lyrs)])

    lyr = sep.join([_s(count, lyr) for count, lyr in zip(counts, lyrs)])
    return lyr


def split_layer_string(lyr, count_sep_ptt=app.escape('*'),
                       sep_ptt=app.escape(';')):
    """ Split a multi-component layer string.

        :param lyr: ChI layer contents to split
        :type lyr: str
        :param count_sep: delimiter for putting numbers before a repeated
            component, e.g. 2*H2O
        :type count_sep: str
        :param sep: delimiter between components, e.g. H2O;CH4
        :type sep: str
    """
    count_ptt = app.UNSIGNED_INTEGER
    group_ptt = (app.STRING_START + app.capturing(count_ptt) + count_sep_ptt +
                 app.capturing(app.zero_or_more(app.WILDCARD)))

    def _expand_group(group_str):
        if apf.has_match(group_ptt, group_str):
            count, part = ap_cast(apf.first_capture(group_ptt, group_str))
            parts = [part] * count
        else:
            parts = [group_str]
        return parts

    parts = tuple(
        itertools.chain(*map(_expand_group, apf.split(sep_ptt, lyr))))
    return parts


# # common multilayer properties
def _bond_stereo_parities(lyr_dct, one_indexed=False):
    """ Parse bond stereo parities from a given layer dictionary
    """
    if 'b' not in lyr_dct:
        bnd_ste_dct = {}
    else:
        lyr = lyr_dct['b']

        # Set up the parser
        integer = pp.Word(pp.nums)
        bond = integer + pp.Suppress('-') + integer
        parity = pp.Or(['+', '-'])
        term = pp.Group(pp.Group(bond) + parity)
        parser = pp.Opt(pp.delimitedList(term, ','))

        # Do the parsing
        lst = ap_cast(parser.parseString(lyr).asList())

        # Interpret the list
        shift = 0 if one_indexed else -1
        bnd_ste_dct = {frozenset({k1+shift, k2+shift}): (p == '+')
                       for (k1, k2), p in lst}
    return bnd_ste_dct


def _atom_stereo_parities(lyr_dct, one_indexed=False):
    """ Parse atom stereo parities from a given layer dictionary
    """
    if 't' not in lyr_dct:
        atm_ste_dct = {}
    else:
        lyr = lyr_dct['t']

        # Set up the parser
        integer = pp.Word(pp.nums)
        parity = pp.Or(['+', '-'])
        term = pp.Group(integer + parity)
        parser = pp.Opt(pp.delimitedList(term, ','))

        # Do the parsing
        lst = ap_cast(parser.parseString(lyr).asList())

        # Interpret the list
        shift = 0 if one_indexed else -1
        atm_ste_dct = {k+shift: (p == '+') for k, p in lst}
    return atm_ste_dct


def _is_inverted_enantiomer(lyr_dct):
    """ Determine enantiomer inversion from a given layer dictionary.
    """
    is_inv = None
    if 'm' in lyr_dct:
        if lyr_dct['m'] == '1':
            is_inv = True
        else:
            assert lyr_dct['m'] == '0'
            is_inv = False
    return is_inv


# # parsing helpers
def _layers(lyrs_str):
    """ Parse the layers of the specified layer of a ChI string to a
        dictionary, keyed by prefix.

        :param lyrs_str: a string containing one or more ChI layers
        :type lyrs_str: str
        :returns: a dictionary of layers, keyed by layer prefix
        :rtype: dict[str: str]
    """
    if lyrs_str:
        ptt = _layer_pattern(key_ptt=app.capturing(app.LOWERCASE_LETTER),
                             val_ptt=app.capturing(NONSLASHES))
        dct = dict(apf.all_captures(ptt, lyrs_str))
    else:
        dct = {}
    return dct


# # parsing patterns
def _version_pattern():
    """ Build the autoparse regex pattern for the ChI string version.

        :rtype: str
    """
    ptt = app.preceded_by('=') + _layer_pattern()
    return ptt


def _formula_pattern():
    """ Build the autoparse regex pattern for the chemical formual layer.

        :rtype: str
    """
    ptt = _layer_pattern(key_ptt=app.not_followed_by(app.LOWERCASE_LETTER))
    return ptt


def _main_layers_pattern():
    """ Build the autoparse regex pattern for the connectivity layer.

        :rtype: str
    """
    c_lyr_ptt = _layer_pattern(key_ptt='c')
    h_lyr_ptt = _layer_pattern(key_ptt='h')
    ptt = (app.one_of_these([c_lyr_ptt, h_lyr_ptt]) +
           app.maybe(SLASH + h_lyr_ptt))
    return ptt


def _charge_layers_pattern():
    """ Build the autoparse regex pattern for the charge layer.

        :rtype: str
    """
    q_lyr_ptt = _layer_pattern(key_ptt='q')
    p_lyr_ptt = _layer_pattern(key_ptt='p')
    ptt = (app.one_of_these([q_lyr_ptt, p_lyr_ptt]) +
           app.maybe(SLASH + p_lyr_ptt))
    return ptt


def _stereo_layers_pattern():
    """ Build the autoparse regex pattern for the stereochemistry layer.

        :rtype: str
    """
    b_lyr_ptt = _layer_pattern(key_ptt='b')
    t_lyr_ptt = _layer_pattern(key_ptt='t')
    m_lyr_ptt = _layer_pattern(key_ptt='m')
    s_lyr_ptt = _layer_pattern(key_ptt='s')
    ptt = (app.one_of_these([b_lyr_ptt, t_lyr_ptt]) +
           app.maybe(SLASH + t_lyr_ptt) +
           app.maybe(SLASH + m_lyr_ptt) +
           app.maybe(SLASH + s_lyr_ptt))
    return ptt


def _isotope_layers_pattern():
    """ Build the autoparse regex pattern for the isotope layer.

        :rtype: str
    """
    i_lyr_ptt = _layer_pattern(key_ptt='i')
    h_lyr_ptt = _layer_pattern(key_ptt='h')
    ptt = (i_lyr_ptt +
           app.maybe(SLASH + h_lyr_ptt) +
           app.maybe(SLASH + _stereo_layers_pattern()))
    return ptt


def _layer_pattern(key_ptt='', val_ptt=NONSLASHES):
    """ Build the autoparse regex pattern for an arbitrary ChI layer.

        :rtype: str
    """
    return key_ptt + val_ptt


# # split/join helpers
def _join_s_layer_strings(s_lyrs):
    s_lyrs = [s_lyr for s_lyr in s_lyrs if s_lyr]
    return max(s_lyrs, key=int)


def _join_m_layer_strings(m_lyrs):
    m_lyrs = [m_lyr if m_lyr else '.' for m_lyr in m_lyrs]
    return ''.join(m_lyrs)


def _split_m_layer_string(m_lyr):
    return tuple(m_lyr)
