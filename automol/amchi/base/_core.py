""" Level 2 AMChI functions (depend on L1)

The parsing functions apply equally well to InChI or AMChI strings, so the
documentation simply refers to "ChI" strings.

Future task: Rewrite all of this to use the pyparsing module, rather than
autoparse. It will be much cleaner.
"""

import itertools
import functools
from collections import abc
import pyparsing as pp
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
NONSLASH = '[^/]'
NONSLASHES = app.one_or_more(NONSLASH)
SLASH = app.escape('/')
SLASH_OR_START = app.one_of_these([SLASH, app.STRING_START])
SLASH_OR_END = app.one_of_these([SLASH, app.STRING_END])


# # constructor
def from_data(fml_str, main_lyr_dct=None,
              char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None):
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
        :rtype: str
    """

    main_dct = dict_.empty_if_none(main_lyr_dct)
    char_dct = dict_.empty_if_none(char_lyr_dct)
    ste_dct = dict_.empty_if_none(ste_lyr_dct)
    iso_dct = dict_.empty_if_none(iso_lyr_dct)

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

    chi = '/'.join(['AMChI=1', fml_str] + main_lyrs + char_lyrs +
                   ste_lyrs + iso_lyrs)

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


# # conversions
def formula(ich):
    """ Generate a formula dictionary from a ChI string.

        :param ich: ChI string
        :type ich: str
        :rtype: dict[str: int]
    """
    sym_ptt = app.UPPERCASE_LETTER + app.zero_or_more(app.LOWERCASE_LETTER)
    num_ptt = app.maybe(app.UNSIGNED_INTEGER)
    ptt = app.capturing(sym_ptt) + app.capturing(num_ptt)

    def _connected_formula(ich):
        fml_str = formula_string(ich)
        fml = {s: int(n) if n else 1
               for s, n in apf.all_captures(ptt, fml_str)}
        return fml

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = split(ich)
    fmls = list(map(_connected_formula, ichs))
    fml = functools.reduce(automol.formula.join, fmls)

    return fml


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
    if 'H' in pool:
        pool.remove('H')

    symbs = [s for symb in pool for s in itertools.repeat(symb, fml[symb])]

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
    chain = pp.delimitedList(integer, delim='-')
    chains = chain + pp.ZeroOrMore(',' + chain)
    side_chain = pp.nestedExpr('(', ')', content=chains)
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
    block = integer + pp.ZeroOrMore(sep + integer) + 'H' + pp.Opt(integer)
    parser = pp.Opt(pp.Group(block) + pp.ZeroOrMore(sep + pp.Group(block)))

    # Do the parsing
    main_lyr_dct = main_layers(chi)
    nhyd_lyr = main_lyr_dct['h'] if 'h' in main_lyr_dct else ''
    nhyd_lsts = ap_cast(parser.parseString(nhyd_lyr).asList())

    # Interpret the list
    shift = 0 if one_indexed else -1
    all_idxs = canonical_indices(chi, one_indexed=one_indexed)
    nhyd_dct = dict_.by_key({}, all_idxs, fill_val=0)
    for nhyd_lst in nhyd_lsts:
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

    return nhyd_dct


# # # charge layers
def charge(chi):
    """ Determine charge from the ChI string

        :param chi: ChI string
        :type chi: str
        :rtype: int
    """
    char_lyr_dct = charge_layers(chi)
    char = int(char_lyr_dct['q'])
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
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the ChI is inverted; returns None if not an
            enantiomer
        :rtype: bool
    """
    ste_lyr_dct = stereo_layers(chi)
    is_inv = _is_inverted_enantiomer(ste_lyr_dct)
    return is_inv


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


def has_multiple_components(chi):
    """ Determine if the ChI string has multiple components.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    return len(split(chi)) > 1


def low_spin_multiplicity(chi):
    """ Guess spin multiplicity based on the number of electrons.

        :param chi: ChI string
        :type chi: str
        :rtype: int
    """

    fml = formula(chi)
    nelec = automol.formula.electron_count(fml) - charge(chi)

    if (nelec % 2) == 0:
        mult = 1
    else:
        mult = 2

    return mult


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
    fml_strs = _split_layer_string(
        fml_str, count_sep_ptt='', sep_ptt=app.escape('.'))
    count = len(fml_strs)

    main_dcts = _split_layers(main_dct, count)
    char_dcts = _split_layers(char_dct, count)
    ste_dcts = _split_layers(ste_dct, count)
    iso_dcts = _split_layers(iso_dct, count)

    chis = tuple(from_data(fml_str=fml_str,
                           main_lyr_dct=main_dct,
                           char_lyr_dct=char_dct,
                           ste_lyr_dct=ste_dct,
                           iso_lyr_dct=iso_dct)
                 for fml_str, main_dct, char_dct, ste_dct, iso_dct
                 in zip(fml_strs, main_dcts, char_dcts, ste_dcts, iso_dcts))
    return chis


def join(chis):
    """ Join separate ChI strings into one multi-component ChI string.

        :param chis: sequence of ChI strings
        :type chis: tuple[str]
        :returns: the joined ChI string
        :rtype: str
    """
    # first, make sure they are completely split up
    chis = list(itertools.chain(*map(split, chis)))
    fml_strs = list(map(formula_string, chis))
    fml_str = _join_layer_strings(fml_strs, count_sep='', sep='.')
    main_dct = _join_layers(list(map(main_layers, chis)))
    char_dct = _join_layers(list(map(charge_layers, chis)))
    ste_dct = _join_layers(list(map(stereo_layers, chis)))
    iso_dct = _join_layers(list(map(isotope_layers, chis)))

    return from_data(fml_str=fml_str,
                     main_lyr_dct=main_dct,
                     char_lyr_dct=char_dct,
                     ste_lyr_dct=ste_dct,
                     iso_lyr_dct=iso_dct)


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
        parser = pp.Opt(pp.delimitedList(term, delim=','))

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
        parser = pp.Opt(pp.delimitedList(term, delim=','))

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
def _join_layers(dcts):
    """ Join all of the components of a ChI layer.

        :param dcts: layer components, grouped by prefix
        :type dct: dict[str: str]
        :rtype: dict[str: str]
    """

    pfxs = sorted(functools.reduce(set.union, map(set, dcts)))
    if 's' in pfxs:
        pfxs.remove('s')
    dcts = [dict_.by_key(dct, pfxs, fill_val='') for dct in dcts]
    lyrs_lst = [[dct[pfx] for dct in dcts] for pfx in pfxs]
    dct = {pfx: (_join_layer_strings(lyrs) if pfx != 'm' else
                 _join_m_layer_strings(lyrs))
           for pfx, lyrs in zip(pfxs, lyrs_lst)}

    return dct


def _join_m_layer_strings(m_lyrs):
    m_lyrs = [m_lyr if m_lyr else '.' for m_lyr in m_lyrs]
    return ''.join(m_lyrs)


def _join_layer_strings(lyrs, count_sep='*', sep=';'):
    """ Join layer strings into one multi-component layer string.

        :param lyrs: layers to join
        :type lyrs: tuple(str)?
        :param count_sep: delimiter for ???
        :type count_sep: str
        :param sep: delimiter for ???
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


def _split_layers(dct, count):
    """ split a multi-component layer dictionary into separate ones
    """
    if dct:
        pfxs = sorted(dct.keys())
        if 's' in pfxs:
            pfxs.remove('s')
        lyrs_lst = [
            _split_layer_string(dct[pfx]) if pfx != 'm'
            else _split_m_layer_string(dct[pfx]) for pfx in pfxs]
        assert all(len(lyrs) == count for lyrs in lyrs_lst)
        dcts = tuple({pfx: lyr for pfx, lyr in zip(pfxs, lyrs) if lyr}
                     for lyrs in zip(*lyrs_lst))
    else:
        return ({},) * count
    return dcts


def _split_m_layer_string(m_lyr):
    return tuple(m_lyr)


def _split_layer_string(lyr, count_sep_ptt=app.escape('*'),
                        sep_ptt=app.escape(';')):
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


if __name__ == '__main__':
    # ACH = ('AMChI=1/C10H14ClFO/c1-7(9(6-12)10(13)5-11)8-3-2-4-8'
    #        '/h2-4,7,9-10,13H,5-6H2,1H3')
    # ICH2 = 'InChI=1S/C3H8FNO2/c1-3(4,2-5)7-6/h6H,2,5H2,1H3'
    # BNDS = bonds(ACH, one_indexed=True)
    # print(BNDS)

    # Atom stereo
    CHI = 'AMChI=1/C3H3Cl2F3/c4-2(7)1(6)3(5)8/h1-3H/t2-,3-/m1/s1'
    print(atom_stereo_parities(CHI, one_indexed=True))
    print(is_inverted_enantiomer(CHI))
    CHI = 'InChI=1S/C3H8O/c1-3(2)4/h3-4H,1-2H3/i1+0,2+1/t3-/m1/s1'
    print(atom_stereo_parities(CHI, one_indexed=True))
    print(atom_isotope_stereo_parities(CHI, one_indexed=True))
    print(is_inverted_enantiomer(CHI))
    print(is_inverted_isotope_enantiomer(CHI))

    # # Bond stereo
    # CHI = 'AMChI=1/C3H5N3/c4-1-3(6)2-5/h1-2,4-6H/b4-1-,5-2+,6-3-'
    # print(bond_stereo_parities(CHI, one_indexed=True))
    # CHI = 'InChI=1S/C3H6/c1-3-2/h3H,1H2,2H3/i1D/b3-1+'
    # print(bond_isotope_stereo_parities(CHI, one_indexed=True))
