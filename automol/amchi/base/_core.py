""" Level 2 AMChI functions (depend on L1)

The parsing functions apply equally well to InChI or AMChI strings, so the
documentation simply refers to "ChI" strings.
"""

import itertools
import functools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
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
def is_chiral(chi):
    """ Determine if the ChI string has chirality information.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)
    return ste_dct['s'] == '1' or iso_dct['s'] == '1'


def is_enantiomer(chi, iso=True):
    """ Is this ChI an enantiomer?

        :param chi: ChI string
        :type chi: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the ChI is enantiomeric
        :rtype: bool
    """
    ste_dct = stereo_layers(chi)
    ret = 'm' in ste_dct
    if iso:
        iso_dct = isotope_layers(chi)
        ret = ret or 'm' in iso_dct
    return ret


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


def charge(chi):
    """ Determine charge from the ChI string

        :param chi: ChI string
        :type chi: str
        :rtype: int
    """
    char_lyr_dct = charge_layers(chi)
    char = int(char_lyr_dct['q'])
    return char


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


def canonical_indices(chi, one_indexed=False):
    """ Determine the list of canonical indices for a ChI string
    """
    idxs = sorted(symbols(chi, one_indexed=one_indexed).keys())
    return idxs


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


def hydrogen_valences(chi, one_indexed=False):
    """ Determine the hydrogen valences of backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
        :returns: a dictionary of hydrogen valences, keyed by canonical index
        :rtype: dict[int: int]
    """
    idxs = canonical_indices(chi)
    print(chi)
    print(idxs)


def bonds(chi, one_indexed=False):
    """ Determine bonds between backbone atoms in a ChI string

        :param chi: ChI string
        :type chi: str
        :param one_indexed: use one-indexing?
        :type one_indexed: bool
    """
    # TODO: use apf.all_captures_with_spans to find all indices in the
    # connection layer and identify what comes before/after
    main_lyr_dct = main_layers(chi)
    conn_lyr = main_lyr_dct['c'] if 'c' in main_lyr_dct else ''
    ptt = app.capturing(app.UNSIGNED_INTEGER)
    lst = apf.all_captures_with_spans(ptt, conn_lyr)
    print(lst)


def stereo_atoms(chi, iso=True, one_indexed=False):
    """ Parse the stereo atoms from the stereochemistry layer.

        :param chi: ChI string
        :type chi: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    assert len(split(chi)) == 1, (
        "Not for multicomponent ChIs. Call inchi.split() first.")

    atm_ptt = (app.capturing(app.UNSIGNED_INTEGER) +
               app.one_of_these(list(map(app.escape, '+-'))))

    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)

    tlyr = ''
    if 't' in ste_dct:
        tlyr += ste_dct['t']

    if iso and 't' in iso_dct:
        tlyr += ',' + iso_dct['t']

    atms = ()
    if tlyr:
        atms = ap_cast(apf.all_captures(atm_ptt, tlyr))

    if not one_indexed:
        atms = tuple(i-1 for i in atms)
        atms = atms if atms is not None else ()

    return atms


def stereo_bonds(chi, iso=True, one_indexed=False):
    """ Parse the stereo bonds from the stereochemistry layer.

        :param chi: ChI string
        :type chi: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    assert len(split(chi)) == 1, (
        "Not for multicomponent ChIs. Call inchi.split() first.")

    bnd_ptt = '-'.join([app.capturing(app.UNSIGNED_INTEGER)]*2)

    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)

    blyr = ''
    if 'b' in ste_dct:
        blyr += ste_dct['b']

    if iso and 'b' in iso_dct:
        blyr += ',' + iso_dct['b']

    bnds = ()
    if blyr:
        bnds = ap_cast(apf.all_captures(bnd_ptt, blyr))

    if not one_indexed:
        bnds = tuple((i-1, j-1) for i, j in bnds)
        bnds = bnds if bnds is not None else ()

    return bnds


def unassigned_stereo_bonds(chi, iso=True, one_indexed=False):
    """ Parse the stereo bonds wth missing assignments from the stereochemistry
    layer.

        :param chi: ChI string
        :type chi: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    assert len(split(chi)) == 1, (
        "Not for multicomponent ChIs. Call inchi.split() first.")

    bnd_ptt = ('-'.join([app.capturing(app.UNSIGNED_INTEGER)]*2) +
               app.escape('?'))

    ste_dct = stereo_layers(chi)
    iso_dct = isotope_layers(chi)

    blyr = ''
    if 'b' in ste_dct:
        blyr += ste_dct['b']

    if iso and 'b' in iso_dct:
        blyr += ',' + iso_dct['b']

    bnds = ()
    if blyr:
        bnds = ap_cast(apf.all_captures(bnd_ptt, blyr))
        bnds = bnds if bnds is not None else ()

    if not one_indexed:
        bnds = tuple((i-1, j-1) for i, j in bnds)

    return bnds


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
    ACH = ('AMChI=1/C10H14ClFO/c1-7(9(6-12)10(13)5-11)8-3-2-4-8'
           '/h2-4,7,9-10,13H,5-6H2,1H3')
    ICH = ('InChI=1S/C10H14ClFO/c1-7(8-3-2-4-8)9(6-12)10(13)5-11'
           '/h2-4,7,9-10,13H,5-6H2,1H3')
    SYMB_DCT = symbols(ACH)
    print(SYMB_DCT)
    # hydrogen_valences(ACH, one_indexed=True)
    print(ACH)
    bonds(ACH)
