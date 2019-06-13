""" InChI chemical identifiers
"""
import itertools
import functools
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from automol import dict_
import automol.convert.inchi

NONSLASH = '[^/]'
NONSLASHES = app.one_or_more(NONSLASH)
SLASH = app.escape('/')
SLASH_OR_START = app.one_of_these([SLASH, app.STRING_START])
SLASH_OR_END = app.one_of_these([SLASH, app.STRING_END])


MAIN_PFXS = ['c', 'h']
CHAR_PFXS = ['q', 'p']
STE_PFXS = ['b', 't', 'm', 's']
ISO_NONSTE_PFXS = ['i', 'h']
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS


# "constructor"
def from_data(fml_slyr, main_dct=None, char_dct=None, ste_dct=None,
              iso_dct=None):
    """ calculate an inchi string from layers
    """
    main_dct = dict_.empty_if_none(main_dct)
    char_dct = dict_.empty_if_none(char_dct)
    ste_dct = dict_.empty_if_none(ste_dct)
    iso_dct = dict_.empty_if_none(iso_dct)

    main_slyrs = [
        pfx + slyr for pfx, slyr
        in zip(MAIN_PFXS, dict_.values_by_key(main_dct, MAIN_PFXS)) if slyr]
    char_slyrs = [
        pfx + slyr for pfx, slyr
        in zip(CHAR_PFXS, dict_.values_by_key(char_dct, CHAR_PFXS)) if slyr]
    ste_slyrs = [
        pfx + slyr for pfx, slyr
        in zip(STE_PFXS, dict_.values_by_key(ste_dct, STE_PFXS)) if slyr]
    iso_slyrs = [
        pfx + slyr for pfx, slyr
        in zip(ISO_PFXS, dict_.values_by_key(iso_dct, ISO_PFXS)) if slyr]

    ich = '/'.join(['InChI=1', fml_slyr] + main_slyrs + char_slyrs +
                   ste_slyrs + iso_slyrs)
    ich = automol.convert.inchi.recalculate(ich)
    return ich


# getters
def version(ich):
    """ version
    """
    ptt = app.capturing(_version_pattern())
    ver = apf.first_capture(ptt, ich)
    return ver


def formula_sublayer(ich):
    """ formula sublayer
    """
    ptt = (_version_pattern() +
           SLASH + app.capturing(_formula_sublayer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def main_sublayers(ich):
    """ main sublayers, by prefix
    """
    return _sublayers(_main_layer(ich))


def charge_sublayers(ich):
    """ charge sublayers, by prefix
    """
    return _sublayers(_charge_layer(ich))


def stereo_sublayers(ich):
    """ stereo sublayers, by prefix
    """
    return _sublayers(_stereo_layer(ich))


def isotope_sublayers(ich):
    """ isotope sublayers, by prefix
    """
    return _sublayers(_isotope_layer(ich))


# setters
def has_stereo(ich):
    """ does this inchi have stereo information?
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return bool(ste_dct or any(pfx in iso_dct for pfx in STE_PFXS))


def without_stereo(ich):
    """ inchi without stereo layer (and sublayers)
    """
    fml_slyr = formula_sublayer(ich)
    main_dct = main_sublayers(ich)
    char_dct = charge_sublayers(ich)
    iso_dct = dict_.by_key(isotope_sublayers(ich), ISO_NONSTE_PFXS)
    return from_data(fml_slyr=fml_slyr, main_dct=main_dct, char_dct=char_dct,
                     iso_dct=iso_dct)


def is_closed(ich):
    """ is this inchi closed?
    """
    return ich == recalculate(ich)


def is_complete(ich):
    """ is this inchi complete? (are all stereo-centers assigned?)
    """
    return is_closed(ich) and not (
        has_stereo(ich) ^ has_stereo(recalculate(ich, force_stereo=True)))


# comparisons
def equivalent(ich1, ich2):
    """ are these inchis equivalent? (only considers up to the isotope layer
    """
    return (formula_sublayer(ich1) == formula_sublayer(ich2) and
            main_sublayers(ich1) == main_sublayers(ich2) and
            charge_sublayers(ich1) == charge_sublayers(ich2) and
            stereo_sublayers(ich1) == stereo_sublayers(ich2) and
            isotope_sublayers(ich1) == isotope_sublayers(ich2))


def same_connectivity(ich1, ich2):
    """ do these inchis have the same connectivity
    """
    return without_stereo(ich1) == without_stereo(ich2)


# transformations/operations
def join(ichs):
    """ join separate inchis into one multi-component inchi
    """
    fml_slyr = '.'.join(map(formula_sublayer, ichs))
    main_dct = _join_sublayers(list(map(main_sublayers, ichs)))
    char_dct = _join_sublayers(list(map(charge_sublayers, ichs)))
    ste_dct = _join_sublayers(list(map(stereo_sublayers, ichs)))
    iso_dct = _join_sublayers(list(map(isotope_sublayers, ichs)))
    return from_data(fml_slyr=fml_slyr, main_dct=main_dct, char_dct=char_dct,
                     ste_dct=ste_dct, iso_dct=iso_dct)


def _join_sublayers(dcts):
    """ join sublayer dictionaries into one multi-component sublayer dictionary
    """
    pfxs = functools.reduce(set.union, map(set, dcts))
    dcts = [dict_.by_key(dct, pfxs, fill_val='') for dct in dcts]
    dct = {pfx: ';'.join([dct[pfx] for dct in dcts]) for pfx in pfxs}
    return dct


def split(ich):
    """ split a multi-component inchi into inchis for each of its components
    """
    fml_slyr = formula_sublayer(ich)
    main_dct = main_sublayers(ich)
    char_dct = charge_sublayers(ich)
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    fml_slyrs = _split_sublayer_string(
        fml_slyr, count_sep_ptt='', sep_ptt=app.escape('.'))
    count = len(fml_slyrs)

    main_dcts = _split_sublayers(main_dct, count)
    char_dcts = _split_sublayers(char_dct, count)
    ste_dcts = _split_sublayers(ste_dct, count)
    iso_dcts = _split_sublayers(iso_dct, count)

    ichs = tuple(from_data(fml_slyr=fml_slyr, main_dct=main_dct,
                           char_dct=char_dct, ste_dct=ste_dct, iso_dct=iso_dct)
                 for fml_slyr, main_dct, char_dct, ste_dct, iso_dct
                 in zip(fml_slyrs, main_dcts, char_dcts, ste_dcts, iso_dcts))
    return ichs


def _split_sublayers(dct, count):
    """ split a multi-component sublayer dictionary into separate ones
    """
    if dct:
        pfxs = sorted(dct.keys())
        slyrs_lst = list(
            map(_split_sublayer_string, dict_.values_by_key(dct, pfxs)))
        assert all(len(slyrs) == count for slyrs in slyrs_lst)
        dcts = tuple({pfx: slyr for pfx, slyr in zip(pfxs, slyrs) if slyr}
                     for slyrs in zip(*slyrs_lst))
    else:
        return (dict(),) * count
    return dcts


def _split_sublayer_string(lyr, count_sep_ptt=app.escape('*'),
                           sep_ptt=app.escape(';')):
    count_ptt = app.UNSIGNED_INTEGER
    group_ptt = (app.STRING_START + app.capturing(count_ptt) + count_sep_ptt +
                 app.capturing(app.one_or_more(app.WILDCARD)))

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


# conversions
def recalculate(ich, force_stereo=False):
    """ recalculate an inchi string
    """
    fml_slyr = formula_sublayer(ich)
    main_dct = main_sublayers(ich)
    char_dct = charge_sublayers(ich)
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    ich = from_data(fml_slyr, main_dct=main_dct, char_dct=char_dct,
                    ste_dct=ste_dct, iso_dct=iso_dct)
    return automol.convert.inchi.recalculate(ich, force_stereo=force_stereo)


def geometry(ich):
    """ inchi => geometry
    """
    return automol.convert.inchi.geometry(ich)


def graph(ich, no_stereo=False):
    """ inchi => graph
    """
    return automol.convert.inchi.graph(ich, no_stereo=no_stereo)


def smiles(ich):
    """ inchi => smiles
    """
    return automol.convert.inchi.smiles(ich)


def inchi_key(ich):
    """ inchi => inchi_key
    """
    return automol.convert.inchi.inchi_key(ich)


def formula(ich):
    """ inchi => formula
    """
    return automol.convert.inchi.formula(ich)


def _sublayers(lyr):
    """ get sublayers from a layer, by prefix
    """
    if lyr:
        ptt = _sublayer_pattern(key_ptt=app.capturing(app.LOWERCASE_LETTER),
                                val_ptt=app.capturing(NONSLASHES))
        dct = dict(apf.all_captures(ptt, lyr))
    else:
        dct = dict()
    return dct


def _main_layer(ich):
    """ main layer
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           SLASH + app.capturing(_main_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _charge_layer(ich):
    """ charge layer
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           SLASH + app.capturing(_charge_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _stereo_layer(ich):
    """ stereo layer
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           app.maybe(SLASH + _charge_layer_pattern()) +
           SLASH + app.capturing(_stereo_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _isotope_layer(ich):
    """ isotope layer
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           app.maybe(SLASH + _charge_layer_pattern()) +
           app.maybe(SLASH + _stereo_layer_pattern()) +
           SLASH + app.capturing(_isotope_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _isotope_layer_pattern():
    i_slyr_ptt = _sublayer_pattern(key_ptt='i')
    h_slyr_ptt = _sublayer_pattern(key_ptt='h')
    ptt = (i_slyr_ptt +
           app.maybe(SLASH + h_slyr_ptt) +
           app.maybe(SLASH + _stereo_layer_pattern()))
    return ptt


def _stereo_layer_pattern():
    b_slyr_ptt = _sublayer_pattern(key_ptt='b')
    t_slyr_ptt = _sublayer_pattern(key_ptt='t')
    m_slyr_ptt = _sublayer_pattern(key_ptt='m')
    s_slyr_ptt = _sublayer_pattern(key_ptt='s')
    ptt = (app.one_of_these([b_slyr_ptt, t_slyr_ptt]) +
           app.maybe(SLASH + t_slyr_ptt) +
           app.maybe(SLASH + m_slyr_ptt) +
           app.maybe(SLASH + s_slyr_ptt))
    return ptt


def _charge_layer_pattern():
    q_slyr_ptt = _sublayer_pattern(key_ptt='q')
    p_slyr_ptt = _sublayer_pattern(key_ptt='p')
    ptt = (app.one_of_these([q_slyr_ptt, p_slyr_ptt]) +
           app.maybe(SLASH + p_slyr_ptt))
    return ptt


def _main_layer_pattern():
    c_slyr_ptt = _sublayer_pattern(key_ptt='c')
    h_slyr_ptt = _sublayer_pattern(key_ptt='h')
    ptt = (app.one_of_these([c_slyr_ptt, h_slyr_ptt]) +
           app.maybe(SLASH + h_slyr_ptt))
    return ptt


def _formula_sublayer_pattern():
    ptt = _sublayer_pattern(key_ptt=app.not_followed_by(app.LOWERCASE_LETTER))
    return ptt


def _version_pattern():
    ptt = app.preceded_by('InChI=') + _sublayer_pattern()
    return ptt


def _sublayer_pattern(key_ptt='',
                      val_ptt=NONSLASHES):
    """ inchi sublayer pattern
    """
    return key_ptt + val_ptt
