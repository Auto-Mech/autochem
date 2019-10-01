""" InChI chemical identifiers
"""
import itertools
import functools
import numpy
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


# "constructor"
def from_data(fml_slyr, main_dct=None, char_dct=None, ste_dct=None,
              iso_dct=None):
    """ calculate an inchi string from layers
    """
    return automol.create.inchi.from_data(
        formula_sublayer=fml_slyr, main_sublayer_dct=main_dct,
        charge_sublayer_dct=char_dct, stereo_sublayer_dct=ste_dct,
        isotope_sublayer_dct=iso_dct)


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
def standard_form(ich, remove_stereo=False):
    """ return an inchi string in standard form

    (eventually we should just designate standard-form as standard inchi
    ordering for all but the hardcoded exceptions, which we can put at the end)
    """
    if remove_stereo:
        fml_slyr = formula_sublayer(ich)
        main_dct = main_sublayers(ich)
        char_dct = charge_sublayers(ich)
        ste_dct = {}
        iso_dct = dict_.by_key(isotope_sublayers(ich),
                               automol.create.inchi.ISO_NONSTE_PFXS)
    else:
        fml_slyr = formula_sublayer(ich)
        main_dct = main_sublayers(ich)
        char_dct = charge_sublayers(ich)
        ste_dct = stereo_sublayers(ich)
        iso_dct = isotope_sublayers(ich)

    ich = from_data(fml_slyr, main_dct=main_dct, char_dct=char_dct,
                    ste_dct=ste_dct, iso_dct=iso_dct)
    return recalculate(ich)


def has_stereo(ich):
    """ does this inchi have stereo information?
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return bool(ste_dct or
                any(pfx in iso_dct for pfx in automol.create.inchi.STE_PFXS))


def has_multiple_components(ich):
    """ does this inchi have multiple components?
    """
    return len(split(ich)) > 1


def is_standard_form(ich):
    """ is this inchi closed?
    """
    return ich == standard_form(ich)


def is_complete(ich):
    """ is this inchi complete? (are all stereo-centers assigned?)
    """
    return equivalent(ich, standard_form(ich)) and not (
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
    return (standard_form(ich1, remove_stereo=True) ==
            standard_form(ich2, remove_stereo=True))


def sorted_(ichs):
    """ sort the inchies in their standard form sort order
    """
    return tuple(ichs[idx] for idx in argsort(ichs))


def argsort(ichs):
    """ determine the sort order for the inchi standard form
    """

    assert not any(map(has_multiple_components, ichs))
    current_ichs = ichs #list(map(standard_form, ichs))
    ref_ichs = list(map(standard_form, split(recalculate(join(ichs)))))
    #ref_ichs = list(map(standard_form, split(recalculate(join(ichs)))))
    idxs = tuple(numpy.argsort(list(map(ref_ichs.index, ichs))))
    return idxs


# transformations/operations
def join(ichs):
    """ join separate inchis into one multi-component inchi

    (fix this for /s [which should be removed in split/join operations] and /m,
    which is joined as /m0110..  with no separators)
    """
    # first, make sure they are completely split up
    ichs = list(itertools.chain(*map(split, ichs)))
    fmls = list(map(formula_sublayer, ichs))
    fml_slyr = _join_sublayer_strings(fmls, count_sep='', sep='.')
    main_dct = _join_sublayers(list(map(main_sublayers, ichs)))
    char_dct = _join_sublayers(list(map(charge_sublayers, ichs)))
    ste_dct = _join_sublayers(list(map(stereo_sublayers, ichs)))
    iso_dct = _join_sublayers(list(map(isotope_sublayers, ichs)))
    return from_data(fml_slyr=fml_slyr, main_dct=main_dct, char_dct=char_dct,
                     ste_dct=ste_dct, iso_dct=iso_dct)


def _join_sublayers(dcts):
    pfxs = sorted(functools.reduce(set.union, map(set, dcts)))
    if 's' in pfxs:
        pfxs.remove('s')
    dcts = [dict_.by_key(dct, pfxs, fill_val='') for dct in dcts]
    slyrs_lst = [[dct[pfx] for dct in dcts] for pfx in pfxs]
    dct = {pfx: (_join_sublayer_strings(slyrs) if pfx != 'm' else
                 _join_m_sublayer_strings(slyrs))
           for pfx, slyrs in zip(pfxs, slyrs_lst)}
    return dct


def _join_m_sublayer_strings(m_slyrs):
    return ''.join(m_slyrs)


def _join_sublayer_strings(slyrs, count_sep='*', sep=';'):
    """ join sublayer strings into one multi-component sublayer string
    """
    def _s(count, slyr):
        if count > 1 and slyr:
            ret = ('{:d}' + count_sep + '{:s}').format(count, slyr)
        elif slyr:
            ret = slyr
        else:
            ret = sep * (count - 1)
        return ret

    counts, slyrs = zip(*[
        (len(list(g)), slyr) for slyr, g in itertools.groupby(slyrs)])

    slyr = sep.join([_s(count, slyr) for count, slyr in zip(counts, slyrs)])
    return slyr


def split(ich):
    """ split a multi-component inchi into inchis for each of its components

    (fix this for /s [which should be removed in split/join operations] and /m,
    which is joined as /m0110..  with no separators)
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
        if 's' in pfxs:
            pfxs.remove('s')
        slyrs_lst = [
            _split_sublayer_string(dct[pfx]) if pfx != 'm'
            else _split_m_sublayer_string(dct[pfx]) for pfx in pfxs]
        assert all(len(slyrs) == count for slyrs in slyrs_lst)
        dcts = tuple({pfx: slyr for pfx, slyr in zip(pfxs, slyrs) if slyr}
                     for slyrs in zip(*slyrs_lst))
    else:
        return (dict(),) * count
    return dcts


def _split_m_sublayer_string(m_slyr):
    return tuple(m_slyr)


def _split_sublayer_string(slyr, count_sep_ptt=app.escape('*'),
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
        itertools.chain(*map(_expand_group, apf.split(sep_ptt, slyr))))
    return parts


# conversions
def recalculate(ich, force_stereo=False):
    """ recalculate an inchi string
    """
    return automol.convert.inchi.recalculate(ich, force_stereo=force_stereo)


def geometry(ich):
    """ inchi => geometry
    """
    return automol.convert.inchi.geometry(ich)


def conformers(ich, nconfs=100):
    """ inchi => conformers
    """
    return automol.convert.inchi.conformers(ich, nconfs)


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

def add_stereo(ich):
    """ inchi => inchi after adding stereoinfo
    """
    gra = automol.inchi.graph(ich)
    stereo_gras = automol.graph.stereomers(gra)
    return list(map(lambda x: automol.graph.inchi(x), stereo_gras))

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
