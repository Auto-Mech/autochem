""" InChI chemical identifiers
"""

import itertools
import functools
import numpy
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from automol.util import dict_
import automol.convert.geom
import automol.convert.inchi


NONSLASH = '[^/]'
NONSLASHES = app.one_or_more(NONSLASH)
SLASH = app.escape('/')
SLASH_OR_START = app.one_of_these([SLASH, app.STRING_START])
SLASH_OR_END = app.one_of_these([SLASH, app.STRING_END])


# "constructor"
def from_data(fml_slyr, main_lyr_dct=None, char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None):
    """ Construct an InChI string from its constituent layers, where
        most layers are input as dictionary of prefixes for some part of the
        sublayer and the corresponding string for that sublayer part.

        :param fml_slyr: sublayer of InChI string containing molecular formula
        :type fml_slyr: str
        :param main_lyr_dct: information for connectivity layer of InChI
        :type main_lyr_dct: dict[str: str]
        :param char_lyr_dct: information for charge layer of InChI
        :type char_lyr_dct: dict[str: str]
        :param ste_lyr_dct: information for stereochemistry layer of InChI
        :type ste_lyr_dct: dict[str: str]
        :param iso_lyr_dct: information for isotope layer of InChI
        :type iso_lyr_dct: dict[str: str]
        :rtype: str
    """
    return automol.create.inchi.from_data(
        fml_slyr=fml_slyr, main_lyr_dct=main_lyr_dct,
        char_lyr_dct=char_lyr_dct, ste_lyr_dct=ste_lyr_dct,
        iso_lyr_dct=iso_lyr_dct)


# getters
def version(ich):
    """ Determine version of InChI the string corresponds to.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = app.capturing(_version_pattern())
    ver = apf.first_capture(ptt, ich)
    return ver


def formula_sublayer(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    ptt = (_version_pattern() +
           SLASH + app.capturing(_formula_sublayer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def main_sublayers(ich):
    """ Parse the InChI string for the sublayers of the connectivity layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return _sublayers(_main_layer(ich))


def charge_sublayers(ich):
    """ Parse the InChI string for the sublayers of the charge layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return _sublayers(_charge_layer(ich))


def stereo_sublayers(ich):
    """ Parse the InChI string for the sublayers of the stereochemisty layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return _sublayers(_stereo_layer(ich))


def isotope_sublayers(ich):
    """ Parse the InChI string for the sublayers of the isotope layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return _sublayers(_isotope_layer(ich))


# setters
def standard_form(ich, stereo=True):
    """ Return an InChI string in standard form.

        Eventually we should just designate standard-form as standard InChI
        ordering for all but the hardcoded exceptions, put at the end.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """

    fml_slyr = formula_sublayer(ich)
    main_dct = main_sublayers(ich)
    char_dct = charge_sublayers(ich)

    if stereo:
        ste_dct = stereo_sublayers(ich)
        iso_dct = isotope_sublayers(ich)
    else:
        ste_dct = {}
        iso_dct = dict_.by_key(isotope_sublayers(ich),
                               automol.create.inchi.ISO_NONSTE_PFXS)

    ich = from_data(fml_slyr,
                    main_lyr_dct=main_dct,
                    char_lyr_dct=char_dct,
                    ste_lyr_dct=ste_dct,
                    iso_lyr_dct=iso_dct)

    return recalculate(ich)


def has_stereo(ich):
    """ Determine if the InChI string has stereochemistry information.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return bool(ste_dct or
                any(pfx in iso_dct for pfx in automol.create.inchi.STE_PFXS))


def has_multiple_components(ich):
    """ Determine if the InChI string has multiple components.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return len(split(ich)) > 1


def is_standard_form(ich):
    """ Determine if the InChI string is closed.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return ich == standard_form(ich)


def is_complete(ich):
    """ Determine if the InChI string is complete
        (has all stereo-centers assigned).

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return equivalent(ich, standard_form(ich)) and not (
        has_stereo(ich) ^ has_stereo(recalculate(ich, stereo=True)))


def is_chiral(ich):
    """ Determine if the InChI string has chirality information.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return ste_dct['s'] == '1' or iso_dct['s'] == '1'


# comparisons
def equivalent(ich1, ich2):
    """ Determine if two InChI strings are equivalent. Currently
        the srings are only checked up to the isotope sublayer.

        :param ich1: InChI string 1
        :type ich1: str
        :param ich2: InChI string 2
        :type ich2: str
        :rtype: bool
    """
    return (formula_sublayer(ich1) == formula_sublayer(ich2) and
            main_sublayers(ich1) == main_sublayers(ich2) and
            charge_sublayers(ich1) == charge_sublayers(ich2) and
            stereo_sublayers(ich1) == stereo_sublayers(ich2) and
            isotope_sublayers(ich1) == isotope_sublayers(ich2))


def same_connectivity(ich1, ich2):
    """ Determine if two InChI strings have the same connectivity.

        :param ich1: InChI string 1
        :type ich1: str
        :param ich2: InChI string 2
        :type ich2: str
        :rtype: bool
    """
    return (standard_form(ich1, stereo=False) ==
            standard_form(ich2, stereo=False))


def sorted_(ichs):
    """ Sort a sequence of InChI strings in their standard form sort order.

        :param ichs: sequence of InChI strings
        :type ichs: tuple(str)
        :rtype: tuple(str)
    """
    return tuple(ichs[idx] for idx in argsort(ichs))


def argsort(ichs):
    """ Determine the sort order for the InChI standard form.

        :param ichs: sequence of InChI strings
        :type ichs: tuple(str)
    """

    assert not any(map(has_multiple_components, ichs))
    ref_ichs = list(map(standard_form, split(recalculate(join(ichs)))))
    idxs = tuple(numpy.argsort(list(map(ref_ichs.index, ichs))))
    return idxs


# transformations/operations
def join(ichs):
    """ Join separate InChI strings into one multi-component InChI string.

        Currently:
        (fix for /s [which should be removed in split/join operations] and /m,
         which is joined as /m0110..  with no separators).

        :param ichs: sequence of InChI strings
        :type ichs: tuple(str)
        :rtype: str
    """

    # first, make sure they are completely split up
    ichs = list(itertools.chain(*map(split, ichs)))
    fmls = list(map(formula_sublayer, ichs))
    fml_slyr = _join_sublayer_strings(fmls, count_sep='', sep='.')
    main_dct = _join_sublayers(list(map(main_sublayers, ichs)))
    char_dct = _join_sublayers(list(map(charge_sublayers, ichs)))
    ste_dct = _join_sublayers(list(map(stereo_sublayers, ichs)))
    iso_dct = _join_sublayers(list(map(isotope_sublayers, ichs)))

    return from_data(fml_slyr=fml_slyr,
                     main_lyr_dct=main_dct,
                     char_lyr_dct=char_dct,
                     ste_lyr_dct=ste_dct,
                     iso_lyr_dct=iso_dct)


def _join_sublayers(dcts):
    """ Join all of the components of an InChI sublayer.

        :param dcts: sublayer components, grouped by prefix
        :type dct: dict[str: str]
        :rtype: dict[str: str]
    """

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
    m_slyrs = [m_slyr if m_slyr else '.' for m_slyr in m_slyrs]
    return ''.join(m_slyrs)


def _join_sublayer_strings(slyrs, count_sep='*', sep=';'):
    """ Join sublayer strings into one multi-component sublayer string.

        :param slyrs: sublayers to join
        :type slyrs: tuple(str)?
        :param count_sep: delimiter for ???
        :type count_sep: str
        :param sep: delimiter for ???
        :type sep: str
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
    """ Split a multi-component InChI into InChIs for each of its components.

        (fix this for /s [which should be removed in split/join operations]
         and /m, which is joined as /m0110..  with no separators)

        :param ich: InChI string
        :type ich: str
        :rtype: tuple(str)
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

    ichs = tuple(from_data(fml_slyr=fml_slyr,
                           main_lyr_dct=main_dct,
                           char_lyr_dct=char_dct,
                           ste_lyr_dct=ste_dct,
                           iso_lyr_dct=iso_dct)
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
def low_spin_multiplicity(ich):
    """ Guess spin multiplicity based on the number of electrons.

        :param ich: InChI string
        :type ich: str
        :rtype: int
    """
    return automol.convert.inchi.low_spin_multiplicity(ich)

def recalculate(ich, stereo=False):
    """ Recalculate an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: str
    """
    return automol.convert.inchi.recalculate(ich, stereo=stereo)


def geometry(ich):
    """ Generate a molecular geometry from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: automol geometry data structure
    """
    return automol.convert.inchi.geometry(ich)


def conformers(ich, nconfs=100):
    """ Generate a molecular geometry for many conformers from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param nconfs: number of conformers to generate
        :type nconfs: int
        :rtype: tuple(automol geometry data structure)
    """
    return automol.convert.inchi.conformers(ich, nconfs)


def graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol graph data structure
    """
    return automol.convert.inchi.graph(ich, stereo=stereo)


def smiles(ich):
    """ Generate a corresponding SMILES string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return automol.convert.inchi.smiles(ich)


def inchi_key(ich):
    """ Generate an InChI key (what?) from an InChI string.

        :param ich: InChI string
        :type ich: str
    """
    return automol.convert.inchi.inchi_key(ich)


def formula(ich):
    """ Generate a formula dictionary from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: int]
    """
    return automol.convert.inchi.formula(ich)


def formula_string(ich):
    """ Generate a formula string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return formula_sublayer(ich)


def add_stereo(ich):
    """ Add stereochemistry to an InChI string converting to/from geometry.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    geo = automol.convert.inchi.geometry(ich)
    ich = automol.convert.geom.inchi(geo, stereo=True)
    return ich


def expand_stereo(ich):
    """ Obtain all possible stereoisomers compatible with an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: list[str]
    """
    gra = automol.inchi.graph(ich)
    stereo_gras = automol.graph.stereomers(gra)

    stereo_ichs = []
    for gra in stereo_gras:
        stereo_ichs.append(automol.graph.inchi(gra, stereo=True))
    return stereo_ichs
    # return list(map(automol.graph.inchi, stereo_gras))


def _sublayers(lyr):
    """ Parse the sublayers of the specified layer of an InChI string,
        organized by prefix.

        :param lyr: layer of the InChI string
        :type lyr: str
        :rtype: dict[str: str]
    """
    if lyr:
        ptt = _sublayer_pattern(key_ptt=app.capturing(app.LOWERCASE_LETTER),
                                val_ptt=app.capturing(NONSLASHES))
        dct = dict(apf.all_captures(ptt, lyr))
    else:
        dct = dict()
    return dct


def _main_layer(ich):
    """ Parse the InChI string for the connectivity layer.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           SLASH + app.capturing(_main_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _charge_layer(ich):
    """ Parse the InChI string for the charge layer.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           SLASH + app.capturing(_charge_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _stereo_layer(ich):
    """ Parse the InChI string for the stereochemisty layer.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           app.maybe(SLASH + _charge_layer_pattern()) +
           SLASH + app.capturing(_stereo_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _isotope_layer(ich):
    """ Parse the InChI string for the isotope layer.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = (_version_pattern() +
           SLASH + _formula_sublayer_pattern() +
           app.maybe(SLASH + _main_layer_pattern()) +
           app.maybe(SLASH + _charge_layer_pattern()) +
           app.maybe(SLASH + _stereo_layer_pattern()) +
           SLASH + app.capturing(_isotope_layer_pattern()))
    lyr = apf.first_capture(ptt, ich)
    return lyr


def _formula_sublayer_pattern():
    """ Build the autoparse regex pattern for the chemical formual sublayer.

        :rtype: str
    """
    ptt = _sublayer_pattern(key_ptt=app.not_followed_by(app.LOWERCASE_LETTER))
    return ptt


def _main_layer_pattern():
    """ Build the autoparse regex pattern for the connectivity layer.

        :rtype: str
    """
    c_slyr_ptt = _sublayer_pattern(key_ptt='c')
    h_slyr_ptt = _sublayer_pattern(key_ptt='h')
    ptt = (app.one_of_these([c_slyr_ptt, h_slyr_ptt]) +
           app.maybe(SLASH + h_slyr_ptt))
    return ptt


def _charge_layer_pattern():
    """ Build the autoparse regex pattern for the charge layer.

        :rtype: str
    """
    q_slyr_ptt = _sublayer_pattern(key_ptt='q')
    p_slyr_ptt = _sublayer_pattern(key_ptt='p')
    ptt = (app.one_of_these([q_slyr_ptt, p_slyr_ptt]) +
           app.maybe(SLASH + p_slyr_ptt))
    return ptt


def _stereo_layer_pattern():
    """ Build the autoparse regex pattern for the stereochemistry layer.

        :rtype: str
    """
    b_slyr_ptt = _sublayer_pattern(key_ptt='b')
    t_slyr_ptt = _sublayer_pattern(key_ptt='t')
    m_slyr_ptt = _sublayer_pattern(key_ptt='m')
    s_slyr_ptt = _sublayer_pattern(key_ptt='s')
    ptt = (app.one_of_these([b_slyr_ptt, t_slyr_ptt]) +
           app.maybe(SLASH + t_slyr_ptt) +
           app.maybe(SLASH + m_slyr_ptt) +
           app.maybe(SLASH + s_slyr_ptt))
    return ptt


def _isotope_layer_pattern():
    """ Build the autoparse regex pattern for the isotope layer.

        :rtype: str
    """
    i_slyr_ptt = _sublayer_pattern(key_ptt='i')
    h_slyr_ptt = _sublayer_pattern(key_ptt='h')
    ptt = (i_slyr_ptt +
           app.maybe(SLASH + h_slyr_ptt) +
           app.maybe(SLASH + _stereo_layer_pattern()))
    return ptt


def _version_pattern():
    """ Build the autoparse regex pattern for the InChI string version.

        :rtype: str
    """
    ptt = app.preceded_by('InChI=') + _sublayer_pattern()
    return ptt


def _sublayer_pattern(key_ptt='',
                      val_ptt=NONSLASHES):
    """ Build the autoparse regex pattern for an arbitrary InChI sublayer.

        :rtype: str
    """
    return key_ptt + val_ptt
