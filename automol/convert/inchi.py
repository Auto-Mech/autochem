""" inchi conversions
"""

import functools
import itertools
import operator
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from phydat import phycon
from automol.util import dict_
import automol.formula
from automol.convert import _rdkit


NONSLASH = '[^/]'
NONSLASHES = app.one_or_more(NONSLASH)
SLASH = app.escape('/')
SLASH_OR_START = app.one_of_these([SLASH, app.STRING_START])
SLASH_OR_END = app.one_of_these([SLASH, app.STRING_END])


def low_spin_multiplicity(ich):
    """ Guess spin multiplicity based on the number of electrons.

        :param ich: InChI string
        :type ich: str
        :rtype: int
    """

    fml = formula(ich)
    nelec = automol.formula.electron_count(fml)

    if (nelec % 2) == 0:
        mult = 1
    else:
        mult = 2

    return mult


def recalculate(ich, stereo=False):
    """ Recalculate an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: force the same stereochem in recalculated InChI
        :type stereo: bool
        :rtype: str
    """

    # for now, just assert that we have no multi-component strings with
    # hardcoded parts -- these are guaranteed to fail
    ichs = split(ich)
    if len(ichs) > 1:
        if any(object_from_hardcoded_inchi_by_key('inchi', ich)
               for ich in ichs):
            ref_ichs = []
            for ich_i in ichs:
                ref_ichs.append(recalculate(ich_i))
            ref_ichs.sort()
            ret = join(ref_ichs)
            return ret
        # raise error.FailedInchiGenerationError

    ret = object_from_hardcoded_inchi_by_key('inchi', ich)
    if ret is None:
        _options = '-SUU' if stereo else ''
        rdm = _rdkit.from_inchi(ich)
        ret = _rdkit.to_inchi(rdm, options=_options, with_aux_info=False)

    return ret


def smiles(ich):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = split(ich)
    smis = list(map(_connected_smiles, ichs))
    smi = '.'.join(smis)
    return smi


def _connected_smiles(ich):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """

    smi = object_from_hardcoded_inchi_by_key('smiles', ich)
    if smi is None:
        ich = standard_form(ich)
        rdm = _rdkit.from_inchi(ich)
        smi = _rdkit.to_smiles(rdm)

    return smi


def inchi_key(ich):
    """ Generate an InChIKey from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return _rdkit.inchi_to_inchi_key(ich)


def formula(ich):
    """ Generate a formula dictionary from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: int]
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = split(ich)
    fmls = list(map(_connected_formula, ichs))
    fml = functools.reduce(automol.formula.join, fmls)

    return fml


def _connected_formula(ich):
    """ Create a combined molecular from the formulas of a
        multi-component InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: int]
    """

    fml = object_from_hardcoded_inchi_by_key('formula', ich)
    if fml is None:
        ich = standard_form(ich)
        rdm = _rdkit.from_inchi(ich)
        fml = _rdkit.to_formula(rdm)

    return fml


# hardcoded inchis which neither RDKit nor Pybel can handle
HARDCODED_INCHI_DCT = {
    'InChI=1S/C': {
        'inchi': 'InChI=1S/C',
        'geom': (('C', (0., 0., 0.)),),
        'graph': ({0: ('C', 0, None)}, {}),
        'smiles': '[C]',
        'formula': {'C': 1},
    },
    'InChI=1S/B': {
        'inchi': 'InChI=1S/B',
        'geom': (('B', (0., 0., 0.)),),
        'graph': ({0: ('B', 0, None)}, {}),
        'smiles': '[B]',
        'formula': {'B': 1},
    },
    'InChI=1S/N': {
        'inchi': 'InChI=1S/N',
        'geom': (('N', (0., 0., 0.)),),
        'graph': ({0: ('N', 0, None)}, {}),
        'smiles': '[N]',
        'formula': {'N': 1},
    },
    'InChI=1S/CH/h1H': {
        'inchi': 'InChI=1S/CH/h1H',
        'geom': (('C', (0., 0., 0.)),
                 ('H', (0., 0., 1.12 * phycon.ANG2BOHR))),
        'graph': ({0: ('C', 1, None)}, {}),
        'smiles': '[CH]',
        'formula': {'C': 1, 'H': 1},
    },
    'InChI=1S/CF/c1-2': {
        'inchi': 'InChI=1S/CF/c1-2',
        'geom': (('C', (0., 0., 0.)),
                 ('F', (0., 0., 1.27 * phycon.ANG2BOHR))),
        'graph': ({0: ('C', 0, None), 1: ('F', 0, None)},
                  {frozenset({0, 1}): (1, None)}),
        'smiles': '[C]F',
        'formula': {'C': 1, 'F': 1},
    },
    'InChI=1S/CCl/c1-2': {
        'inchi': 'InChI=1S/CCl/c1-2',
        'geom': (('C', (0., 0., 0.)),
                 ('Cl', (0., 0., 1.65 * phycon.ANG2BOHR))),
        'graph': ({0: ('C', 0, None), 1: ('Cl', 0, None)},
                  {frozenset({0, 1}): (1, None)}),
        'smiles': '[C]Cl',
        'formula': {'C': 1, 'Cl': 1},
    },
    'InChI=1S/CBr/c1-2': {
        'inchi': 'InChI=1S/CBr/c1-2',
        'geom': (('C', (0., 0., 0.)),
                 ('Br', (0., 0., 1.8 * phycon.ANG2BOHR))),
        'graph': ({0: ('C', 0, None), 1: ('Br', 0, None)},
                  {frozenset({0, 1}): (1, None)}),
        'smiles': '[C]Br',
        'formula': {'C': 1, 'Br': 1},
    },
    'InChI=1S/CI/c1-2': {
        'inchi': 'InChI=1S/CI/c1-2',
        'geom': (('C', (0., 0., 0.)),
                 ('I', (0., 0., 1.8 * phycon.ANG2BOHR))),
        'graph': ({0: ('C', 0, None), 1: ('I', 0, None)},
                  {frozenset({0, 1}): (1, None)}),
        'smiles': '[C]I',
        'formula': {'C': 1, 'I': 1},
    },
}


def object_from_hardcoded_inchi_by_key(key, ich):
    """ Obtains the requested structural identifier object
        for certain hardcoded InChI string.

        InChI strings: C, B, N, CH, CF, CCl, CBr, CI

        :param key: key for structural identifier
        :type key: str
        :param ich: InChI string
        :type ich: str
        :rtype: obj
    """

    obj = None
    for ich_, obj_dct in HARDCODED_INCHI_DCT.items():
        if equivalent(ich, ich_):
            obj = obj_dct[key]
    return obj


def object_to_hardcoded_inchi_by_key(key, obj, comp=operator.eq):
    """ Convert a structural identifier to an InChI string object if that
        InChI <=> relation is hardoded in automol.

        InChI strings: C, B, N, CH, CF, CCl, CBr, CI

        :param key: key for structural identifier
        :type key: str
        :param obj: obj for structural identifier
        :type obj: str
        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ich = None
    for ich_, obj_dct in HARDCODED_INCHI_DCT.items():
        obj_ = obj_dct[key]
        if comp(obj, obj_):
            ich = ich_
    return ich


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


def formula_sublayer(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    ptt = (version_pattern() +
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
    ptt = (version_pattern() +
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
    ptt = (version_pattern() +
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
    ptt = (version_pattern() +
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
    ptt = (version_pattern() +
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


def version_pattern():
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
    fmls = list(map(automol.convert.inchi.formula_sublayer, ichs))
    fml_slyr = _join_sublayer_strings(fmls, count_sep='', sep='.')
    main_dct = _join_sublayers(list(map(
        automol.convert.inchi.main_sublayers, ichs)))
    char_dct = _join_sublayers(list(map(
        automol.convert.inchi.charge_sublayers, ichs)))
    ste_dct = _join_sublayers(list(map(
        automol.convert.inchi.stereo_sublayers, ichs)))
    iso_dct = _join_sublayers(list(map(
        automol.convert.inchi.isotope_sublayers, ichs)))

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
