""" Level 3 InChI functions (depend on extern and L1-2)
"""

import operator
import functools
import itertools
import numpy
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
from phydat import phycon
import automol.formula
from automol.util import dict_
from automol.extern import rdkit_


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
    'InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5+': {
        'inchi': 'InChI=1S/C7H13/c1-6(2)5-7(3)4/h5,7H,1H2,2-4H3/b6-5+',
        'geom': (('C', (2.9504879811871962, 2.292895414872005, -1.3368243864930003)), ('C', (2.5673838037731307, -0.008139050422346866, -0.17152288150330522)), ('H', (1.3969346848176885, 3.447419811768214, -2.0135183044843195)), ('H', (4.841962103311073, 3.0346431547322097, -1.6187582963284302)), ('C', (4.829134642371465, -1.49268143841626, 0.7751467594015465)), ('C', (0.16754878746146742, -1.0464660775902797, 0.22737940631958772)), ('H', (6.602497102201977, -0.7242545451168418, 0.04279095838486706)), ('H', (4.710550548546736, -3.478611631195021, 0.2053414202444985)), ('H', (4.926744665929738, -1.4471654949584831, 2.8436504249583128)), ('C', (-2.309165956812197, 0.1601694069415546, -0.45936407520691613)), ('H', (0.059305274995243025, -2.84879803892406, 1.2211693681627305)), ('C', (-4.338569299107115, -1.849291655551407, -0.9380128055241295)), ('C', (-3.1609089852569254, 2.0152228374494827, 1.6138601262112435)), ('H', (-2.064962318797658, 1.2437024624826896, -2.2124932374121506)), ('H', (-4.671871634210865, -2.974326876259847, 0.7713597482461291)), ('H', (-6.134982638215961, -0.9730483278901176, -1.4766508916940015)), ('H', (-3.7713018621365566, -3.142000385645594, -2.450961556635925)), ('H', (-1.7341147379308777, 3.4801498682611434, 1.9233556915864756)), ('H', (-4.944132036010076, 2.9327283350601436, 1.0936960026438471)), ('H', (-3.4520912147503457, 1.0193031542629487, 3.407049592550058))),
        'graph': ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None), 3: ('H', 0, None), 4: ('C', 0, None), 5: ('C', 0, None), 6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None), 9: ('C', 0, None), 10: ('H', 0, None), 11: ('C', 0, None), 12: ('C', 0, None), 13: ('H', 0, None), 14: ('H', 0, None), 15: ('H', 0, None), 16: ('H', 0, None), 17: ('H', 0, None), 18: ('H', 0, None), 19: ('H', 0, None)}, {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None), frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None), frozenset({1, 5}): (1, False), frozenset({4, 6}): (1, None), frozenset({4, 7}): (1, None), frozenset({8, 4}): (1, None), frozenset({9, 5}): (1, None), frozenset({10, 5}): (1, None), frozenset({9, 11}): (1, None), frozenset({9, 12}): (1, None), frozenset({9, 13}): (1, None), frozenset({11, 14}): (1, None), frozenset({11, 15}): (1, None), frozenset({16, 11}): (1, None), frozenset({17, 12}): (1, None), frozenset({18, 12}): (1, None), frozenset({19, 12}): (1, None)}),
        'smiles': '[CH2]C(C)=CC(C)C',
        'formula': {'C': 7, 'H': 13},
    },
    'InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5+': {
        'inchi': 'InChI=1S/C7H13/c1-5-6-7(2,3)4/h5-6H,1H2,2-4H3/b6-5+',
        'geom': (('C', (5.753492286913037, -0.37920000323886965, 0.00022487740892948157)), ('C', (3.3144889171990086, 0.5684220596332129, 0.00027778974044230075)), ('H', (7.3841936523543605, 0.8619853440447102, -0.00043841646110621616)), ('H', (6.105984790820811, -2.402233078764869, 0.0007899055204413722)), ('C', (1.1576991367869782, -0.9122161541855046, 0.0014343021292224916)), ('H', (3.096112166141102, 2.6144833377240415, -0.00012094247202930099)), ('C', (-1.55486287657445, -0.03241447222997813, 0.0001285013765311323)), ('H', (1.4091309772313934, -2.9602427474468094, 0.0018746083164541657)), ('C', (-2.866181629552146, -1.1002420139424332, -2.372545481333927)), ('C', (-2.870580911972212, -1.1060415834214632, 2.3675906194329763)), ('C', (-1.7937960684250869, 2.8544974529184404, 0.0034619782618387413)), ('H', (-4.873624464543119, -0.5850554289894927, -2.389365933576627)), ('H', (-1.9909209594760948, -0.36018935841676386, -4.096007833095599)), ('H', (-2.7384531510063264, -3.166933432144885, -2.4268259745615777)), ('H', (-2.743417461537904, -3.172903076975206, 2.4168538897975367)), ('H', (-1.998343803696893, -0.3706093082725383, 4.094586759049254)), ('H', (-4.877917922300159, -0.5904487073515493, 2.3823096962241674)), ('H', (-0.9125204000917033, 3.6852210576697018, 1.6830165434984983)), ('H', (-0.9108857869931823, 3.6894011318592144, -1.6731370553146045)), ('H', (-3.7908794247963002, 3.400690664137893, 0.0030254515268579828))),
        'graph': ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None), 3: ('H', 0, None), 4: ('C', 0, None), 5: ('H', 0, None), 6: ('C', 0, None), 7: ('H', 0, None), 8: ('C', 0, None), 9: ('C', 0, None), 10: ('C', 0, None), 11: ('H', 0, None), 12: ('H', 0, None), 13: ('H', 0, None), 14: ('H', 0, None), 15: ('H', 0, None), 16: ('H', 0, None), 17: ('H', 0, None), 18: ('H', 0, None), 19: ('H', 0, None)}, {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None), frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, False), frozenset({1, 5}): (1, None), frozenset({4, 6}): (1, None), frozenset({4, 7}): (1, None), frozenset({8, 6}): (1, None), frozenset({9, 6}): (1, None), frozenset({10, 6}): (1, None), frozenset({8, 11}): (1, None), frozenset({8, 12}): (1, None), frozenset({8, 13}): (1, None), frozenset({9, 14}): (1, None), frozenset({9, 15}): (1, None), frozenset({16, 9}): (1, None), frozenset({17, 10}): (1, None), frozenset({10, 18}): (1, None), frozenset({10, 19}): (1, None)}),
        'smiles': '[CH2]C=CC(C)(C)C',
        'formula': {'C': 7, 'H': 13},
    },
    'InChI=1S/C8H15/c1-7(2)6-8(3,4)5/h6H,1H2,2-5H3/b7-6-': {
        'inchi': 'InChI=1S/C8H15/c1-7(2)6-8(3,4)5/h6H,1H2,2-5H3/b7-6-',
        'geom': (('C', (-5.013415064947736, -2.0050050882891317, -0.03441380247071251)), ('C', (-2.956257315048219, -0.38735417147022017, -0.028551872029542322)), ('H', (-6.938091675095906, -1.2975539878798614, 0.018453175615095692)), ('H', (-4.765566145237189, -4.041653640549056, -0.09874385923354789)), ('C', (-3.4718445219397567, 2.429972368560465, 0.011336467026621512)), ('C', (-0.5176923617952974, -1.408040605257004, -0.0855025482724649)), ('H', (-5.442137231030354, 2.7997917713125617, 0.5130946581320586)), ('H', (-3.1290104082591976, 3.2783384564190006, -1.8447449744935553)), ('H', (-2.281311393722948, 3.4232917885983687, 1.37801663657573)), ('C', (2.042058838157109, -0.12047948912856385, -0.0067198661021280375)), ('H', (-0.43926494813654665, -3.4667384819487697, -0.08574821266877442)), ('C', (4.051410296534795, -2.1075000538817137, -0.6961694354402874)), ('C', (2.276604085944433, 2.082793774517478, -1.8945146914597384)), ('C', (2.590754046766668, 0.8422055806895431, 2.691075827314973)), ('H', (5.9521119513292815, -1.290420271756258, -0.6051583355121131)), ('H', (3.9868043397576423, -3.7156248613021927, 0.6081951253957237)), ('H', (3.7522590919703185, -2.827077526838047, -2.615186315842712)), ('H', (4.24332710238404, 2.7304048074380014, -1.974402863413468)), ('H', (1.116969869604986, 3.704794840877194, -1.3594500773931069)), ('H', (1.7191764529091338, 1.4888774197257137, -3.7976389751470627)), ('H', (4.4749376054985275, 1.7028359911022999, 2.7808945100579834)), ('H', (2.5270872838738683, -0.7134812164756067, 4.054832590547998)), ('H', (1.2129887437656235, 2.2646515682009123, 3.288346445979425))),
        'graph': ({0: ('C', 0, None), 1: ('C', 0, None), 2: ('H', 0, None), 3: ('H', 0, None), 4: ('C', 0, None), 5: ('C', 0, None), 6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None), 9: ('C', 0, None), 10: ('H', 0, None), 11: ('C', 0, None), 12: ('C', 0, None), 13: ('C', 0, None), 14: ('H', 0, None), 15: ('H', 0, None), 16: ('H', 0, None), 17: ('H', 0, None), 18: ('H', 0, None), 19: ('H', 0, None), 20: ('H', 0, None), 21: ('H', 0, None), 22: ('H', 0, None)}, {frozenset({0, 1}): (1, None), frozenset({0, 2}): (1, None), frozenset({0, 3}): (1, None), frozenset({1, 4}): (1, None), frozenset({1, 5}): (1, True), frozenset({4, 6}): (1, None), frozenset({4, 7}): (1, None), frozenset({8, 4}): (1, None), frozenset({9, 5}): (1, None), frozenset({10, 5}): (1, None), frozenset({9, 11}): (1, None), frozenset({9, 12}): (1, None), frozenset({9, 13}): (1, None), frozenset({11, 14}): (1, None), frozenset({11, 15}): (1, None), frozenset({16, 11}): (1, None), frozenset({17, 12}): (1, None), frozenset({18, 12}): (1, None), frozenset({19, 12}): (1, None), frozenset({20, 13}): (1, None), frozenset({21, 13}): (1, None), frozenset({13, 22}): (1, None)}),
        'smiles': '[CH2]C(C)=CC(C)(C)C',
        'formula': {'C': 8, 'H': 15},
    },
}


# # "constructor"
def from_data(fml_slyr, main_lyr_dct=None,
              char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None):
    """ Build an InChI string from each of the various layers.

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

    main_dct = dict_.empty_if_none(main_lyr_dct)
    char_dct = dict_.empty_if_none(char_lyr_dct)
    ste_dct = dict_.empty_if_none(ste_lyr_dct)
    iso_dct = dict_.empty_if_none(iso_lyr_dct)

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

    return ich


# # recalculate/standardize
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
        if any(hardcoded_object_from_inchi_by_key('inchi', ich)
               for ich in ichs):
            ref_ichs = []
            for ich_i in ichs:
                ref_ichs.append(recalculate(ich_i))
            ref_ichs.sort()
            ret = join(ref_ichs)
            return ret
        # raise error.FailedInchiGenerationError

    ret = hardcoded_object_from_inchi_by_key('inchi', ich)
    if ret is None:
        _options = '-SUU' if stereo else ''
        rdm = rdkit_.from_inchi(ich)
        ret = rdkit_.to_inchi(rdm, options=_options, with_aux_info=False)

    return ret


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
        iso_dct = dict_.by_key(isotope_sublayers(ich), ISO_NONSTE_PFXS)
    ich = from_data(fml_slyr,
                    main_lyr_dct=main_dct,
                    char_lyr_dct=char_dct,
                    ste_lyr_dct=ste_dct,
                    iso_lyr_dct=iso_dct)
    return recalculate(ich)


# # getters
def version(ich):
    """ Determine version of InChI the string corresponds to.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    ptt = app.capturing(version_pattern())
    ver = apf.first_capture(ptt, ich)
    return ver


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


def formula_string(ich):
    """ Generate a formula string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return formula_sublayer(ich)


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


# # conversions
def inchi_key(ich):
    """ Generate an InChIKey from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return rdkit_.inchi_to_inchi_key(ich)


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

    smi = hardcoded_object_from_inchi_by_key('smiles', ich)
    if smi is None:
        ich = standard_form(ich)
        rdm = rdkit_.from_inchi(ich)
        smi = rdkit_.to_smiles(rdm)

    return smi


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

    fml = hardcoded_object_from_inchi_by_key('formula', ich)
    if fml is None:
        ich = standard_form(ich)
        rdm = rdkit_.from_inchi(ich)
        fml = rdkit_.to_formula(rdm)

    return fml


# # properties
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


def has_multiple_components(ich):
    """ Determine if the InChI string has multiple components.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return len(split(ich)) > 1


def is_chiral(ich):
    """ Determine if the InChI string has chirality information.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return ste_dct['s'] == '1' or iso_dct['s'] == '1'


def has_stereo(ich):
    """ Determine if the InChI string has stereochemistry information.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)
    return bool(ste_dct or
                any(pfx in iso_dct for pfx in STE_PFXS))


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


# # comparisons
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


# # sort
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


# # split/join
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


# # hardcoded inchi workarounds
def hardcoded_object_from_inchi_by_key(key, ich):
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


def hardcoded_object_to_inchi_by_key(key, obj, comp=operator.eq):
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


# # helpers
def version_pattern():
    """ Build the autoparse regex pattern for the InChI string version.

        :rtype: str
    """
    ptt = app.preceded_by('InChI=') + _sublayer_pattern()
    return ptt


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
        return ({},) * count
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
        dct = {}
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


def _sublayer_pattern(key_ptt='',
                      val_ptt=NONSLASHES):
    """ Build the autoparse regex pattern for an arbitrary InChI sublayer.

        :rtype: str
    """
    return key_ptt + val_ptt
