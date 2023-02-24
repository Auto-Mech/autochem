""" Level 3 InChI functions (depend on extern and L1-2)
"""

import itertools
import warnings
import numpy
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast
import automol.formula
from automol.util import dict_
from automol.extern import rdkit_
import automol.amchi.base


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
        :param stereo: force stereochem in recalculated InChI
        :type stereo: bool
        :rtype: str
    """
    _options = '-SUU' if stereo else ''
    rdm = rdkit_.from_inchi(ich)
    if rdm is not None:
        ret = rdkit_.to_inchi(rdm, options=_options, with_aux_info=False)

    return ret


def standard_form(ich, stereo=True, ste_dct=None, iso_dct=None):
    """ Return an InChI string in standard form.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :param ste_dct: a dictionary to overwrite stereo information; layers
            not overwritten will be left in tact; if the attempted overwrite
            fails, the function will return None
        :type ste_dct: dict
        :param iso_dct: a dictionary to overwrite isotope stereo information;
            layers not overwritten will be left in tact; if the attempted
            overwrite fails, the function will return None
        :type iso_dct: dict
        :rtype: str
    """
    fml_slyr = formula_sublayer(ich)
    main_dct = main_sublayers(ich)
    char_dct = charge_sublayers(ich)

    extra_ste_dct = ste_dct
    extra_iso_dct = iso_dct

    if stereo:
        ste_dct = stereo_sublayers(ich)
        iso_dct = isotope_sublayers(ich)
    else:
        ste_dct = {}
        iso_dct = dict_.by_key(isotope_sublayers(ich), ISO_NONSTE_PFXS)

    if extra_ste_dct is not None:
        ste_dct.update(extra_ste_dct)

    if extra_iso_dct is not None:
        iso_dct.update(extra_iso_dct)

    ich = from_data(fml_slyr,
                    main_lyr_dct=main_dct,
                    char_lyr_dct=char_dct,
                    ste_lyr_dct=ste_dct,
                    iso_lyr_dct=iso_dct)
    ich = recalculate(ich)

    if ich is not None:
        recalc_ste_dct = stereo_sublayers(ich)
        if 's' in recalc_ste_dct:
            recalc_ste_dct.pop('s')
        if 's' in ste_dct:
            ste_dct.pop('s')

        recalc_iso_dct = isotope_sublayers(ich)
        if 's' in recalc_iso_dct:
            recalc_iso_dct.pop('s')
        if 's' in iso_dct:
            iso_dct.pop('s')

        # If we were attempting to force special stereo and it failed, return
        # None
        if extra_ste_dct is not None and recalc_ste_dct != ste_dct:
            ich = None

        if extra_iso_dct is not None and recalc_iso_dct != iso_dct:
            ich = None

    return ich


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
    return automol.amchi.base.formula_string(ich)


def formula_string(ich):
    """ Generate a formula string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return automol.amchi.base.formula_string(ich)


def main_sublayers(ich):
    """ Parse the InChI string for the sublayers of the connectivity layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.amchi.base.main_layers(ich)


def charge_sublayers(ich):
    """ Parse the InChI string for the sublayers of the charge layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.amchi.base.charge_layers(ich)


def stereo_sublayers(ich):
    """ Parse the InChI string for the sublayers of the stereochemisty layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.amchi.base.stereo_layers(ich)


def isotope_sublayers(ich):
    """ Parse the InChI string for the sublayers of the isotope layer,
        organized by prefix.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.amchi.base.isotope_layers(ich)


def stereo_atoms(ich, iso=True, one_indexed=False):
    """ Parse the stereo atoms from the stereochemistry layer.

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    if len(split(ich)) > 1:
        raise NotImplementedError("Multicomponent InChIs not implemented."
                                  "Call inchi.split() first")

    atm_ptt = (app.capturing(app.UNSIGNED_INTEGER) +
               app.one_of_these(list(map(app.escape, '+-'))))

    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)

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


def stereo_bonds(ich, iso=True, one_indexed=False):
    """ Parse the stereo bonds from the stereochemistry layer.

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    if len(split(ich)) > 1:
        raise NotImplementedError("Multicomponent InChIs not implemented."
                                  "Call inchi.split() first")

    bnd_ptt = '-'.join([app.capturing(app.UNSIGNED_INTEGER)]*2)

    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)

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


def unassigned_stereo_bonds(ich, iso=True, one_indexed=False):
    """ Parse the stereo bonds wth missing assignments from the stereochemistry
    layer.

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    if len(split(ich)) > 1:
        raise NotImplementedError("Multicomponent InChIs not implemented."
                                  "Call inchi.split() first")

    bnd_ptt = ('-'.join([app.capturing(app.UNSIGNED_INTEGER)]*2) +
               app.escape('?'))

    ste_dct = stereo_sublayers(ich)
    iso_dct = isotope_sublayers(ich)

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


def is_enantiomer(ich, iso=True):
    """ Is this InChI an enantiomer?

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the InChI is enantiomeric
        :rtype: bool
    """
    ste_dct = stereo_sublayers(ich)
    ret = 'm' in ste_dct
    if iso:
        iso_dct = isotope_sublayers(ich)
        ret = ret or 'm' in iso_dct
    return ret


def are_enantiomers(ich_a, ich_b):
    """ Are these InChI enantiomers of eachother?

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the InChI is enantiomeric
        :rtype: bool
    """
    ste_dct_a = stereo_sublayers(ich_a)
    ste_dct_b = stereo_sublayers(ich_b)
    enant = False
    if main_sublayers(ich_a) == main_sublayers(ich_b):
        if (len(ste_dct_b.keys()) == len(ste_dct_a.keys())
                and 'm' in ste_dct_a.keys()):
            if ste_dct_a['m'] != ste_dct_b['m']:
                if 'b' in ste_dct_a.keys():
                    if ste_dct_a['b'] == ste_dct_b['b']:
                        if 't' in ste_dct_a.keys():
                            if ste_dct_a['t'] == ste_dct_b['t']:
                                enant = True
                        else:
                            enant = True
                elif 't' in ste_dct_a.keys():
                    if ste_dct_a['t'] == ste_dct_b['t']:
                        enant = True
                else:
                    enant = True
    return enant


def are_diastereomers(ich_a, ich_b):
    """ Are these InChI diastereomers of each other?

        Checks if main layer is the same, if so then checks
        if the stereo layers differ in any way.

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :returns: whether or not the InChI is enantiomeric
        :rtype: bool
    """

    diast = False
    if ich_a != ich_b:  # chk not same InChIs
        if main_sublayers(ich_a) == main_sublayers(ich_b):
            ste_dct_a = stereo_sublayers(ich_a)
            ste_dct_b = stereo_sublayers(ich_b)
            # b-lyr are diastereomers; t-lyr may be, need check
            if len(ste_dct_a.keys()) == len(ste_dct_b.keys()):
                if 'b' in ste_dct_a.keys():
                    if ste_dct_a['b'] != ste_dct_b['b']:
                        diast = True
                elif 't' in ste_dct_a.keys():
                    if ste_dct_a['t'] != ste_dct_b['t']:
                        diast = True

    return diast


def reflect(ich):
    """ If this is an enantiomer, flip to the other enantiomer by changing the
        m-layer

        :param ich: InChI string
        :type ich: str
        :returns: the other enantiomer
        :rtype: bool
    """
    ste_upd_dct = stereo_sublayers(ich)
    iso_upd_dct = isotope_sublayers(ich)

    refl_trans = str.maketrans('01', '10')
    if 'm' in ste_upd_dct:
        ste_upd_dct['m'] = ste_upd_dct['m'].translate(refl_trans)

    if 'm' in iso_upd_dct:
        iso_upd_dct['m'] = iso_upd_dct['m'].translate(refl_trans)

    ich = standard_form(ich, ste_dct=ste_upd_dct, iso_dct=iso_upd_dct)

    return ich


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
    rdm = rdkit_.from_inchi(ich)
    smi = rdkit_.to_smiles(rdm)
    return smi


def formula(ich):
    """ Generate a formula dictionary from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: int]
    """
    rdm = rdkit_.from_inchi(ich)
    fml = rdkit_.to_formula(rdm)
    return fml


def without_stereo(ich):
    """ Remove all stereo layers
    """

    return standard_form(ich, stereo=False)


def connectivity(ich, parse_connection_layer=True, parse_h_layer=True):
    """ Return the 'c' and 'h' layers of the connectivity string

        The user may also specify what combination of the two layers
        that they wish to return
    """

    # Read the two sublayers that are requested to be parsed
    conn_slyrs = main_sublayers(ich)

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


# # properties
def is_standard_form(ich):
    """ Determine if the InChI string is closed.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return ich == standard_form(ich)


def has_multiple_components(ich):
    """ Determine if the InChI string has multiple components.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return automol.amchi.base.has_multiple_components(ich)


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
    fml_dct1 = formula_sublayer(ich1)
    fml_dct2 = formula_sublayer(ich2)
    conn_dct1 = main_sublayers(ich1)
    conn_dct2 = main_sublayers(ich2)
    chg_dct1 = charge_sublayers(ich1)
    chg_dct2 = charge_sublayers(ich2)
    ste_dct1 = stereo_sublayers(ich1)
    ste_dct2 = stereo_sublayers(ich2)
    iso_dct1 = isotope_sublayers(ich1)
    iso_dct2 = isotope_sublayers(ich2)
    # Stereo layers get dropped upon split/joins, so remove these from the
    # equivalence test
    for dct in (ste_dct1, ste_dct2, iso_dct1, iso_dct2):
        if 's' in dct:
            dct.pop('s')
    return (fml_dct1 == fml_dct2 and conn_dct1 == conn_dct2 and
            chg_dct1 == chg_dct2 and ste_dct1 == ste_dct2 and
            iso_dct1 == iso_dct2)


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
    fml_slyrs = automol.amchi.base.split_layer_string(
        fml_slyr, count_sep_ptt='', sep_ptt=app.escape('.'))
    count = len(fml_slyrs)

    main_dcts, _ = automol.amchi.base.split_layers(main_dct, count)
    char_dcts, _ = automol.amchi.base.split_layers(char_dct, count)
    ste_dcts, warn_msg1 = automol.amchi.base.split_layers(ste_dct, count)
    iso_dcts, warn_msg2 = automol.amchi.base.split_layers(iso_dct, count)

    warn_msg = warn_msg1 if warn_msg1 else warn_msg2
    if warn_msg:
        warnings.warn(f"\n{ich}\n{warn_msg}")

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
    fml_slyr = automol.amchi.base.join_layer_strings(
        fmls, count_sep='', sep='.')
    main_dct, _ = automol.amchi.base.join_layers(
        list(map(main_sublayers, ichs)))
    char_dct, _ = automol.amchi.base.join_layers(
        list(map(charge_sublayers, ichs)))
    ste_dct, warn_msg1 = automol.amchi.base.join_layers(
        list(map(stereo_sublayers, ichs)))
    iso_dct, warn_msg2 = automol.amchi.base.join_layers(
        list(map(isotope_sublayers, ichs)))

    warn_msg = warn_msg1 if warn_msg1 else warn_msg2
    if warn_msg:
        warnings.warn(f"\n{ichs}\n{warn_msg}")

    return from_data(fml_slyr=fml_slyr,
                     main_lyr_dct=main_dct,
                     char_lyr_dct=char_dct,
                     ste_lyr_dct=ste_dct,
                     iso_lyr_dct=iso_dct)


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


# # helpers
def version_pattern():
    """ Build the autoparse regex pattern for the InChI string version.

        :rtype: str
    """
    return automol.amchi.base.version_pattern()
