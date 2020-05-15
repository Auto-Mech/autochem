""" inchi "constructor"
"""
from automol import dict_


MAIN_PFXS = ('c', 'h')
CHAR_PFXS = ('q', 'p')
STE_PFXS = ('b', 't', 'm', 's')
ISO_NONSTE_PFXS = ('i', 'h')
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS


def from_data(formula_sublayer, main_sublayer_dct=None,
              charge_sublayer_dct=None, stereo_sublayer_dct=None,
              isotope_sublayer_dct=None):
    """ calculate an inchi string from layers
    """
    main_dct = dict_.empty_if_none(main_sublayer_dct)
    char_dct = dict_.empty_if_none(charge_sublayer_dct)
    ste_dct = dict_.empty_if_none(stereo_sublayer_dct)
    iso_dct = dict_.empty_if_none(isotope_sublayer_dct)

    fml_slyr = formula_sublayer
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
