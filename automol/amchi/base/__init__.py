""" Level 2 AMChI functions (depend on L1)
"""

from automol.util import dict_


MAIN_PFXS = ('c', 'h')
CHAR_PFXS = ('q', 'p')
STE_PFXS = ('b', 't', 'm', 's')
ISO_NONSTE_PFXS = ('i', 'h')
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS


# # constructor
def from_data(fml_lyr, main_lyrs_dct=None,
              char_lyrs_dct=None, ste_lyrs_dct=None,
              iso_lyrs_dct=None):
    """ Build an InChI string from each of the various layers.

        :param fml_lyr: sublayer of InChI string containing molecular formula
        :type fml_lyr: str
        :param main_lyrs_dct: information for connectivity layer of InChI
        :type main_lyrs_dct: dict[str: str]
        :param char_lyrs_dct: information for charge layer of InChI
        :type char_lyrs_dct: dict[str: str]
        :param ste_lyrs_dct: information for stereochemistry layer of InChI
        :type ste_lyrs_dct: dict[str: str]
        :param iso_lyrs_dct: information for isotope layer of InChI
        :type iso_lyrs_dct: dict[str: str]
        :rtype: str
    """

    main_dct = dict_.empty_if_none(main_lyrs_dct)
    char_dct = dict_.empty_if_none(char_lyrs_dct)
    ste_dct = dict_.empty_if_none(ste_lyrs_dct)
    iso_dct = dict_.empty_if_none(iso_lyrs_dct)

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

    ich = '/'.join(['AMChI=1', fml_lyr] + main_lyrs + char_lyrs +
                   ste_lyrs + iso_lyrs)

    return ich
