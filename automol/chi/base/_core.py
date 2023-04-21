""" Level 3 ChI functions (depend on inchi.base and extern and L1-2)
"""
import automol.inchi.base
import automol.amchi.base


# # "constructor"
def from_data(fml_slyr, main_lyr_dct=None,
              char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None, amchi=False):
    """ Build a ChI string from each of the various layers.

        :param fml_slyr: sublayer of ChI string containing molecular formula
        :type fml_slyr: str
        :param main_lyr_dct: information for connectivity layer of ChI
        :type main_lyr_dct: dict[str: str]
        :param char_lyr_dct: information for charge layer of ChI
        :type char_lyr_dct: dict[str: str]
        :param ste_lyr_dct: information for stereochemistry layer of ChI
        :type ste_lyr_dct: dict[str: str]
        :param iso_lyr_dct: information for isotope layer of ChI
        :type iso_lyr_dct: dict[str: str]
        :param amchi: create an AMChI instead of an InChI?
        :type amchi: bool
        :returns: the ChI string
        :rtype: str
    """
    fun_ = (automol.amchi.base.from_data if amchi else
            automol.inchi.base.from_data)
    return fun_(fml_slyr, main_lyr_dct=main_lyr_dct,
                char_lyr_dct=char_lyr_dct, ste_lyr_dct=ste_lyr_dct,
                iso_lyr_dct=iso_lyr_dct)


# # recalculate/standardize
def recalculate(chi, stereo=False):
    """ Recalculate a ChI string.

        Does nothing for AMChI.

        :param chi: ChI string
        :type chi: str
        :param stereo: force stereochem in recalculated ChI
        :type stereo: bool
        :rtype: str
    """
    pfx = automol.amchi.base.prefix(chi)
    if pfx == 'AMChI':
        ret = chi
    elif pfx == 'InChI':
        ret = automol.inchi.base.recalculate(chi, stereo=stereo)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


def standard_form(chi, stereo=True, ste_dct=None, iso_dct=None):
    """ Return a ChI string in standard form.

        :param chi: ChI string
        :type chi: str
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
    pfx = prefix(chi)
    if pfx == 'AMChI':
        ret = automol.amchi.base.standard_form(
            chi, stereo=stereo, ste_dct=ste_dct, iso_dct=iso_dct)
    elif pfx == 'InChI':
        ret = automol.inchi.base.standard_form(
            chi, stereo=stereo, ste_dct=ste_dct, iso_dct=iso_dct)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


# # getters
def prefix(chi):
    """ Determine the chemical identifier prefix (InChI or AMChI).

        :param chi: ChI string
        :type chi: str
        :rtype: str
    """
    return automol.amchi.base.prefix(chi)


def are_enantiomers(chi_a, chi_b, log=False):
    """ Assess if ChI string for two species are enantiomers of one another.

        :param chi: ChI string
        :type chi: str
    """
    pfx_a, pfx_b = prefix(chi_a), prefix(chi_b)
    if pfx_a == pfx_b:
        if all(pfx == 'AMChI' for pfx in (pfx_a, pfx_b)):
            ret = automol.amchi.are_enantiomers(chi_a, chi_b)
        elif all(pfx == 'InChI' for pfx in (pfx_a, pfx_b)):
            ret = automol.inchi.are_enantiomers(chi_a, chi_b)
        else:
            if log:
                print("ChI string '{chi_a}' or '{chi_b}' has unknown prefix")
    else:
        if log:
            print("ChI string prefixes of '{chi_a}' and '{chi_b}' don't match")
        ret = False

    return ret


def are_diastereomers(chi_a, chi_b, log=False):
    """ Assess if ChI string for two species are diastereomers of one another.

        :param chi: ChI string
        :type chi: str
    """
    pfx_a, pfx_b = prefix(chi_a), prefix(chi_b)
    if pfx_a == pfx_b:
        if all(pfx == 'AMChI' for pfx in (pfx_a, pfx_b)):
            ret = automol.amchi.are_diastereomers(chi_a, chi_b)
        elif all(pfx == 'InChI' for pfx in (pfx_a, pfx_b)):
            ret = automol.inchi.are_diastereomers(chi_a, chi_b)
        else:
            if log:
                print("ChI string '{chi_a}' or '{chi_b}' has unknown prefix")
    else:
        if log:
            print("ChI string prefixes of '{chi_a}' and '{chi_b}' don't match")
        ret = False

    return ret


def reflect(chi):
    """ If this is an enantiomer, flip to the other enantiomer by changing the
        m-layer

        :param chi: ChI string
        :type chi: str
        :returns: the other enantiomer
        :rtype: bool
    """
    pfx = prefix(chi)
    if pfx == 'AMChI':
        ret = automol.amchi.base.reflect(chi)
    elif pfx == 'InChI':
        ret = automol.inchi.base.reflect(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


def racemic(chi):
    """ If chiral, convert the ChI into a racemic mixture

        This drops the /m layer and replaces /s1 with /s3, indicating a racemic
        mixture. The chirality of the species is still implied by the presence
        of the /s layer.

        :param chi: ChI string
        :type chi: str
    """
    pfx = prefix(chi)
    if pfx == 'AMChI':
        ret = automol.amchi.base.racemic(chi)
    elif pfx == 'InChI':
        ret = automol.inchi.base.racemic(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


# # properties
def is_standard_form(chi):
    """ Determine whether the ChI string is in standard form.

        :param chi: ChI string
        :type chi: str
        :rtype: bool
    """
    return chi == standard_form(chi)


# # split/join
def split(chi):
    """ Split a multi-component ChI into ChIs for each of its components.

        :param chi: ChI string
        :type chi: str
        :returns: the split ChI strings
        :rtype: tuple[str]
    """
    pfx = prefix(chi)
    if pfx == 'AMChI':
        ret = automol.amchi.base.split(chi)
    elif pfx == 'InChI':
        ret = automol.inchi.base.split(chi)
    else:
        raise ValueError(f"ChI string '{chi}' has unknown prefix '{pfx}'.")
    return ret


# # sort
def sorted_(chis):
    """ Sort a sequence of ChI strings in the standard sort order (see argsort)

        :param chis: sequence of ChI strings
        :type chis: tuple(str)
        :rtype: tuple(str)
    """
    pfxs = list(map(prefix, chis))
    if 'AMChI' in pfxs:
        ret = automol.amchi.base.sorted_(chis)
    elif set(pfxs) == {'InChI'}:
        ret = automol.inchi.base.sorted_(chis)
    else:
        raise ValueError(f"One of these has an unknown prefix: {chis}.")
    return ret


def argsort(chis):
    """ Determine the standard sort order for multiple ChIs.

        Follows the sort order for multicomponent InChIs as much as possible.

        :param chis: sequence of ChI strings
        :type chis: tuple(str)
    """
    pfxs = list(map(prefix, chis))
    if 'AMChI' in pfxs:
        ret = automol.amchi.base.argsort(chis)
    elif set(pfxs) == {'InChI'}:
        ret = automol.inchi.base.argsort(chis)
    else:
        raise ValueError(f"One of these has an unknown prefix: {chis}.")
    return ret
