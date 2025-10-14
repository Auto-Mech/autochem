""" Level 3 InChI functions (depend on extern and L1-2)
"""

import numpy

from ...amchi import base as amchi_base
from ...extern import rdkit_
from ...util import dict_

MAIN_PFXS = ("c", "h")
CHAR_PFXS = ("q", "p")
STE_PFXS = ("b", "t", "m", "s")
ISO_NONSTE_PFXS = ("i", "h")
ISO_PFXS = ISO_NONSTE_PFXS + STE_PFXS


# # "constructor"
def from_data(
    fml_lyr, main_lyr_dct=None, char_lyr_dct=None, ste_lyr_dct=None, iso_lyr_dct=None
):
    """Build an InChI string from layers

    :param fml_lyr: The formula layer
    :type fml_lyr: str
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
    return amchi_base.from_data(
        fml_lyr=fml_lyr,
        main_lyr_dct=main_lyr_dct,
        char_lyr_dct=char_lyr_dct,
        ste_lyr_dct=ste_lyr_dct,
        iso_lyr_dct=iso_lyr_dct,
        pfx="InChI",
    )


# # recalculate/standardize
def recalculate(ich, stereo=False, racem=False):
    """Recalculate an InChI string.

    :param ich: InChI string
    :type ich: str
    :param stereo: force stereochem in recalculated InChI
    :type stereo: bool
    :rtype: str
    """
    _options = "-SUU" if stereo else ""
    _options += "-SRac" if racem else ""
    rdm = rdkit_.from_inchi(ich)
    if rdm is not None:
        ret = rdkit_.to_inchi(rdm, options=_options, with_aux_info=False)
    else:
        raise ValueError(f"Invalid InChI: {ich}")

    return ret


def standard_form(ich, stereo=True, racem=False, ste_dct=None, iso_dct=None):
    """Return an InChI string in standard form.

    :param ich: InChI string
    :type ich: str
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param racem: parameter to designate a racemic mixture, if chiral
    :type racem: bool
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

    extra_ste_dct = ste_dct
    extra_iso_dct = iso_dct

    if stereo:
        ste_dct = stereo_layers(ich)
        iso_dct = isotope_layers(ich)
    else:
        ste_dct = {}
        iso_dct = dict_.by_key(isotope_layers(ich), ISO_NONSTE_PFXS)

    if extra_ste_dct is not None:
        ste_dct.update(extra_ste_dct)

    if extra_iso_dct is not None:
        iso_dct.update(extra_iso_dct)

    ich = amchi_base.standard_form(
        ich,
        stereo=stereo,
        racem=racem,
        ste_dct=ste_dct,
        iso_dct=iso_dct,
        racem_m_layer=False,
    )

    ich = recalculate(ich, racem=(racem and is_enantiomer(ich)))

    if ich is not None:
        recalc_ste_dct = stereo_layers(ich)
        if "s" in recalc_ste_dct:
            recalc_ste_dct.pop("s")
        if "s" in ste_dct:
            ste_dct.pop("s")

        recalc_iso_dct = isotope_layers(ich)
        if "s" in recalc_iso_dct:
            recalc_iso_dct.pop("s")
        if "s" in iso_dct:
            iso_dct.pop("s")

        # If we were attempting to force special stereo and it failed, return
        # None
        if extra_ste_dct is not None and recalc_ste_dct != ste_dct:
            ich = None

        if extra_iso_dct is not None and recalc_iso_dct != iso_dct:
            ich = None

    return ich


# # getters
def version(ich):
    """Determine version of InChI the string corresponds to.

    :param ich: InChI string
    :type ich: str
    :rtype: str
    """
    return amchi_base.version(ich)


def formula_layer(ich):
    """Parse the InChI string for the formula sublayer.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: str]
    """
    return amchi_base.formula_layer(ich)


def main_layers(ich):
    """Parse the InChI string for the sublayers of the connectivity layer,
    organized by prefix.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: str]
    """
    return amchi_base.main_layers(ich)


def charge_layers(ich):
    """Parse the InChI string for the sublayers of the charge layer,
    organized by prefix.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: str]
    """
    return amchi_base.charge_layers(ich)


def stereo_layers(ich):
    """Parse the InChI string for the sublayers of the stereochemisty layer,
    organized by prefix.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: str]
    """
    return amchi_base.stereo_layers(ich)


def isotope_layers(ich):
    """Parse the InChI string for the sublayers of the isotope layer,
    organized by prefix.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: str]
    """
    return amchi_base.isotope_layers(ich)


# # setters
def reflect(ich):
    """If this is an enantiomer, flip to the other enantiomer by changing the
    m-layer

    :param ich: InChI string
    :type ich: str
    :returns: the other enantiomer
    :rtype: bool
    """
    ich = amchi_base.reflect(ich)
    return recalculate(ich)


def stereo_bonds(ich, iso=True, one_indexed=False):
    """Parse the stereo bonds from the stereochemistry layer.

    :param ich: InChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    """
    bnd_ste_dct = amchi_base.bond_stereo_parities(ich, one_indexed=one_indexed)
    if iso:
        bnd_ste_dct.update(
            amchi_base.bond_isotope_stereo_parities(ich, one_indexed=one_indexed)
        )
    bnds = tuple(sorted(tuple(sorted(k, reverse=True)) for k in bnd_ste_dct))
    return bnds


def unassigned_stereo_bonds(ich, iso=True, one_indexed=False):
    """Parse the stereo bonds wth missing assignments from the stereochemistry
    layer.

        :param ich: InChI string
        :type ich: str
        :param iso: Include isotope stereochemistry?
        :type iso: bool
        :param one_indexed: Return indices in one-indexing?
        :type one_indexed: bool
    """
    bnd_ste_dct = amchi_base.bond_stereo_parities(ich, one_indexed=one_indexed)
    if iso:
        bnd_ste_dct.update(
            amchi_base.bond_isotope_stereo_parities(ich, one_indexed=one_indexed)
        )
    bnd_ste_dct = dict_.filter_by_value(bnd_ste_dct, lambda x: x is None)
    bnds = tuple(sorted(tuple(sorted(k, reverse=True)) for k in bnd_ste_dct))
    return bnds


def stereo_atoms(ich, iso=True, one_indexed=False):
    """Parse the stereo atoms from the stereochemistry layer.

    :param ich: InChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :param one_indexed: Return indices in one-indexing?
    :type one_indexed: bool
    """
    atm_ste_dct = amchi_base.atom_stereo_parities(ich, one_indexed=one_indexed)
    if iso:
        atm_ste_dct.update(
            amchi_base.atom_isotope_stereo_parities(ich, one_indexed=one_indexed)
        )
    return tuple(sorted(atm_ste_dct))


def is_enantiomer(ich, iso=True):
    """Is this InChI an enantiomer? (I.e., is it chiral?)

    Determined based on whether or not the InChI has an s-layer.

    :param ich: InChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :returns: whether or not the InChI is an enantiomer
    :rtype: bool
    """
    return amchi_base.is_enantiomer(ich, iso=iso)


def are_enantiomers(ich_a, ich_b):
    """Are these InChI enantiomers of eachother?

    :param ich: InChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :returns: whether or not the InChI is enantiomeric
    :rtype: bool
    """
    return amchi_base.are_enantiomers(ich_a, ich_b)


def are_diastereomers(ich_a, ich_b):
    """Are these InChI diastereomers of each other?

    Checks if main layer is the same, if so then checks
    if the stereo layers differ in any way.

    :param ich: InChI string
    :type ich: str
    :param iso: Include isotope stereochemistry?
    :type iso: bool
    :returns: whether or not the InChI is enantiomeric
    :rtype: bool
    """
    return amchi_base.are_diastereomers(ich_a, ich_b)


# # conversions
def inchi_key(ich):
    """Generate an InChIKey from an InChI string.

    :param ich: InChI string
    :type ich: str
    :rtype: str
    """
    return rdkit_.inchi_to_inchi_key(ich)


def smiles(ich):
    """Convert a SMILES string into an InChI string.

    :param smi: SMILES string
    :type smi: str
    :rtype: str
    """
    rdm = rdkit_.from_inchi(ich)
    smi = rdkit_.to_smiles(rdm)
    return smi


def formula(ich):
    """Generate a formula dictionary from an InChI string.

    :param ich: InChI string
    :type ich: str
    :rtype: dict[str: int]
    """
    return amchi_base.formula(ich)


def formula_string(ich):
    """Generate a formula string from an InChI string.

    :param ich: InChI string
    :type ich: str
    :rtype: str
    """
    return amchi_base.formula_string(ich)


def without_stereo(ich):
    """Remove all stereo layers"""

    return standard_form(ich, stereo=False)


def racemic(ich):
    """If chiral, convert the InChI into a racemic mixture

    This drops the /m layer and replaces /s1 with /s3, indicating a racemic
    mixture. The chirality of the species is still implied by the presence
    of the /s layer.

    :param ich: InChI string
    :type ich: str
    """
    return standard_form(ich, racem=True)


def connectivity(ich, parse_connection_layer=True, parse_h_layer=True):
    """Return the 'c' and 'h' layers of the connectivity string

    The user may also specify what combination of the two layers
    that they wish to return
    """
    return amchi_base.connectivity(
        ich, parse_connection_layer=parse_connection_layer, parse_h_layer=parse_h_layer
    )


# # properties
def is_standard_form(ich):
    """Determine if the InChI string is closed.

    :param ich: InChI string
    :type ich: str
    :rtype: bool
    """
    return ich == standard_form(ich)


def has_multiple_components(ich):
    """Determine if the InChI string has multiple components.

    :param ich: InChI string
    :type ich: str
    :rtype: bool
    """
    return amchi_base.has_multiple_components(ich)


def has_stereo(ich):
    """Determine if the InChI string has stereochemistry information.

    :param ich: InChI string
    :type ich: str
    :rtype: bool
    """
    return amchi_base.has_stereo(ich)


def low_spin_multiplicity(ich):
    """Guess spin multiplicity based on the number of electrons.

    :param ich: InChI string
    :type ich: str
    :rtype: int
    """
    return amchi_base.low_spin_multiplicity(ich)


# # comparisons
def same_connectivity(ich1, ich2):
    """Determine if two InChI strings have the same connectivity.

    :param ich1: InChI string 1
    :type ich1: str
    :param ich2: InChI string 2
    :type ich2: str
    :rtype: bool
    """
    return standard_form(ich1, stereo=False) == standard_form(ich2, stereo=False)


def equivalent(ich1, ich2):
    """Determine if two InChI strings are equivalent. Currently
    the strings are only checked up to the isotope sublayer.

    :param ich1: InChI string 1
    :type ich1: str
    :param ich2: InChI string 2
    :type ich2: str
    :rtype: bool
    """
    return amchi_base.equivalent(ich1, ich2)


# # split/join
def split(ich):
    """Split a multi-component InChI into InChIs for each of its components.

    (fix this for /s [which should be removed in split/join operations]
     and /m, which is joined as /m0110..  with no separators)

    :param ich: InChI string
    :type ich: str
    :rtype: tuple(str)
    """
    return amchi_base.split(ich)


def join(ichs):
    """Join separate InChI strings into one multi-component InChI string.

    Currently:
    (fix for /s [which should be removed in split/join operations] and /m,
     which is joined as /m0110..  with no separators).

    :param ichs: sequence of InChI strings
    :type ichs: tuple(str)
    :rtype: str
    """
    return amchi_base.join(ichs)


# # sort
def sorted_(ichs):
    """Sort a sequence of InChI strings in their standard form sort order.

    :param ichs: sequence of InChI strings
    :type ichs: tuple(str)
    :rtype: tuple(str)
    """
    return tuple(ichs[idx] for idx in argsort(ichs))


def argsort(ichs):
    """Determine the sort order for the InChI standard form.

    :param ichs: sequence of InChI strings
    :type ichs: tuple(str)
    """

    assert not any(map(has_multiple_components, ichs))
    ref_ichs = list(map(standard_form, split(recalculate(join(ichs)))))
    idxs = tuple(numpy.argsort(list(map(ref_ichs.index, ichs))))
    return idxs
