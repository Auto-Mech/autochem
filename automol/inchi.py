""" InChI chemical identifiers
"""

import numpy
import autoparse.pattern as app
import autoparse.find as apf
import automol.convert.geom
import automol.convert.inchi


# "constructor"
def from_data(fml_slyr, main_lyr_dct=None, char_lyr_dct=None, ste_lyr_dct=None,
              iso_lyr_dct=None):
    """ Construct an InChI string from its constituent layers, where
        most layers are input as dictionary of prefixes for some part of the
        sublayer and the corresponding string for that sublayer part.

        :param fml_slyr: sublayer of InChI string containing molecular formula
    """
    return automol.convert.inchi.from_data(
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
    ptt = app.capturing(automol.convert.inchi.version_pattern())
    ver = apf.first_capture(ptt, ich)
    return ver


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
    return automol.convert.inchi.standard_form(ich, stereo=stereo)


def has_stereo(ich):
    """ Determine if the InChI string has stereochemistry information.

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """
    return automol.convert.inchi.has_stereo(ich)


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
    ste_dct = automol.convert.inchi.stereo_sublayers(ich)
    iso_dct = automol.convert.inchi.isotope_sublayers(ich)
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
    return automol.convert.inchi.equivalent(ich1, ich2)


def same_connectivity(ich1, ich2):
    """ Determine if two InChI strings have the same connectivity.

        :param ich1: InChI string 1
        :type ich1: str
        :param ich2: InChI string 2
        :type ich2: str
        :rtype: bool
    """
    return automol.convert.inchi.same_connectivity(ich1, ich2)


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
    return automol.convert.inchi.join(ichs)


def split(ich):
    """ Split a multi-component InChI into InChIs for each of its components.

        (fix this for /s [which should be removed in split/join operations]
         and /m, which is joined as /m0110..  with no separators)

        :param ich: InChI string
        :type ich: str
        :rtype: tuple(str)
    """
    return automol.convert.inchi.split(ich)


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
    return automol.convert.geom.geometry(ich)


def conformers(ich, nconfs=100):
    """ Generate a molecular geometry for many conformers from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param nconfs: number of conformers to generate
        :type nconfs: int
        :rtype: tuple(automol geometry data structure)
    """
    return automol.convert.geom.conformers(ich, nconfs)


def graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol graph data structure
    """
    return automol.convert.geom.inchi_graph(ich, stereo=stereo)


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


def formula_sublayer(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.convert.inchi.formula_sublayer(ich)


def main_sublayers(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.convert.inchi.main_sublayers(ich)


def charge_sublayers(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.convert.inchi.charge_sublayers(ich)


def isotope_sublayers(ich):
    """ Parse the InChI string for the formula sublayer.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: str]
    """
    return automol.convert.inchi.isotope_sublayers(ich)


def formula_string(ich):
    """ Generate a formula string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return automol.convert.inchi.formula_sublayer(ich)


def stereo_sublayers(ich):
    """ Generate a formula string from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return automol.convert.inchi.stereo_sublayers(ich)


def add_stereo(ich):
    """ Add stereochemistry to an InChI string converting to/from geometry.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    geo = automol.convert.geom.geometry(ich)
    ich = automol.convert.geom.inchi(geo, stereo=True)
    return ich


def expand_stereo(ich):
    """ Obtain all possible stereoisomers compatible with an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: list[str]
    """
    gra = automol.inchi.graph(ich)
    sgrs = automol.graph.stereomers(gra)
    ste_ichs = [automol.graph.stereo_inchi(sgr) for sgr in sgrs]
    return ste_ichs
