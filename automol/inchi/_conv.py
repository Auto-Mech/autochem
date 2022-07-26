""" Level 4 functions depending on other basic types (graph, geom)
"""

import automol.formula
import automol.geom
import automol.graph
import automol.amchi
from automol.inchi.base import standard_form
from automol.inchi.base import has_stereo
from automol.inchi.base import equivalent
from automol.inchi.base import hardcoded_object_from_inchi_by_key


# # conversions
def graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """
    gra = hardcoded_object_from_inchi_by_key('graph', ich)
    if gra is None:
        gra = automol.amchi.graph(ich, stereo=stereo)
    return gra


def geometry(ich, check=True):
    """ Generate a molecular geometry from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry data structure
    """
    geo = hardcoded_object_from_inchi_by_key('geom', ich)
    if geo is None:
        geo = automol.amchi.geometry(ich, check=check)
    return geo


def conformers(ich, nconfs=1):
    """ Generate a list of molecular geometries for various conformers
        of a species from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param nconfs: number of conformer structures to generate
        :type: int
        :rtype: automol molecular geometry data structure
    """
    geo = hardcoded_object_from_inchi_by_key('geom', ich)
    if geo is None:
        geos = automol.amchi.conformers(ich, nconfs=nconfs)
    else:
        geos = [geo] * nconfs
    return geos


def zmatrix(ich, check=True):
    """ Generate a z-matrix from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol z-matrix data structure
    """
    zma = automol.amchi.zmatrix(ich, check=check)
    return zma


def amchi(ich, stereo=True):
    """ Convert an InChI to an AMChI string

        Only for good InChIs, where this can be done validly.

        :param ich: InChI string
        :type ich: str
        :returns: AMChI string
        :rtype: str
    """
    gra = graph(ich, stereo=stereo)

    if stereo:
        assert not is_bad(ich, gra=gra), (
            "Don't use this function with bad InChIs if stereo=True")

    ach = automol.graph.amchi(gra)
    return ach


# # derived properties
def is_complete(ich):
    """ Determine if the InChI string is complete
        (has all stereo-centers assigned).

        :param ich: InChI string
        :type ich: str
        :rtype: bool
    """

    gra = graph(ich, stereo=False)
    ste_atm_keys = automol.graph.stereogenic_atom_keys(gra)
    ste_bnd_keys = automol.graph.stereogenic_bond_keys(gra)
    graph_has_stereo = bool(ste_atm_keys or ste_bnd_keys)

    _complete = equivalent(ich, standard_form(ich)) and not (
        has_stereo(ich) ^ graph_has_stereo)

    return _complete


def is_bad(ich, gra=None):
    """ Determine if the InChI string is bad, i.e. one of the InChI failure
        cases (resonance bond stereo, vinyl bond stereo, etc.

        :param ich: InChI string
        :type ich: str
        :param gra: A graph version of the InChI, to avoid recalculating
        :type gra: automol graph
        :returns: True if it is a bad InChI
        :rtype: bool
    """
    gra = graph(ich) if gra is None else gra
    ret = automol.graph.base.inchi_is_bad(gra, ich)
    return ret


# # derived transformations
def add_stereo(ich):
    """ Add stereochemistry to an InChI string converting to/from geometry.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    geo = geometry(ich)
    ich = automol.geom.inchi(geo, stereo=True)
    return ich


def expand_stereo(ich):
    """ Obtain all possible stereoisomers of an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: list[str]
    """
    gra = graph(ich, stereo=False)
    sgrs = automol.graph.expand_stereo(gra)
    ste_ichs = [automol.graph.inchi(sgr, stereo=True) for sgr in sgrs]
    return ste_ichs
