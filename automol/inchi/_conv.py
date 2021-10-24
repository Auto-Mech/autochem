""" Level 4 functions depending on other basic types (graph, geom, zmat)
"""

import functools
from automol import error
import automol.formula
import automol.geom
import automol.graph
from automol.extern import rdkit_
from automol.extern import pybel_
from automol.inchi.base import standard_form
from automol.inchi.base import split
from automol.inchi.base import has_stereo
from automol.inchi.base import same_connectivity
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

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = split(ich)
    gras = [_inchi_connected_graph(ich, stereo=stereo) for ich in ichs]
    for idx, gra in enumerate(gras):
        if idx == 0:
            num = 0
        else:
            num = max(map(max, map(automol.graph.atom_keys, gras[:idx]))) + 1
        gras[idx] = automol.graph.transform_keys(gra, num.__add__)
    gra = functools.reduce(automol.graph.union, gras)
    return gra


def _inchi_connected_graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string where
        all all atoms are connected by at least one bond.

        :param ich: InChI string
        :type ich: str
        :param remove_stereo: parameter to include stereochemistry information
        :type remove_stereo: bool
        :rtype: automol molecular graph
    """

    gra = hardcoded_object_from_inchi_by_key('graph', ich)
    if gra is None:
        ich = standard_form(ich)
        if not stereo or not has_stereo(ich):
            rdm = rdkit_.from_inchi(ich)
            gra = rdkit_.to_connectivity_graph(rdm)
        else:
            geo = geometry(ich)
            gra = automol.geom.graph(geo, stereo=stereo)

    gra = automol.graph.implicit(gra)
    return gra


def geometry(ich, check=True):
    """ Generate a molecular geometry from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry data structure
    """

    # rdkit fails for multi-component inchis, so we split it up and space out
    # the geometries
    ichs = split(ich)
    geos = [_connected_geometry(ich, check=check) for ich in ichs]
    geos = [automol.geom.translate(geo, [50. * idx, 0., 0.])
            for idx, geo in enumerate(geos)]
    geo = functools.reduce(automol.geom.join, geos)
    return geo


def _connected_geometry(ich, check=True):
    """ Generate a molecular geometry from an InChI string where
        all atoms are connected by at least one bond.

        :param ich: InChI string
        :type ich: str
        :param check: check stereo and connectivity?
        :type check: bool
        :rtype: automol molecular geometry data structure
    """
    # print("inchi in:", ich)

    geo = hardcoded_object_from_inchi_by_key('geom', ich)
    if geo is None:

        def _gen1(ich):
            rdm = rdkit_.from_inchi(ich)
            geo, = rdkit_.to_conformers(rdm, nconfs=1)
            return geo

        def _gen2(ich):
            pbm = pybel_.from_inchi(ich)
            geo = pybel_.to_geometry(pbm)
            return geo

        def _gen3(ich):
            if has_stereo(ich):
                raise ValueError

            gra = graph(ich, stereo=False)
            gra = automol.graph.explicit(gra)
            geo = automol.graph.embed.geometry(gra)
            return geo

        for gen_ in [_gen1, _gen1, _gen1, _gen2, _gen3]:
            success = False
            try:
                geo = gen_(ich)
                geo_ich = automol.geom.inchi(geo)
                # Check connectivity
                same_conn = same_connectivity(ich, geo_ich)
                conn = automol.geom.connected(geo)
                _has_stereo = has_stereo(ich)
                ich_equiv = equivalent(ich, geo_ich)
                checks_pass = ((same_conn and conn) and
                               (not _has_stereo or ich_equiv))
                # print('original inchi', ich)
                # print('geometry inchi', geo_ich)
                if not check or checks_pass:
                    success = True
                    break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

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
        ich = standard_form(ich)

        def _gen1(ich):
            rdm = rdkit_.from_inchi(ich)
            geos = rdkit_.to_conformers(rdm, nconfs)
            return geos

        for gen_ in [_gen1]:
            success = False
            try:
                geos = gen_(ich)
                for geo in geos:
                    geo_ich = automol.geom.inchi(geo)
                    if same_connectivity(ich, geo_ich) and (
                            not has_stereo(ich) or
                            equivalent(ich, geo_ich)):
                        success = True  # fix
                        break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

    return geos


# # derived properties
def is_complete(ich):
    """ Determine if the InChI string is complete
        (has all stereo-centers assigned).

        Currently only checks species that does not have any
        resonance structures.

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
    """ Obtain all possible stereoisomers compatible with an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: list[str]
    """
    gra = graph(ich)
    sgrs = automol.graph.stereomers(gra)
    ste_ichs = [automol.graph.stereo_inchi(sgr) for sgr in sgrs]
    return ste_ichs
