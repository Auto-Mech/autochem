""" inchi conversions
"""

import functools
import operator
from phydat import phycon
from automol import error
import automol.inchi
import automol.geom
import automol.graph
import automol.convert.geom
from automol.convert import _rdkit
from automol.convert import _pybel


def geometry(ich):
    """ Generate a molecular geometry from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: automol molecular geometry data structure
    """

    # rdkit fails for multi-component inchis, so we split it up and space out
    # the geometries
    ichs = automol.inchi.split(ich)
    geos = list(map(_connected_geometry, ichs))
    geos = [automol.geom.translate(geo, [50. * idx, 0., 0.])
            for idx, geo in enumerate(geos)]
    geo = functools.reduce(automol.geom.join, geos)
    return geo


def _connected_geometry(ich):
    """ Generate a molecular geometry from an InChI string where
        all atoms are connected by at least one bond.

        :param ich: InChI string
        :type ich: str
        :rtype: automol molecular geometry data structure
    """

    geo = object_from_hardcoded_inchi_by_key('geom', ich)
    if geo is None:
        ich = automol.inchi.standard_form(ich)

        def _gen1(ich):
            rdm = _rdkit.from_inchi(ich)
            geo, = _rdkit.to_conformers(rdm, nconfs=1)
            return geo

        def _gen2(ich):
            pbm = _pybel.from_inchi(ich)
            geo = _pybel.to_geometry(pbm)
            return geo

        def _gen3(ich):
            if automol.inchi.has_stereo(ich):
                raise ValueError

            gra = automol.convert.inchi.graph(ich, stereo=False)
            gra = automol.graph.explicit(gra)
            geo = automol.graph.embed.geometry(gra)
            return geo

        for gen_ in [_gen1, _gen1, _gen1, _gen2, _gen3]:
            success = False
            try:
                geo = gen_(ich)
                geo_ich = automol.convert.geom.inchi(geo)
                # Check connectivity
                same_conn = automol.inchi.same_connectivity(ich, geo_ich)
                conn = automol.geom.connected(geo)
                has_stereo = automol.inchi.has_stereo(ich)
                ich_equiv = automol.inchi.equivalent(ich, geo_ich)
                if (same_conn and conn) and (not has_stereo or ich_equiv):
                    success = True
                    break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

    return geo


def conformers(ich, nconfs):
    """ Generate a list of molecular geometries for various conformers
        of a species from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param nconfs: number of conformer structures to generate
        :type: int
        :rtype: automol molecular geometry data structure
    """

    geo = object_from_hardcoded_inchi_by_key('geom', ich)
    if geo is None:
        ich = automol.inchi.standard_form(ich)

        def _gen1(ich):
            rdm = _rdkit.from_inchi(ich)
            geos = _rdkit.to_conformers(rdm, nconfs)
            return geos

        # def _gen2(ich):
        #     pbm = _pybel.from_inchi(ich)
        #     geos = _pybel.to_conformers(pbm)
        #     return geos

        for gen_ in [_gen1]:
            success = False
            try:
                geos = gen_(ich)
                for geo in geos:
                    geo_ich = automol.convert.geom.inchi(geo)
                    if automol.inchi.same_connectivity(ich, geo_ich) and (
                            not automol.inchi.has_stereo(ich) or
                            automol.inchi.equivalent(ich, geo_ich)):
                        success = True  # fix
                        break
            except (RuntimeError, TypeError, ValueError):
                continue

        if not success:
            raise error.FailedGeometryGenerationError

    return geos


def recalculate(ich, stereo=False):
    """ Recalculate an InChI string

        :param ich: InChI string
        :type ich: str
        :param stereo: force the same stereochem in recalculated InChI
        :type stereo: bool
        :rtype: str
    """

    # for now, just assert that we have no multi-component strings with
    # hardcoded parts -- these are guaranteed to fail
    ichs = automol.inchi.split(ich)
    if len(ichs) > 1:
        if any(object_from_hardcoded_inchi_by_key('inchi', ich)
               for ich in ichs):
            ref_ichs = []
            for ich_i in ichs:
                ref_ichs.append(recalculate(ich_i))
            ref_ichs.sort()
            ret = automol.inchi.join(ref_ichs)
            return ret
        # raise error.FailedInchiGenerationError

    ret = object_from_hardcoded_inchi_by_key('inchi', ich)
    if ret is None:
        _options = '-SUU' if stereo else ''
        rdm = _rdkit.from_inchi(ich)
        ret = _rdkit.to_inchi(rdm, options=_options, with_aux_info=False)

    return ret


def graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
    gras = [_connected_graph(ich, stereo=stereo) for ich in ichs]
    for idx, gra in enumerate(gras):
        if idx == 0:
            num = 0
        else:
            num = max(map(max, map(automol.graph.atom_keys, gras[:idx]))) + 1
        gras[idx] = automol.graph.transform_keys(gra, num.__add__)
    gra = functools.reduce(automol.graph.union, gras)
    return gra


def _connected_graph(ich, stereo=True):
    """ Generate a molecular graph from an InChI string where
        all all atoms are connected by at least one bond.

        :param ich: InChI string
        :type ich: str
        :param remove_stereo: parameter to include stereochemistry information
        :type remove_stereo: bool
        :rtype: automol molecular graph
    """

    gra = object_from_hardcoded_inchi_by_key('graph', ich)
    if gra is None:
        ich = automol.inchi.standard_form(ich)
        if not stereo or not automol.inchi.has_stereo(ich):
            rdm = _rdkit.from_inchi(ich)
            gra = _rdkit.to_connectivity_graph(rdm)
        else:
            geo = geometry(ich)
            gra = automol.convert.geom.graph(geo)

    gra = automol.graph.implicit(gra)
    return gra


def smiles(ich):
    """ Convert a SMILES string into an InChI string.

        :param smi: SMILES string
        :type smi: str
        :rtype: str
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
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
        ich = automol.inchi.standard_form(ich)
        rdm = _rdkit.from_inchi(ich)
        smi = _rdkit.to_smiles(rdm)
    return smi


def inchi_key(ich):
    """ InChI => InChIKey
    """
    ick = _rdkit.inchi_to_inchi_key(ich)
    return ick


def formula(ich):
    """ Generate a formula dictionary from an InChI string.

        :param ich: InChI string
        :type ich: str
        :rtype: dict[str: int]
    """

    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
    fmls = list(map(_connected_formula, ichs))
    fml = functools.reduce(automol.formula.join, fmls)
    return fml


def _connected_formula(ich):
    fml = object_from_hardcoded_inchi_by_key('formula', ich)
    if fml is None:
        ich = automol.inchi.standard_form(ich)
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
        if automol.inchi.equivalent(ich, ich_):
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
