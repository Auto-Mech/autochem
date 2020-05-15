""" inchi conversions
"""
import automol.inchi
import automol.graph
from automol import error
from automol.convert import _rdkit
from automol.convert import _pybel
import automol.convert.smiles


INCHI_BLACKLIST = (
    'InChI=1S/C',
    'InChI=1S/B',
    'InChI=1S/N',
    'InChI=1S/CH/h1H',
    'InChI=1S/CF/c1-2',
    'InChI=1S/CCl/c1-2',
    'InChI=1S/CBr/c1-2',
    'InChI=1S/CI/c1-2',
)


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    smi = smiles(ich)
    ich = automol.convert.smiles.inchi(smi)

    # if forcing stereo, we have to go through RDKit and avoid the blacklist
    if force_stereo:
        def _force_stereo(ich):
            if ich not in INCHI_BLACKLIST:
                _options = '-SUU' if force_stereo else ''
                rdm = _rdkit.from_inchi(ich)
                ich = _rdkit.to_inchi(rdm, options=_options)
            return ich

        ichs = automol.inchi.split(ich)
        ichs = list(map(_force_stereo, ichs))
        ich = automol.inchi.join(ichs)

    return ich


def inchi_key(ich):
    """ InChI => InChIKey
    """
    ick = _rdkit.inchi_to_inchi_key(ich)
    return ick


def geometry(ich):
    """ InChI => geometry
    """
    def _gen1(ich):
        smi = smiles(ich)
        rdm = _rdkit.from_smiles(smi)
        geo, = _rdkit.to_conformers(rdm, nconfs=1)
        return geo

    def _gen2(ich):
        pbm = _pybel.from_inchi(ich)
        geo = _pybel.to_geometry(pbm)
        return geo

    def _gen3(ich):
        gra = _graph_without_stereo(ich)
        geo = automol.graph.heuristic_geometry(gra)
        return geo

    for gen_ in [_gen1, _gen2, _gen3]:
        success = False
        try:
            geo = gen_(ich)
            geo_ich = automol.convert.geom.inchi(geo)
            if automol.inchi.same_connectivity(ich, geo_ich) and (
                    not automol.inchi.has_stereo(ich) or
                    automol.inchi.equivalent(ich, geo_ich)):
                success = True
                break
        except (RuntimeError, TypeError, ValueError):
            continue

    if not success:
        raise error.FailedGeometryGenerationError

    return geo


def conformers(ich, nconfs):
    """ InChI => conformers (several geometries)
    """
    try:
        smi = smiles(ich)
        rdm = _rdkit.from_smiles(smi)
        geos = _rdkit.to_conformers(rdm, nconfs=nconfs)
    except (RuntimeError, TypeError, ValueError):
        raise error.FailedGeometryGenerationError

    return geos


def graph(ich, no_stereo=False):
    """ InChI => graph
    """
    gra = _graph_without_stereo(ich)

    if automol.inchi.has_stereo(ich) and not no_stereo:
        geo = geometry(ich)
        gra = automol.graph.set_stereo_from_geometry(gra, geo)

    gra = automol.graph.implicit(gra)
    return gra


def smiles(ich):
    """ InChI => SMILES
    """
    # using OpenBabel because RDKit fails on molecules with unusual valences
    if automol.inchi.has_isotope(ich):
        rdm = _rdkit.from_inchi(ich)
        smi = _rdkit.to_smiles(rdm)
    else:
        pbm = _pybel.from_inchi(ich)
        smi = _pybel.to_smiles(pbm)
    return smi


def formula(ich):
    """ InChI => formula
    """
    smi = smiles(ich)

    rdm = _rdkit.from_smiles(smi)
    fml = _rdkit.to_formula(rdm)
    return fml


def _graph_without_stereo(ich):
    smi = smiles(ich)

    rdm = _rdkit.from_smiles(smi)
    gra = _rdkit.to_connectivity_graph(rdm)
    return gra
