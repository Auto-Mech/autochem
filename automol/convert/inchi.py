""" inchi conversions
"""
import functools
from qcelemental import constants as qcc
from automol import error
import automol.inchi
import automol.geom
import automol.convert.geom
from automol.convert import _rdkit
from automol.convert import _pybel

HARDCODED_INCHI_LST = [
    'InChI=1S/C',
    'InChI=1S/N',
    'InChI=1S/CH/h1H',
    'InChI=1S/CF/c1-2',
    'InChI=1S/CCl/c1-2',
    'InChI=1S/CBr/c1-2',
    'InChI=1S/CI/c1-2',
]
ANG2BOHR = qcc.conversion_factor('angstrom', 'bohr')
HARDCODED_INCHI_TO_GEOM_DCT = {
    'InChI=1S/C': (('C', (0., 0., 0.)),),
    'InChI=1S/N': (('N', (0., 0., 0.)),),
    'InChI=1S/CH/h1H': (('C', (0., 0., 0.)),
                        ('H', (0., 0., 1.12 * ANG2BOHR))),
    'InChI=1S/CF/c1-2': (('C', (0., 0., 0.)),
                         ('F', (0., 0., 1.27 * ANG2BOHR))),
    'InChI=1S/CCl/c1-2': (('C', (0., 0., 0.)),
                          ('Cl', (0., 0., 1.65 * ANG2BOHR))),
    'InChI=1S/CBr/c1-2': (('C', (0., 0., 0.)),
                          ('Br', (0., 0., 1.8 * ANG2BOHR))),
    'InChI=1S/CI/c1-2': (('C', (0., 0., 0.)),
                         ('I', (0., 0., 1.8 * ANG2BOHR))),
}
HARDCODED_INCHI_TO_GRAPH_DCT = {
    'InChI=1S/C': ({0: ('C', 0, None)}, {}),
    'InChI=1S/N': ({0: ('N', 0, None)}, {}),
    'InChI=1S/CH/h1H': ({0: ('C', 1, None)}, {}),
    'InChI=1S/CF/c1-2': ({0: ('C', 0, None), 1: ('F', 0, None)},
                         {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CCl/c1-2': ({0: ('C', 0, None), 1: ('Cl', 0, None)},
                          {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CBr/c1-2': ({0: ('C', 0, None), 1: ('Br', 0, None)},
                          {frozenset({0, 1}): (1, None)}),
    'InChI=1S/CI/c1-2': ({0: ('C', 0, None), 1: ('I', 0, None)},
                         {frozenset({0, 1}): (1, None)}),
}
HARDCODED_INCHI_TO_SMILES_DCT = {
    'InChI=1S/C': '[C]',
    'InChI=1S/N': '[N]',
    'InChI=1S/CH/h1H': '[CH]',
    'InChI=1S/CF/c1-2': '[C]F',
    'InChI=1S/CCl/c1-2': '[C]Cl',
    'InChI=1S/CBr/c1-2': '[C]Br',
    'InChI=1S/CI/c1-2': '[C]I',
}
HARDCODED_INCHI_TO_FORMULA_DCT = {
    'InChI=1S/C': {'C': 1},
    'InChI=1S/N': {'N': 1},
    'InChI=1S/CH/h1H': {'C': 1, 'H': 1},
    'InChI=1S/CF/c1-2': {'C': 1, 'F': 1},
    'InChI=1S/CCl/c1-2': {'C': 1, 'Cl': 1},
    'InChI=1S/CBr/c1-2': {'C': 1, 'Br': 1},
    'InChI=1S/CI/c1-2': {'C': 1, 'I': 1},
}


def geometry(ich):
    """ InChI => geometry
    """
    # rdkit fails for multi-component inchis, so we split it up and space out
    # the geometries
    ichs = automol.inchi.split(ich)
    geos = list(map(_connected_geometry, ichs))
    geos = [automol.geom.translated(geo, [50. * idx, 0., 0.])
            for idx, geo in enumerate(geos)]
    geo = functools.reduce(automol.geom.join, geos)
    return geo


def _connected_geometry(ich):
    if ich in HARDCODED_INCHI_TO_GEOM_DCT:
        geo = HARDCODED_INCHI_TO_GEOM_DCT[ich]
    else:
        try:
            rdm = _rdkit.from_inchi(ich)
            geo = _rdkit.to_geometry(rdm)
            geo_ich = automol.convert.geom.inchi(geo)
            if not automol.inchi.same_connectivity(ich, geo_ich) or (
                    automol.inchi.has_stereo(ich) and not
                    automol.inchi.equivalent(ich, geo_ich)):
                raise error.FailedGeometryGenerationError
        except (error.FailedGeometryGenerationError, RuntimeError):
            pbm = _pybel.from_inchi(ich)
            geo = _pybel.to_geometry(pbm)
            geo_ich = automol.convert.geom.inchi(geo)
            if not automol.inchi.same_connectivity(ich, geo_ich) or (
                    automol.inchi.has_stereo(ich) and not
                    automol.inchi.equivalent(ich, geo_ich)):
                raise error.FailedGeometryGenerationError
    return geo


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdkit.from_inchi(ich)
    ich = _rdkit.to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def graph(ich, no_stereo=False):
    """ inchi => graph
    """
    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
    gras = [_connected_graph(ich, no_stereo=no_stereo) for ich in ichs]
    for idx, gra in enumerate(gras):
        if idx == 0:
            num = 0
        else:
            num = max(map(max, map(automol.graph.atom_keys, gras[:idx]))) + 1
        gras[idx] = automol.graph.transform_keys(gra, num.__add__)
    gra = functools.reduce(automol.graph.union, gras)
    return gra


def _connected_graph(ich, no_stereo=False):
    if ich in HARDCODED_INCHI_TO_GRAPH_DCT:
        gra = HARDCODED_INCHI_TO_GRAPH_DCT[ich]
    else:
        if no_stereo or not automol.inchi.has_stereo(ich):
            rdm = _rdkit.from_inchi(ich)
            gra = _rdkit.to_connectivity_graph(rdm)
        else:
            geo = geometry(ich)
            gra = automol.convert.geom.graph(geo)

    gra = automol.graph.implicit(gra)
    return gra


def smiles(ich):
    """ InChI => SMILEs
    """
    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
    smis = list(map(_connected_smiles, ichs))
    smi = '.'.join(smis)
    return smi


def _connected_smiles(ich):
    if ich in HARDCODED_INCHI_TO_SMILES_DCT:
        smi = HARDCODED_INCHI_TO_SMILES_DCT[ich]
    else:
        rdm = _rdkit.from_inchi(ich)
        smi = _rdkit.to_smiles(rdm)
    return smi


def inchi_key(ich):
    """ InChI => InChIKey
    """
    ick = _rdkit.inchi_to_inchi_key(ich)
    return ick


def formula(ich):
    """ InChI => formula
    """
    # split it up to handle hard-coded molecules in multi-component inchis
    ichs = automol.inchi.split(ich)
    fmls = list(map(_connected_formula, ichs))
    fml = functools.reduce(automol.formula.join, fmls)
    return fml


def _connected_formula(ich):
    if ich in HARDCODED_INCHI_TO_FORMULA_DCT:
        fml = HARDCODED_INCHI_TO_FORMULA_DCT[ich]
    else:
        rdm = _rdkit.from_inchi(ich)
        fml = _rdkit.to_formula(rdm)
    return fml
