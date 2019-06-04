""" inchi conversions
"""
from qcelemental import constants as qcc
from automol import error
import automol.inchi
import automol.convert.geom
from automol.convert import _rdkit
from automol.convert import _pybel

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


def geometry(ich):
    """ InChI => geometry
    """
    if ich in HARDCODED_INCHI_TO_GEOM_DCT:
        geo = HARDCODED_INCHI_TO_GEOM_DCT[ich]
    else:
        try:
            rdm = _rdkit.from_inchi(ich)
            geo = _rdkit.to_geometry(rdm)
            geo_ich = automol.convert.geom.inchi(geo)
            if not automol.inchi.same_connectivity(ich, geo_ich):
                raise error.FailedGeometryGenerationError
        except (error.FailedGeometryGenerationError, RuntimeError):
            pbm = _pybel.from_inchi(ich)
            geo = _pybel.to_geometry(pbm)
            geo_ich = automol.convert.geom.inchi(geo)
            if not automol.inchi.same_connectivity(ich, geo_ich):
                raise error.FailedGeometryGenerationError
    return geo


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdkit.from_inchi(ich)
    ich = _rdkit.to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def graph(ich):
    """ inchi => graph
    """
    if ich in HARDCODED_INCHI_TO_GRAPH_DCT:
        gra = HARDCODED_INCHI_TO_GRAPH_DCT[ich]
    else:
        if not automol.inchi.has_known_stereo_elements(ich):
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
