""" rdkit interface
"""
from rdkit import RDLogger
import rdkit.Chem as _rd_chem
import rdkit.Chem.AllChem as _rd_all_chem
import automol.create
from automol.convert import _util

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


# geometry
def to_geometry(rdm):
    """ cartesian geometry from an rdkit molecule object
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    natms = len(rdm.GetAtoms())
    if natms == 1:
        syms = [str(atms[0].GetSymbol()).title()]
        xyzs = [(0., 0., 0.)]
    else:
        _rd_all_chem.EmbedMolecule(rdm)
        _rd_all_chem.MMFFOptimizeMolecule(rdm)
        syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
        xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
    geo = automol.create.geom.from_data(syms, xyzs, angstrom=True)
    return geo


def to_conformers(rdm, nconfs):
    """ list of cartesian geometries for conformers
        from an rdkit molecule object
        currently not removing redundant conformers
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    natms = len(rdm.GetAtoms())
    geos = []
    if natms == 1:
        syms = [str(atms[0].GetSymbol()).title()]
        xyzs = [(0., 0., 0.)]
        geos.append(
            automol.create.geom.from_data(syms, xyzs, angstrom=True))
    else:
        cids = _rd_all_chem.EmbedMultipleConfs(rdm, numConfs=nconfs)
        res = _rd_all_chem.MMFFOptimizeMoleculeConfs(rdm)
        energies = list(zip(*res))[1]
        for cid in cids:
            syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
            xyzs = tuple(map(tuple, rdm.GetConformer(cid).GetPositions()))
            geos.append(
                automol.create.geom.from_data(syms, xyzs, angstrom=True))
        # Sort geometries using the energies
        geos = [x for _, x in sorted(zip(energies, geos), key=lambda pair: pair[0])]
    return geos


# inchi
def from_inchi(ich):
    """ rdkit molecule object from an InChI string
    """
    rdm = _rd_chem.inchi.MolFromInchi(ich, treatWarningAsError=False)
    assert rdm is not None
    return rdm


def to_inchi(rdm, options='', with_aux_info=False):
    """ InChI string from an rdkit molecule object
    """
    if with_aux_info:
        ret = _rd_chem.inchi.MolToInchiAndAuxInfo(rdm, options=options)
    else:
        ret = _rd_chem.inchi.MolToInchi(rdm, options=options)
    return ret


# smiles
def from_smiles(smi):
    """ rdkit molecule object from a SMILES string
    """
    # print('rd smi', smi)
    rdm = _rd_chem.MolFromSmiles(smi)
    assert rdm is not None
    return rdm


def to_smiles(rdm):
    """ SMILES string from an rdkit molecule object
    """
    smi = _rd_chem.MolToSmiles(rdm)
    return smi


# connectivity graph
def to_connectivity_graph(rdm):
    """ connection graph from an rdkit molecule object
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    bnds = rdm.GetBonds()
    syms = [rda.GetSymbol() for rda in atms]
    idx = {rda.GetIdx(): idx for idx, rda in enumerate(atms)}
    bnds = [(idx[rdb.GetBeginAtomIdx()], idx[rdb.GetEndAtomIdx()])
            for rdb in bnds]
    sym_dct = dict(enumerate(syms))
    gra = automol.create.graph.from_data(atom_symbols=sym_dct, bond_keys=bnds)
    return gra


# molfile
def from_molfile(mfl):
    """ rdkit molecule object from a mol block string
    """
    rdm = _rd_chem.rdmolfiles.MolFromMolBlock(mfl, removeHs=False)
    assert rdm is not None
    return rdm


# inchi key
def inchi_to_inchi_key(ich):
    """ InChI-Key from an InChI string
    """
    ick = _rd_chem.inchi.InchiToInchiKey(ich)
    return ick


# formula
def to_formula(rdm):
    """ formula from an rdkit molecule object
    """
    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    syms = [rda.GetSymbol() for rda in atms]
    fml = _util.formula(syms)
    return fml
