""" RDKit interface
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
    """ Generate a molecular geometry from an RDKit molecule object.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :rtype: automol geometry data structure
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
    """ Generate molecular geometries for a set of conformers
        from am RDKit molecule object.
   
        Currently not removing redundant conformers.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :param nconfs: number of conformers to generate
        :type nconfs: int
        :rtype: automol geometry data structure
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
def from_inchi(ich, print_debug=True):
    """ Generate an RDKit molecule object from an InChI string.

        :param ich: InChI string
        :type ich: str
        :param print_debug: control the printing of a debug message
        :type print_debug: bool
        :rtype: RDKit molecule object
    """

    rdm = _rd_chem.inchi.MolFromInchi(ich, treatWarningAsError=False)
    if rdm is None and print_debug:
        print('rdm fails for {} by returning {}'.format(ich, rdm))

    return rdm


def to_inchi(rdm, options='', with_aux_info=False):
    """ Generate an InChI string from an RDKit molecule object.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :param options: 
        :type options: str
        :param with_aux_info: include auxiliary information
        :type with_aux_info: bool
        :rtype: str
    """

    if with_aux_info:
        ret = _rd_chem.inchi.MolToInchiAndAuxInfo(rdm, options=options)
    else:
        ret = _rd_chem.inchi.MolToInchi(rdm, options=options)

    return ret


# smiles
def from_smiles(smi, print_debug=True):
    """ Generate an RDKit molecule object from a SMILES string.

        :param smi: SMILES string
        :type smi: str
        :param print_debug: control the printing of a debug message
        :type print_debug: bool
        :rtype: RDKit molecule object
    """

    rdm = _rd_chem.MolFromSmiles(smi)
    if rdm is not None and print_debug:
        print('rdm fails for {} by returning {}'.format(smi, rdm))

    return rdm


def to_smiles(rdm):
    """ Generate a SMILES string from an RDKit molecule object.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :rtype: str
    """
    return _rd_chem.MolToSmiles(rdm)


# connectivity graph
def to_connectivity_graph(rdm):
    """ Generate a connectivity graph from an RDKit molecule object.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :rtype: automol molecular graph object
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
def from_molfile(mfl, print_debug=True):
    """ Generate an RDKit molecule object from a MOLFile string.

        :param mfl: MOLFile block string
        :type mfl: str
        :param print_debug: control the printing of a debug message
        :type print_debug: bool
        :rtype: RDKit molecule object
    """

    rdm = _rd_chem.rdmolfiles.MolFromMolBlock(mfl, removeHs=False)
    if rdm is None and print_debug:
        print('Warning: rdm fails for {} by returning {}'.format(mfl, rdm))

    return rdm


# inchi key
def inchi_to_inchi_key(ich):
    """ Convert an InChI string into an InChIKey using RDKit.

        :param ich: InChI string
        :type ich: str
        :rtype: str
    """
    return _rd_chem.inchi.InchiToInchiKey(ich)


# formula
def to_formula(rdm):
    """ Generate a molecular formula dictionary from an RDKit molecule object.

        :param rdm: molecule object
        :type rdm: RDKit molecule object
        :rtype: dict[str:int]
    """

    rdm = _rd_chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    syms = [rda.GetSymbol() for rda in atms]
    fml = _util.formula(syms)

    return fml
