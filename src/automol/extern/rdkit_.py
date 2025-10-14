"""RDKit interface"""

import numbers

import rdkit
from rdkit import RDLogger
from rdkit.Chem import AllChem, Draw, rdDistGeom

from .. import util
from ..geom import base as geom_base
from ..graph import base as graph_base

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)

ATOM_STEREO_TAG_FROM_BOOL = {
    None: rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    False: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    True: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
}
ATOM_STEREO_BOOL_FROM_TAG = dict(map(reversed, ATOM_STEREO_TAG_FROM_BOOL.items()))

BOND_STEREO_TAG_FROM_BOOL = {
    None: rdkit.Chem.rdchem.BondStereo.STEREONONE,
    False: rdkit.Chem.rdchem.BondStereo.STEREOZ,
    True: rdkit.Chem.rdchem.BondStereo.STEREOE,
}
BOND_STEREO_BOOL_FROM_TAG = dict(map(reversed, BOND_STEREO_TAG_FROM_BOOL.items()))


def turn_3d_visualization_on():
    """Turn 3D drawing in RDKit on"""
    if hasattr(Draw, "IPythonConsole"):
        Draw.IPythonConsole.ipython_3d = True


def turn_3d_visualization_off():
    """Turn 3D drawing in RDKit on"""
    if hasattr(Draw, "IPythonConsole"):
        Draw.IPythonConsole.ipython_3d = False


# inchi
def from_inchi(ich, print_debug=False):
    """Generate an RDKit molecule object from an InChI string.

    :param ich: InChI string
    :type ich: str
    :param print_debug: control the printing of a debug message
    :type print_debug: bool
    :rtype: RDKit molecule object
    """

    rdm = rdkit.Chem.inchi.MolFromInchi(ich)
    rdm = rdkit.Chem.AddHs(rdm)

    if rdm is None and print_debug:
        print(f"rdm fails for {ich} by returning {rdm}")

    return rdm


def to_inchi(rdm, options="", with_aux_info=False):
    """Generate an InChI string from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param options:
    :type options: str
    :param with_aux_info: include auxiliary information
    :type with_aux_info: bool
    :rtype: str
    """

    if with_aux_info:
        ret = rdkit.Chem.inchi.MolToInchiAndAuxInfo(rdm, options=options)
    else:
        ret = rdkit.Chem.inchi.MolToInchi(rdm, options=options)

    return ret


# smiles
def from_smiles(smi, print_debug=False):
    """Generate an RDKit molecule object from a SMILES string.

    :param smi: SMILES string
    :type smi: str
    :param print_debug: control the printing of a debug message
    :type print_debug: bool
    :rtype: RDKit molecule object
    """

    rdm = rdkit.Chem.MolFromSmiles(smi)
    if rdm is None and print_debug:
        print(f"rdm fails for {smi} by returning {rdm}")

    return rdm


def to_smiles(rdm):
    """Generate a SMILES string from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :rtype: str
    """
    return rdkit.Chem.MolToSmiles(rdm)


# geometry
def from_geometry_with_graph(geo, gra):
    """Generate an RDKit molecule object from a molecular geometry.

    :param geo: automol geometry data structure
    :type geo: str
    :param geo: automol graph data structure
    :type geo: str
    :rtype: RDKit molecule object
    """
    natms = geom_base.count(geo)
    rdm = from_graph(gra, stereo=False, exp=True)
    rdc = rdkit.Chem.Conformer(natms)
    rdc.Set3D(True)
    xyzs = geom_base.coordinates(geo)
    for i, xyz in enumerate(xyzs):
        rdc.SetAtomPosition(i, rdkit.Geometry.Point3D(*xyz))
    rdm.AddConformer(rdc)
    return rdm


def to_geometry(rdm, seed: int = -1):
    """Generate a molecular geometry from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :rtype: automol geometry data structure
    """
    rdm = rdkit.Chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    natms = len(rdm.GetAtoms())
    if natms == 1:
        syms = [str(atms[0].GetSymbol()).title()]
        xyzs = [(0.0, 0.0, 0.0)]
    else:
        ps = rdDistGeom.ETKDGv3()
        ps.randomSeed = seed
        ps.maxIterations = 10000

        ret = AllChem.EmbedMolecule(rdm, ps)
        if ret < 0:
            ps.useBasicKnowledge = False
            ret = AllChem.EmbedMolecule(rdm, ps)
        if ret < 0:
            ps.useRandomCoords = False
            ret = AllChem.EmbedMolecule(rdm, ps)
        if ret < 0:
            return None

        AllChem.MMFFOptimizeMolecule(rdm)
        syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
        xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
    geo = geom_base.from_data(syms, xyzs, angstrom=True)

    return geo


def canonicalize_geometry(rdm):
    """Canonicalize a geometry using RDKit

    :param rdm: An RDKit molecule
    :return: The geometry
    """
    rdkit.Chem.rdMolTransforms.CanonicalizeMol(rdm)
    atms = rdm.GetAtoms()
    syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
    xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
    geo = geom_base.from_data(syms, xyzs, angstrom=False)
    return geo


# molfile
def from_molfile(mfl, print_debug=False):
    """Generate an RDKit molecule object from a MOLFile string.

    :param mfl: MOLFile block string
    :type mfl: str
    :param print_debug: control the printing of a debug message
    :type print_debug: bool
    :rtype: RDKit molecule object
    """

    rdm = rdkit.Chem.rdmolfiles.MolFromMolBlock(mfl, removeHs=False)
    if rdm is None and print_debug:
        print(f"Warning: rdm fails for {mfl} by returning {rdm}")

    return rdm


def to_molfile(rdm):
    """Generate a MOLFile string from an RDKit molecule object

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :return: MOLFile block string
    :rtype: str
    """
    mfl = rdkit.Chem.MolToMolBlock(rdm)
    return mfl


# smarts (for reactions)
def from_smarts(smr, print_debug=False):
    """Generate an RDKit reaction object from a SMARTS string.

    :param smr: SMARTS string
    :type smr: str
    :param print_debug: control the printing of a debug message
    :type print_debug: bool
    :rtype: RDKit reaction object
    """

    rdm = AllChem.ReactionFromSmarts(smr)
    if rdm is None and print_debug:
        print(f"rdm fails for {smr} by returning {rdm}")

    return rdm


# inchi key
def inchi_to_inchi_key(ich):
    """Convert an InChI string into an InChIKey using RDKit.

    :param ich: InChI string
    :type ich: str
    :rtype: str
    """
    return rdkit.Chem.inchi.InchiToInchiKey(ich)


# formula
def to_formula(rdm):
    """Generate a molecular formula dictionary from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :rtype: dict[str:int]
    """

    rdm = rdkit.Chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    syms = [rda.GetSymbol() for rda in atms]
    fml = util.formula_from_symbols(syms)

    return fml


# graph
BOND_ORDER_DCT = {
    1: rdkit.Chem.BondType.SINGLE,
    1.5: rdkit.Chem.BondType.ONEANDAHALF,
    2: rdkit.Chem.BondType.DOUBLE,
    3: rdkit.Chem.BondType.TRIPLE,
}
BOND_TYPE_DCT = dict(map(reversed, BOND_ORDER_DCT.items()))


def from_graph(
    gra, stereo=False, exp=False, local_stereo=False, label=False, label_dct=None
):
    """Generate an RDKit rdmecule object from a connected rdmecular graph

    :param gra: A molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param local_stereo: Does the graph have local stereo assignments? defaults to True
    :type local_stereo: bool, optional
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    gra = graph_base.without_bonds_by_orders(gra, ords=[0], skip_dummies=False)
    gra = graph_base.smiles_graph(gra, res_stereo=True, exp=exp, dummy=True)
    rdm, idx_from_key = _from_graph_without_stereo(
        gra, label=label, label_dct=label_dct
    )

    # If there's not stereo, return early
    if not stereo or not graph_base.has_stereo(gra):
        return rdm

    # Otherwise, handle stereo
    gra = gra if local_stereo else graph_base.to_local_stereo(gra)
    egra = graph_base.explicit(gra)
    bad_bkeys = graph_base.bad_stereo_bond_keys_from_kekule(egra)
    egra = graph_base.without_stereo(egra, bnd_keys=bad_bkeys)
    erdm, idx_from_key = _from_graph_without_stereo(egra)
    key_from_idx = dict(map(reversed, idx_from_key.items()))
    nkeys_dct = graph_base.stereocenter_candidates(egra, strict=False)

    # Set atom stereo
    atm_ste_dct = graph_base.atom_stereo_parities(egra)
    for rda in rdm.GetAtoms():
        idx = rda.GetIdx()
        key = key_from_idx[idx]
        par0 = atm_ste_dct[key]
        if par0 is not None:
            erda = erdm.GetAtoms()[idx]

            nkeys0 = nkeys_dct[key]
            nidxs1 = [b.GetOtherAtomIdx(erda.GetIdx()) for b in erda.GetBonds()]
            nkeys1 = list(map(key_from_idx.__getitem__, nidxs1))

            par1 = util.is_odd_permutation(nkeys0, nkeys1) ^ par0
            rda.SetChiralTag(ATOM_STEREO_TAG_FROM_BOOL[par1])

    # Set bond stereo
    bnd_ste_dct = graph_base.bond_stereo_parities(gra)
    for rdb in rdm.GetBonds():
        idxs = (rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx())
        keys = tuple(map(key_from_idx.get, idxs))
        bkey = frozenset(keys)
        par = bnd_ste_dct[bkey]
        if par is not None:
            bnkeys = nkeys_dct[bkey]
            nkeys1, nkeys2 = [bnkeys[keys.index(k)] for k in sorted(keys)]

            nidx1 = idx_from_key[nkeys1[-1]]
            nidx2 = idx_from_key[nkeys2[-1]]

            rdb.SetStereo(BOND_STEREO_TAG_FROM_BOOL[par])
            rdb.SetStereoAtoms(nidx1, nidx2)

    return rdm


def _from_graph_without_stereo(gra, label=False, label_dct=None):
    """Generate an RDKit rdmecule object from a connected rdmecular graph

    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    if label_dct is not None:
        label = True

    # Add atoms
    keys = sorted(graph_base.atom_keys(gra))
    symb_dct = graph_base.atom_symbols(gra, dummy_symbol="He")
    rad_dct = graph_base.atom_unpaired_electrons(gra, bond_order=True)

    erdm = rdkit.Chem.EditableMol(rdkit.Chem.Mol())
    for key in keys:
        rda = rdkit.Chem.Atom(symb_dct[key])
        rda.SetNumRadicalElectrons(rad_dct[key])
        erdm.AddAtom(rda)

    # Get key <=> index mappings
    idx_from_key = dict(map(reversed, enumerate(keys)))

    # Add bonds
    bkeys = sorted(graph_base.bond_keys(gra), key=sorted)
    ord_dct = graph_base.bond_orders(gra)
    for bkey in bkeys:
        idx1, idx2 = map(idx_from_key.__getitem__, bkey)
        erdm.AddBond(idx1, idx2, BOND_ORDER_DCT[ord_dct[bkey]])

    rdm = erdm.GetMol()

    if label:
        label_dct = {k: k for k in keys} if label_dct is None else label_dct
        # Re-index the label dict to use indices
        label_dct = util.dict_.transform_keys(label_dct, idx_from_key.get)
        for idx, rda in enumerate(rdm.GetAtoms()):
            if idx in label_dct and isinstance(label_dct[idx], numbers.Integral):
                rda.SetProp("molAtomMapNumber", str(label_dct[idx]))

    rdm.UpdatePropertyCache()
    rdkit.Chem.SanitizeMol(rdm)
    return rdm, idx_from_key


def to_graph(exp_rdm, stereo=True, order=False):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param stereo: Include stereochemistry information?, defaults to True
    :type stereo: bool, optional
    :param order: Include bond orders?, defaults to False
    :type order: bool, optional
    :rtype: automol molecular graph object
    """
    gra = _to_graph_without_stereo(exp_rdm, order=order)

    if not stereo:
        return gra

    exp_rdm = rdkit.Chem.AddHs(exp_rdm)
    exp_gra = _to_graph_without_stereo(exp_rdm)
    nkeys_dct = graph_base.stereocenter_candidates(exp_gra, strict=False)

    # Assign atom stereo
    for exp_rda in exp_rdm.GetAtoms():
        par0 = ATOM_STEREO_BOOL_FROM_TAG[exp_rda.GetChiralTag()]
        if par0 is not None:
            key = exp_rda.GetIdx()

            nkeys0 = nkeys_dct[key]
            nkeys1 = [b.GetOtherAtomIdx(exp_rda.GetIdx()) for b in exp_rda.GetBonds()]

            par1 = util.is_odd_permutation(nkeys0, nkeys1) ^ par0
            gra = graph_base.set_atom_stereo_parities(gra, {key: par1})

    # Assign bond stereo
    for exp_rdb in exp_rdm.GetBonds():
        keys = (exp_rdb.GetBeginAtomIdx(), exp_rdb.GetEndAtomIdx())
        bkey = frozenset(keys)

        par0 = BOND_STEREO_BOOL_FROM_TAG[exp_rdb.GetStereo()]
        if par0 is not None:
            bnkeys = nkeys_dct[bkey]
            nkeys1, nkeys2 = [bnkeys[keys.index(k)] for k in sorted(keys)]

            nkey1 = nkeys1[-1]
            nkey2 = nkeys2[-1]
            nkey1_, nkey2_ = exp_rdb.GetStereoAtoms()

            par1 = (nkey1 != nkey1_) ^ (nkey2 != nkey2_) ^ par0
            gra = graph_base.set_bond_stereo_parities(gra, {bkey: par1})

    gra = graph_base.from_local_stereo(gra)
    return gra


def _to_graph_without_stereo(rdm, order=False):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param order: Include bond orders?, defaults to False
    :type order: bool, optional
    :rtype: automol molecular graph object
    """

    def _get_order(rdb):
        return BOND_TYPE_DCT[rdb.GetBondType()] if order else 1

    rdm.UpdatePropertyCache()
    rdkit.Chem.SanitizeMol(rdm)

    rdas = rdm.GetAtoms()
    rdbs = rdm.GetBonds()
    sym_dct = {a.GetIdx(): a.GetSymbol() for a in rdas}
    hyd_dct = {b.GetIdx(): b.GetImplicitValence() for b in rdas}
    ord_dct = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): _get_order(b) for b in rdbs}
    gra = graph_base.from_data(
        atm_symb_dct=sym_dct,
        bnd_keys=ord_dct.keys(),
        atm_imp_hyd_dct=hyd_dct,
        bnd_ord_dct=ord_dct,
    )
    return gra


def to_connectivity_graph(rdm):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :rtype: automol molecular graph object
    """

    rdm = rdkit.Chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    bnds = rdm.GetBonds()
    syms = [rda.GetSymbol() for rda in atms]
    idx = {rda.GetIdx(): idx for idx, rda in enumerate(atms)}
    bnds = [(idx[rdb.GetBeginAtomIdx()], idx[rdb.GetEndAtomIdx()]) for rdb in bnds]
    sym_dct = dict(enumerate(syms))
    gra = graph_base.from_data(atm_symb_dct=sym_dct, bnd_keys=bnds)

    return gra


# draw operations
def to_svg_string(rdm, image_size=300):
    """Convert the RDKit molecule to an SVG string"""
    rdd = Draw.MolDraw2DSVG(image_size, image_size)
    # opts = rdd.drawOptions()
    # opts.maxFontSize = 80
    # opts.minFontSize = 30
    # opts.bondLineWidth = 4
    rdd.DrawMolecule(rdm)
    rdd.FinishDrawing()
    svg_str = rdd.GetDrawingText()
    return svg_str


def to_grid_svg_string(rdms, image_size=300):
    """Convert a sequence of RDKit molecules to an SVG string"""
    num = len(rdms)
    rdd = Draw.MolDraw2DSVG(
        image_size * num, image_size // 2, image_size, image_size // 2
    )
    opts = rdd.drawOptions()
    # opts.maxFontSize = 80
    # opts.minFontSize = 30
    # opts.bondLineWidth = 7
    opts.padding = 0
    rdd.DrawMolecules(list(rdms))
    rdd.FinishDrawing()
    svg_str = rdd.GetDrawingText()
    return svg_str


def draw(rdm, filename=None, highlight_radicals=False, image_size=600):
    """Convert the RdKit molecule object to a PNG image."""
    cdraw = Draw.rdMolDraw2D.MolDraw2DCairo(image_size, image_size)
    opts = cdraw.drawOptions()
    # opts.maxFontSize = 90
    # opts.minFontSize = 70
    # opts.bondLineWidth = 20
    opts.clearBackground = False
    highlights = ()
    highlight_radii = {}
    highlight_colors = {}
    if highlight_radicals:
        atms = rdm.GetAtoms()
        for i, _ in enumerate(atms):
            if rdm.GetAtomWithIdx(i).GetNumRadicalElectrons() > 0:
                highlights += (i,)
                highlight_radii[i] = 0.8
                highlight_colors[i] = (0.67, 0.85, 0.9)
    cdraw.DrawMolecule(
        rdm,
        highlightAtoms=highlights,
        highlightAtomRadii=highlight_radii,
        highlightAtomColors=highlight_colors,
    )
    opts.includeRadicals = False
    cdraw.FinishDrawing()
    draw_out = cdraw.GetDrawingText()
    if filename:
        cdraw.WriteDrawingText(filename)
    return draw_out


def draw_grid(rdms, img_per_row=3, sub_img_size=(200, 200), legends=None):
    """Draw a grid"""
    # Set various options for drawing the grid
    if legends is None:
        legends = tuple(i + 1 for i in range(len(rdms)))
    else:
        assert len(rdms) == len(
            legends
        ), "User provided ichs and legends are not the same length"

    return Draw.MolsToGridImage(
        rdms, molsPerRow=img_per_row, subImgSize=sub_img_size, legends=legends
    )


def draw_mult(rdms, sub_img_size=(200, 200)):
    """Draw multiple"""
    return Draw.MolsToImage(rdms, subImgSize=sub_img_size, legends=None)
