""" RDKit interface
"""
import pyparsing as pp
import rdkit
from lxml import etree
from rdkit import RDLogger
from rdkit.Chem import AllChem, Draw

import automol.geom.base
import automol.graph.base
from automol import util

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


def turn_3d_visualization_on():
    """Turn 3D drawing in RDKit on"""
    rdkit.Chem.Draw.IPythonConsole.ipython_3d = True


def turn_3d_visualization_off():
    """Turn 3D drawing in RDKit on"""
    rdkit.Chem.Draw.IPythonConsole.ipython_3d = False


# inchi
def from_inchi(ich, print_debug=False):
    """Generate an RDKit molecule object from an InChI string.

    :param ich: InChI string
    :type ich: str
    :param print_debug: control the printing of a debug message
    :type print_debug: bool
    :rtype: RDKit molecule object
    """

    rdm = rdkit.Chem.inchi.MolFromInchi(ich, treatWarningAsError=False)
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
    natms = automol.geom.base.count(geo)
    rdm = from_graph(gra)
    rdc = rdkit.Chem.Conformer(natms)
    rdc.Set3D(True)
    xyzs = automol.geom.base.coordinates(geo)
    for i, xyz in enumerate(xyzs):
        rdc.SetAtomPosition(i, rdkit.Geometry.Point3D(*xyz))
    rdm.AddConformer(rdc)
    return rdm


def to_geometry(rdm):
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
        AllChem.EmbedMolecule(rdm)
        AllChem.MMFFOptimizeMolecule(rdm)
        syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
        xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
    geo = automol.geom.base.from_data(syms, xyzs, angstrom=True)

    return geo


def to_conformers(rdm, nconfs):
    """Generate molecular geometries for a set of conformers
    from am RDKit molecule object.

    Currently not removing redundant conformers.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :param nconfs: number of conformers to generate
    :type nconfs: int
    :rtype: automol geometry data structure
    """

    rdm = rdkit.Chem.AddHs(rdm)
    atms = rdm.GetAtoms()
    natms = len(rdm.GetAtoms())
    geos = []
    if natms == 1:
        syms = [str(atms[0].GetSymbol()).title()]
        xyzs = [(0.0, 0.0, 0.0)]
        geos.append(automol.geom.base.from_data(syms, xyzs, angstrom=True))
    else:
        cids = AllChem.EmbedMultipleConfs(rdm, numConfs=nconfs)
        res = AllChem.MMFFOptimizeMoleculeConfs(rdm)
        energies = list(zip(*res))[1]
        for cid in cids:
            syms = tuple(str(rda.GetSymbol()).title() for rda in atms)
            xyzs = tuple(map(tuple, rdm.GetConformer(cid).GetPositions()))
            geos.append(automol.geom.base.from_data(syms, xyzs, angstrom=True))
        # Sort geometries using the energies
        geos = [x for _, x in sorted(zip(energies, geos), key=lambda pair: pair[0])]

    return geos


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


def from_graph(gra, stereo=False, label=False, label_dct=None):
    """Generate an RDKit rdmecule object from a connected rdmecular graph

    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    if stereo:
        raise NotImplementedError("Stereo is currently not implemented!")

    if label_dct is not None:
        label = True

    gra = automol.graph.base.without_bonds_by_orders(gra, ords=[0], skip_dummies=False)
    kgr = automol.graph.base.kekule(gra, max_stereo_overlap=True)
    keys = sorted(automol.graph.base.atom_keys(kgr))
    key_map = dict(map(reversed, enumerate(keys)))
    symb_dct = automol.graph.base.atom_symbols(kgr, dummy_symbol="He")
    rad_dct = automol.graph.base.atom_unpaired_electrons(kgr, bond_order=True)
    if label:
        label_dct = {k: k for k in keys} if label_dct is None else label_dct
        # Re-index the label dict to use indices
        label_dct = util.dict_.transform_keys(label_dct, key_map.__getitem__)

    def idx_bond_key_(bnd_key):
        key1, key2 = bnd_key
        return frozenset({key_map[key1], key_map[key2]})

    bnd_keys = list(map(idx_bond_key_, automol.graph.base.bond_keys(kgr)))
    ord_dct = util.dict_.transform_keys(
        automol.graph.base.bond_orders(kgr), idx_bond_key_
    )

    erdm = rdkit.Chem.EditableMol(rdkit.Chem.Mol())
    for key in keys:
        atm = rdkit.Chem.Atom(symb_dct[key])
        atm.SetNumRadicalElectrons(rad_dct[key])
        erdm.AddAtom(atm)

    for bnd_key in bnd_keys:
        erdm.AddBond(*bnd_key, BOND_ORDER_DCT[ord_dct[bnd_key]])

    rdm = erdm.GetMol()

    if label:
        for idx, rda in enumerate(rdm.GetAtoms()):
            if idx in label_dct:
                rda.SetProp("molAtomMapNumber", str(label_dct[idx]))

    rdm.UpdatePropertyCache()
    return rdm


def to_graph(rdm):
    """Generate a connectivity graph from an RDKit molecule object.

    :param rdm: molecule object
    :type rdm: RDKit molecule object
    :rtype: automol molecular graph object
    """

    rdm.UpdatePropertyCache()
    atms = rdm.GetAtoms()
    bnds = rdm.GetBonds()
    sym_dct = {rda.GetIdx(): rda.GetSymbol() for rda in atms}
    hyd_dct = {rda.GetIdx(): rda.GetImplicitValence() for rda in atms}
    ord_dct = {
        (rdb.GetBeginAtomIdx(), rdb.GetEndAtomIdx()): BOND_TYPE_DCT[rdb.GetBondType()]
        for rdb in bnds
    }
    gra = automol.graph.base.from_data(
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
    gra = automol.graph.base.from_data(atm_symb_dct=sym_dct, bnd_keys=bnds)

    return gra


# draw operations
def to_svg_string(rdm, image_size=300):
    """Convert the RDKit molecule to an SVG string"""
    rdd = Draw.MolDraw2DSVG(image_size, image_size)
    opts = rdd.drawOptions()
    opts.maxFontSize = 80
    opts.minFontSize = 30
    opts.bondLineWidth = 7
    rdd.DrawMolecule(rdm)
    rdd.FinishDrawing()
    svg_str = rdd.GetDrawingText()
    return svg_str


def to_reagent_svg_strings(rdr, image_size=300):
    """Convert RDKit reaction to a pair of reactant and product SVG strings"""
    rdd = Draw.MolDraw2DSVG(image_size, image_size)
    opts = rdd.drawOptions()
    opts.maxFontSize = 30
    opts.minFontSize = 12
    opts.bondLineWidth = 2
    rdd.DrawReaction(rdr)
    rdd.FinishDrawing()
    svg_str = rdd.GetDrawingText()
    r_svg_str = _reagent_svg_string(
        svg_str, product=False, width=image_size, height=image_size
    )
    p_svg_str = _reagent_svg_string(
        svg_str, product=True, width=image_size, height=image_size
    )
    return r_svg_str, p_svg_str


# SVG parsers to help extract info rom the reaction string
def _reagent_svg_string(svg_str, product=False, width=300, height=300, padding=15):
    svg_root = _etree_from_svg_string(svg_str)

    # 1. Remove elements not on the left/right half
    for child in svg_root:
        if product:
            xmin = _parse_x_min(child.attrib["d"]) if "d" in child.attrib else None
            if xmin is None or xmin < width / 2.0:
                svg_root.remove(child)
        else:
            xmax = _parse_x_max(child.attrib["d"]) if "d" in child.attrib else None
            if xmax is None or xmax > width / 2.0:
                svg_root.remove(child)

    # 2. For products, remove the two remaining leftmost elements (arrow head parts)
    if product:
        children = sorted(svg_root, key=lambda c: _parse_x_bounds(c.attrib["d"])[0])
        svg_root.remove(children.pop(0))
        svg_root.remove(children.pop(0))

    # 3. Determine the bounding box
    xmax = ymax = 0
    xmin = width
    ymin = height
    for child in svg_root:
        xmin_, xmax_ = _parse_x_bounds(child.attrib["d"])
        ymin_, ymax_ = _parse_y_bounds(child.attrib["d"])
        xmin = min(xmin, xmin_)
        xmax = max(xmax, xmax_)
        ymin = min(ymin, ymin_)
        ymax = max(ymax, ymax_)

    # 4. Adjust the viewBox accordingly, with the desired padding
    xmin -= padding
    ymin -= padding
    xmax += padding
    ymax += padding
    svg_root.attrib["viewBox"] = f"{xmin} {ymin} {xmax-xmin} {ymax-ymin}"
    return _etree_to_svg_string(svg_root)


def _parse_xy_coordinates(string):
    command = pp.Suppress(pp.Char("ML"))
    number = pp.Combine(pp.Word(pp.nums) + pp.Opt("." + pp.Word(pp.nums)))
    parser = command + pp.Group(number + pp.Suppress(pp.Opt(",")) + number)
    coords = [tuple(map(float, r[0])) for r in parser.search_string(string)]
    return coords


def _parse_x_coordinates(string):
    coords = _parse_xy_coordinates(string)
    if not coords:
        return ()

    xvals, _ = zip(*coords)
    return xvals


def _parse_y_coordinates(string):
    coords = _parse_xy_coordinates(string)
    if not coords:
        return ()

    _, yvals = zip(*coords)
    return yvals


def _parse_x_bounds(string):
    xvals = _parse_x_coordinates(string)
    xmin = min(xvals, default=None)
    xmax = max(xvals, default=None)
    return xmin, xmax


def _parse_y_bounds(string):
    yvals = _parse_y_coordinates(string)
    ymin = min(yvals, default=None)
    ymax = max(yvals, default=None)
    return ymin, ymax


def _parse_x_min(string):
    return _parse_x_bounds(string)[0]


def _parse_x_max(string):
    return _parse_x_bounds(string)[1]


def _etree_from_svg_string(svg_str):
    parser = etree.XMLParser(ns_clean=True, recover=True, encoding="utf-8")
    svg_root = etree.fromstring(svg_str.encode("utf-8"), parser=parser)
    return svg_root


def _etree_to_svg_string(svg_root):
    return etree.tostring(svg_root, pretty_print=True).decode("utf-8")


# end SVG parsers


def draw(rdm, filename=None, highlight_radicals=False, image_size=600):
    """Convert the RdKit molecule object to a PNG image."""
    cdraw = Draw.rdMolDraw2D.MolDraw2DCairo(image_size, image_size)
    opts = cdraw.drawOptions()
    opts.maxFontSize = 90
    opts.minFontSize = 70
    opts.bondLineWidth = 20
    opts.clearBackground = False
    # s.setBackgroundColour = 'black'
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
    draw_out = draw.GetDrawingText()
    if filename:
        draw.WriteDrawingText(filename)
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
