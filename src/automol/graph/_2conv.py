"""graph conversions"""

import functools
import itertools
from collections.abc import Sequence
from typing import Any, Dict, Optional
from xml.etree import ElementTree

import IPython.display as ipd
import more_itertools as mit
import numpy
from phydat import phycon

from .. import error, geom
from ..extern import rdkit_
from ..smiles import base as smiles_base
from ..util import dict_, vector
from ._0embed import clean_geometry, geometry_matches
from .base import (
    align_with_geometry,
    amchi,
    atom_keys,
    atom_stereo_parities,
    connected_components,
    explicit,
    geometry_correct_linear_vinyls,
    geometry_dihedrals_near_value,
    has_dummy_atoms,
    has_stereo,
    inchi_is_bad,
    is_ts_graph,
    nonoverlapping_keys_for_sequence,
    smiles,
    standard_keys,
    stereo_corrected_geometry,
    to_local_stereo,
    ts,
    union_from_sequence,
    without_stereo,
)

NAMESPACES = {"": "http://www.w3.org/2000/svg"}
for namespace in NAMESPACES.items():
    ElementTree.register_namespace(*namespace)


# # conversions
def geometry(
    gra,
    check=True,
    log=False,
    rct_geos: Sequence[object] | None = None,
    geo_idx_dct: Dict[int, int] | None = None,
):
    """Convert a molecular graph to a molecular geometry.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param check: Check stereo and connectivity? defaults to True
    :type check: bool, optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :param rct_geos: For a TS graph, reactant geometries
    :param geo_idx_dct: For a TS graph, specify which graph
        keys correspond to which geometry indices, defaults to None
    :rtype: automol molecular geometry data structure
    """
    is_ts = is_ts_graph(gra)

    # If this is a TS graph and reactant geometries were provided, try to build the TS
    # from them. If this fails, try to build the TS from scratch.
    if is_ts and rct_geos is not None:
        try:
            return ts_geometry_from_reactants(
                gra, rct_geos=rct_geos, geo_idx_dct=geo_idx_dct, check=check, log=log
            )
        except error.FailedGeometryGenerationError:
            pass

    if is_ts:
        tsg = gra
        gra = ts.reactants_graph(tsg)

    gra = explicit(gra)
    gras = connected_components(gra)

    geos = [_connected_geometry(g, check=check, log=log) for g in gras]
    geos = [geom.translate(g, [50.0 * i, 0.0, 0.0]) for i, g in enumerate(geos)]
    geo = functools.reduce(geom.join, geos)

    # Work out how the geometry needs to be re-ordered to match the input graph
    all_keys = itertools.chain(*(sorted(atom_keys(g)) for g in gras))
    gra_key_dct = dict(enumerate(all_keys))

    if is_ts:
        geo_idx_dct = dict(map(reversed, gra_key_dct.items()))

        if log:
            print(f"Building geometry for TS graph:\n{tsg}")
            print(f"... using these reactant geometries:\n{geos}")

        geo = ts_geometry_from_reactants(
            tsg, geos, geo_idx_dct=geo_idx_dct, check=check, log=log
        )
    else:
        geo = geom.reorder(geo, gra_key_dct)

    return geo


def _connected_geometry(gra, check=True, log=False, ntries: int = 4):
    """Generate a geometry for a connected molecular graph.

    :param gra: connected molecular graph
    :type gra: automol graph data structure
    :param check: Check stereo and connectivity? defaults to True
    :type check: bool, optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :rtype: automol molecular geometry
    """
    orig_gra = gra

    # Normalize the graph
    gra = standard_keys(gra)
    gra = to_local_stereo(gra)

    # Determine if stereochemistry and/or vinyl radical groups are present
    stereo = has_stereo(gra)

    # Try geometry generation methods until one works
    for seed in range(ntries):
        try:
            rdm = rdkit_.from_graph(gra, stereo=stereo, exp=True, local_stereo=True)
            geo = rdkit_.to_geometry(rdm, seed=seed)
        except ValueError:
            rdm = rdkit_.from_graph(gra, stereo=False, exp=True)
            geo = rdkit_.to_geometry(rdm, seed=seed)

        if geo is None:
            continue

        if log:
            print(f"Try {seed}...")
            print("Raw geometry:")
            print(geom.round_(geo))

        geo = _clean_and_validate_connected_geometry(
            gra, geo, stereo=stereo, check=check, log=log
        )

        if geo is not None:
            return geo

    raise error.FailedGeometryGenerationError(f"\nseed = {seed}\norig_gra = {orig_gra}")


def _clean_and_validate_connected_geometry(
    gra,
    geo,
    stereo: bool = True,
    local_stereo: bool = True,
    check: bool = True,
    log: bool = False,
):
    """Validate and clen up a connected geometry

    :param gra: Connected molecular graph with standard keys
    :type gra: automol graph data structure
    :param geo: Molecular geometry
    :type geo: automol geom data structure
    :param stereo: Take stereochemistry into consideration? defaults to True
    :param local_stereo: Does the graph have local stereo assignments? defaults to True
    :param check: Check stereo and connectivity? defaults to True
    :param log: Log information to the screen? defaults to False
    """
    gra = gra if stereo else without_stereo(gra)
    gra = gra if local_stereo else to_local_stereo(gra)

    if stereo and geo is not None:
        geo = stereo_corrected_geometry(gra, geo, local_stereo=local_stereo)

    geo = geometry_correct_linear_vinyls(gra, geo)

    if not geometry_matches(gra, geo, stereo=stereo, local_stereo=True):
        geo = clean_geometry(gra, geo, stereo=stereo, local_stereo=True)

    # Remove perfectly planar dihedral angles, to avoid symmetry problems
    geo = perturb_geometry_planar_dihedrals(gra, geo, ang=5.0, degree=True)

    if log:
        print("Cleaned geometry:")
        print(None if geo is None else geom.round_(geo))

    matches = geometry_matches(gra, geo, stereo=stereo, local_stereo=local_stereo)

    return geo if matches or not check else None


def inchi(gra, stereo: bool = True, local_stereo: bool = False):
    """Generate an InChI string from a molecular graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    smi = smiles(gra, stereo=stereo, res_stereo=False, local_stereo=local_stereo)
    rdm = rdkit_.from_smiles(smi)
    ich = rdkit_.to_inchi(rdm)
    return ich


def chi(gra, stereo=True):
    """Generate a ChI string from a molecular graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: ChI string
    :rtype: str
    """
    ret = inchi(gra, stereo=stereo)
    if inchi_is_bad(gra, ret):
        ret = amchi(gra, stereo=stereo)

    return ret


def rdkit_molecule(gra, stereo=True, exp=False, label=False, label_dct=None):
    """Convert a molecular graph to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(gra))

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_off()
    rdm = rdkit_.from_graph(
        gra, stereo=stereo, exp=exp, label=label, label_dct=label_dct
    )
    return rdm


def svg_string(
    gra, stereo=True, exp=False, label=False, label_dct=None, image_size=200
):
    """Get an SVG string for visualizing the graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    :param image_size: The image size, defaults to 150
    :type image_size: int, optional
    :return: The SVG string, in svg+xml format
    :rtype: str
    """
    # If this is not a TS graph, use RDKit directly
    if not is_ts_graph(gra):
        rdm = rdkit_molecule(
            gra, stereo=stereo, exp=exp, label=label, label_dct=label_dct
        )
        return rdkit_.to_svg_string(rdm, image_size=image_size)

    rgra = ts.reactants_graph(gra, stereo=stereo)
    pgra = ts.products_graph(gra, stereo=stereo)
    return reaction_svg_string(
        [rgra],
        [pgra],
        stereo=stereo,
        exp=exp,
        label=label,
        label_dct=label_dct,
        image_size=image_size,
    )


def reaction_svg_string(
    rgras, pgras, stereo=True, exp=False, label=False, label_dct=None, image_size=200
):
    """Get an SVG string for visualizing the graph.

    :param rgra: reactants molecular graph
    :param pgra: products molecular graph
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    :param image_size: The image size, defaults to 150
    :type image_size: int, optional
    :return: The SVG string, in svg+xml format
    :rtype: str
    """
    height = int(image_size)
    width = int(2.5 * image_size)

    #  1. Prepare new SVG canvas with reaction arrow
    root = ElementTree.Element(
        "svg",
        attrib={
            "version": "1.1",
            "baseProfile": "full",
            "xmlns:rdkit": "http://www.rdkit.org/xml",
            "xmlns:xlink": "http://www.w3.org/1999/xlink",
            "xml:space": "preserve",
            "width": f"{width}px",
            "height": f"{height}px",
            "viewBox": f"0 0 {width} {height}",
        },
    )
    canvas = ElementTree.Element(
        "rect",
        attrib={
            "style": "opacity:1.0;fill:#FFFFFF;stroke:none",
            "width": f"{width}",
            "height": f"{height}",
            "x": "0.0",
            "y": "0.0",
        },
    )
    arrow = reaction_arrow_svg_element(
        x_start=int(1.1 * image_size),
        x_end=int(1.4 * image_size),
        y_mid=int(0.5 * image_size),
        stroke_width=3,
    )
    root.append(canvas)
    root.append(arrow)

    #   2. Form combined reactant and product graphs
    rgras, rkey_dcts = nonoverlapping_keys_for_sequence(rgras)
    pgras, pkey_dcts = nonoverlapping_keys_for_sequence(pgras)
    rgra = union_from_sequence(rgras)
    pgra = union_from_sequence(pgras)
    rlabel_dct = dict_.merge_sequence(map(dict_.invert, rkey_dcts))
    plabel_dct = dict_.merge_sequence(map(dict_.invert, pkey_dcts))

    if not label and not label_dct:
        rlabel_dct = plabel_dct = None

    if label_dct:
        rlabel_dct = dict_.transform_keys(label_dct, rlabel_dct.get)
        plabel_dct = dict_.transform_keys(label_dct, plabel_dct.get)

    #   3. Generate SVGs for reactants and products
    rsvg = svg_string(
        rgra,
        stereo=stereo,
        exp=exp,
        label=label,
        label_dct=rlabel_dct,
        image_size=image_size,
    )
    psvg = svg_string(
        pgra,
        stereo=stereo,
        exp=exp,
        label=label,
        label_dct=plabel_dct,
        image_size=image_size,
    )

    #   3. Add contents from reactant/product SVGs to combined SVG
    r_group = ElementTree.Element("g")
    p_group = ElementTree.Element("g", attrib={"transform": "translate(300 0)"})
    r_group.extend(ElementTree.fromstring(rsvg).iterfind("path", namespaces=NAMESPACES))
    p_group.extend(ElementTree.fromstring(psvg).iterfind("path", namespaces=NAMESPACES))
    root.append(r_group)
    root.append(p_group)
    return ElementTree.tostring(
        root, encoding="unicode", method="xml", xml_declaration=True
    )


def rdkit_reaction(rgras, pgras, stereo=True, res_stereo=False):
    """Convert reactant and product graphs to an RDKit reaction object.

    This is mainly useful for quick visualization with IPython, which can be
    done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_reaction(pgras, rgras))

        :param rgras: reactant graphs
        :param pgras: product graphs
        :param stereo: Include stereo?
        :type stereo: bool
        :param res_stereo: allow resonant double-bond stereo?
        :type res_stereo: bool
        :param rxn: the reaction object
        :returns: the RDKit reaction
    """
    rdkit_.turn_3d_visualization_off()
    rsmis = [
        smiles(g, stereo=stereo, res_stereo=res_stereo, exp_singles=True) for g in rgras
    ]
    psmis = [
        smiles(g, stereo=stereo, res_stereo=res_stereo, exp_singles=True) for g in pgras
    ]
    rxn_smi = smiles_base.reaction(rsmis, psmis)
    return rdkit_.from_smarts(rxn_smi)


def display(gra, stereo=True, exp=False, label=False, label_dct=None):
    """Display graph to IPython using the RDKit visualizer.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    return ipd.display(
        ipd.SVG(
            svg_string(
                gra=gra, stereo=stereo, exp=exp, label=label, label_dct=label_dct
            )
        )
    )


def display_reaction(rgras, pgras, stereo=True, exp=False, label=False, label_dct=None):
    """Display reaction to IPython using the RDKit visualizer.

    :param rgras: reactant graphs
    :param pgras: product graphs
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param exp: Include explicit hydrogens that aren't needed for stereochemistry?
    :type exp: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    return ipd.display(
        ipd.SVG(
            reaction_svg_string(
                rgras, pgras, stereo=stereo, exp=exp, label=label, label_dct=label_dct
            )
        )
    )


def reaction_arrow_svg_element(
    x_start: int, x_end: int, y_mid: int, stroke_width: int = 3
) -> ElementTree.Element:
    """Generate SVG string for reaction arrow.

    All units are in pixels (px).

    :param x_start: Start x-value
    :param x_end: End x-value
    :param y_mid: Middle y-value
    :param stroke_width: Stroke width, defaults to 3
    :return: SVG string
    """
    x1 = x_start
    x2 = x_end - 3 * stroke_width
    assert x1 < x2, f"{x1} !< {x2}"
    return ElementTree.fromstring(f"""
        <g>
            <defs>
                <marker
                id="arrow"
                viewBox="0 0 12 10"
                refX="9"
                refY="5"
                markerWidth="6"
                markerHeight="6"
                orient="auto-start-reverse">
                    <polygon points="0 0, 12 5, 0 10" />
                </marker>
            </defs>
            <line
                x1="{x1}"
                y1="{y_mid}"
                x2="{x2}"
                y2="{y_mid}"
                stroke="black"
                stroke-width="{stroke_width}px"
                marker-end="url(#arrow)" />
        </g>
    """)


# # TS geometry helpers
def ts_geometry_from_reactants(
    tsg,
    rct_geos,
    geo_idx_dct: Optional[Dict[int, int]] = None,
    check: bool = True,
    log: bool = False,
) -> Any:
    """Generate a TS geometry from reactants.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param rct_geos: Reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :param check: Check stereo and connectivity? defaults to True
    :param log: Print optimization log?, defaults to False
    :return: TS geometry
    :rtype: automol geom data structure
    """
    assert not has_dummy_atoms(tsg), (
        f"Cannot convert graph->geom with dummy atoms:\n{tsg}\n{rct_geos}"
    )

    # 0. Join geometries for bimolecular reactions, yielding a single starting structure
    ts_geo = ts_join_reactant_geometries(tsg, rct_geos, geo_idx_dct=geo_idx_dct)

    # 1. Align the graph and the geometry keys/indices
    tsg, ts_geo, *_, idx_dct = align_with_geometry(tsg, ts_geo, (), geo_idx_dct)

    # 2. Convert to local stereo
    tsg0 = tsg
    tsg = to_local_stereo(tsg)

    # 3. Correct the stereochemistry against the TS graph, so it is consistent with
    # both reactants and products
    ts_geo = stereo_corrected_geometry(tsg, ts_geo, local_stereo=True)

    if log:
        print(f"Raw TS stereo-corrected geometry before cleaning:\n{ts_geo}")

    # 4. Determine the new, aligned keys of the reactants
    key_dct = dict(map(reversed, idx_dct.items()))
    counts = list(map(geom.count, rct_geos))
    rcts_keys = [
        list(map(key_dct.get, range(c0, c0 + c)))
        for c0, c in mit.pairwise([0] + counts)
    ]

    # 5. Embed the TS structure, using distances from the *original* reactant
    # geometries along with forming/breaking bond distances
    ts_geo = clean_geometry(
        tsg,
        ts_geo,
        local_stereo=True,
        none_if_failed=False,
        geos=rct_geos,
        geos_keys=rcts_keys,
        relax_angles=ts.has_reacting_ring(tsg),
        log=log,
    )

    if log:
        print(f"TS geometry after cleaning:\n{ts_geo}")

    if check and not geometry_matches(tsg, ts_geo, local_stereo=True, log=log):
        raise error.FailedGeometryGenerationError(f"\ntsg0 = {tsg0}")

    return ts_geo


def clean_ts_geometry(
    tsg,
    ts_geo,
    geo_idx_dct: Dict[int, int] | None = None,
    check: bool = True,
    log: bool = False,
) -> Any:
    """Clean a TS geometry based on its graph.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param rct_geos: Reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :param check: Check stereo and connectivity? defaults to True
    :param log: Print optimization log?, defaults to False
    :return: TS geometry
    :rtype: automol geom data structure
    """
    assert not has_dummy_atoms(tsg), (
        f"Cannot convert graph->geom with dummy atoms:\n{tsg}"
    )

    # 1. Align the graph and the geometry keys/indices
    tsg, ts_geo, *_ = align_with_geometry(tsg, ts_geo, (), geo_idx_dct)

    # 2. Convert to local stereo
    tsg0 = tsg
    tsg = to_local_stereo(tsg)

    # 3. Correct the stereochemistry against the TS graph, so it is consistent with
    # both reactants and products
    ts_geo = stereo_corrected_geometry(tsg, ts_geo, local_stereo=True)

    if log:
        print(f"Raw TS stereo-corrected geometry before cleaning:\n{ts_geo}")

    # 4. Clean the geometry
    ts_geo = clean_geometry(
        tsg,
        ts_geo,
        local_stereo=True,
        none_if_failed=False,
        relax_angles=ts.has_reacting_ring(tsg),
        log=log,
    )

    if log:
        print(f"TS geometry after cleaning:\n{ts_geo}")

    if check and not geometry_matches(tsg, ts_geo, local_stereo=True, log=log):
        raise error.FailedGeometryGenerationError(f"\ntsg0 = {tsg0}")

    return ts_geo


def ts_join_reactant_geometries(tsg, rct_geos, geo_idx_dct=None):
    """Join two reactant geometries at a single forming bond

    If stereochemistry is present in TS graph, the join will take stereochemistry into
    account.

    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param rct_geos: Reactant geometries
    :type rct_geos: List[automol geom data structure]
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices, defaults to None
    :type geo_idx_dct: dict[int: int], optional
    :param fdist_factor: Set the forming bond distance to this times the average
        van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :returns: The joined geometry
    :rtype: automol geom data structure
    """
    # Early return if there is only one reactant, since there is nothing to join
    if len(rct_geos) == 1:
        return rct_geos[0]
    assert len(rct_geos) == 2, f"Cannot join {len(rct_geos)} reactants, only 2"

    # Combine the geometries into one data structure
    geo = sum(rct_geos, ())

    # Align the graph and the geometry keys/indices
    tsg, geo, *_, idx_dct = align_with_geometry(tsg, geo, (), geo_idx_dct)

    # 0. Identify the keys/indices of atoms to be joined
    #   a. Identify the keys of each fragment
    key_dct = dict(map(reversed, idx_dct.items()))
    len1, len2 = map(geom.count, rct_geos)
    keys1 = tuple(map(key_dct.get, range(len1)))
    keys2 = tuple(map(key_dct.get, range(len1, len1 + len2)))
    #   b. If there are multiple forming keys, join on the one with atom stereochemistry
    par_dct = atom_stereo_parities(tsg)
    frm_keys = ts.forming_bond_keys(tsg)
    frm_key, *_ = sorted(frm_keys, key=lambda bk: sum(par_dct[k] is None for k in bk))
    #   c. Identify which forming key belongs to which fragment
    (fkey1,) = frm_key & set(keys1)
    (fkey2,) = frm_key & set(keys2)

    # 1. Find the atoms on either side of the forming bond
    fxyz1 = geom.coordinates(geo)[fkey1]
    fxyz2 = geom.coordinates(geo)[fkey2]
    # 2. Translate to place them right on top of each other
    geo = geom.translate(geo, numpy.subtract(fxyz1, fxyz2), idxs=keys2)
    fxyz2 = geom.coordinates(geo)[fkey2]
    # 3. Find the reacting electron directions
    rvec1 = geom.ts_reacting_electron_direction(geo, tsg, fkey1)
    rvec2 = geom.ts_reacting_electron_direction(geo, tsg, fkey2)
    # 4. Rotate to align them antiparallel
    rot_ang = vector.angle(rvec1, numpy.negative(rvec2))
    rot_axis = vector.unit_perpendicular(rvec1, rvec2)
    geo = geom.rotate(geo, rot_axis, rot_ang, orig_xyz=fxyz2, idxs=keys2)
    # 5. Translate the second reagent to the appropriate distance away
    fdist = ts.heuristic_bond_distance(tsg, fkey1, fkey2, angstrom=True)
    fvec = numpy.multiply(rvec1, fdist)
    geo = geom.translate(geo, fvec, idxs=keys2, angstrom=True)

    # Restore the original atom ordering of the geometry
    return geom.reorder(geo, idx_dct)


def perturb_geometry_planar_dihedrals(
    gra, geo, geo_idx_dct=None, ang=5.0 * phycon.DEG2RAD, degree=False
):
    """Remove symmetry from a geometry by perturbing planar dihedrals?

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param geo_idx_dct: If they don't already match, specify which graph
        keys correspond to which geometry indices.
    :type geo_idx_dct: dict[int: int]
    :param ang: The angle of rotation
    :type ang: float
    :param degree: Is the angle of rotation in degrees?, default True
    :type degree: bool
    :returns: a geometry in which all planar dihedrals have been shifted by the
        given amount
    :rtype: automol geometry data structure
    """
    if geo is None:
        return None

    ang = ang * phycon.DEG2RAD if degree else ang

    # Align the graph and the geometry keys/indices
    gra, geo, *_, idx_dct = align_with_geometry(gra, geo, (), geo_idx_dct)

    # Keep checking for planar dihedrals and rotating the corresponding bonds
    # until none are left
    cis_keys_lst = geometry_dihedrals_near_value(gra, geo, 0.0, tol=ang)
    trans_keys_lst = geometry_dihedrals_near_value(gra, geo, numpy.pi, tol=ang)

    # Otherwise, adjust the dihedral to give it a non-planar value
    dih_dct = {}
    dih_dct.update({k: 0.0 + ang for k in cis_keys_lst})
    dih_dct.update({k: numpy.pi - ang for k in trans_keys_lst})
    for dih_keys, dih_val in dih_dct.items():
        geo = geom.set_dihedral_angle(geo, dih_keys, dih_val, degree=False, gra=gra)

    # Restore the original atom ordering of the geometry
    return geom.reorder(geo, idx_dct)
