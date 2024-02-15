""" graph conversions
"""
import functools
import itertools

import IPython.display as ipd
import ipywidgets
import numpy
from phydat import phycon
from automol import error, geom
from automol.extern import rdkit_
from automol.graph._0embed import clean_geometry, embed_geometry, geometry_matches
from automol.graph.base import (
    amchi,
    angle_keys,
    atom_keys,
    atom_stereo_parities,
    bond_keys,
    connected_components,
    explicit,
    geometry_correct_linear_vinyls,
    geometry_dihedrals_near_value,
    has_stereo,
    inchi_is_bad,
    is_ts_graph,
    relabel,
    smiles,
    standard_keys,
    stereo_corrected_geometry,
    string,
    to_local_stereo,
    ts,
    with_explicit_stereo_hydrogens,
    without_stereo,
)
from automol.smiles import base as smiles_base
from automol.util import vector


# # conversions
def geometry(gra, fdist_factor=1.1, bdist_factor=0.9, check=True, log=False):
    """Convert a molecular graph to a molecular geometry.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param fdist_factor: For a TS graph, set the forming bond distance to this times
        the average van der Waals radius, defaults to 1.1
    :type fdist_factor: float, optional
    :param bdist_factor: For a TS graph, set the breaking bond distance to this times
        the average van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param check: Check stereo and connectivity? defaults to True
    :type check: bool, optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
    :rtype: automol molecular geometry data structure
    """
    is_ts = is_ts_graph(gra)
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
            tsg,
            geos,
            geo_idx_dct=geo_idx_dct,
            fdist_factor=fdist_factor,
            bdist_factor=bdist_factor,
            check=check,
            log=log,
        )
    else:
        geo = geom.reorder(geo, gra_key_dct)

    return geo


def _connected_geometry(gra, check=True, log=False):
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

    # Define geometry generation methods
    def method1_(gra_):
        rdm = rdkit_.from_graph(gra_, stereo=stereo, local_stereo=True)
        geo = rdkit_.to_geometry(rdm)
        return geo

    def method2_(gra_):
        rdm = rdkit_.from_graph(gra_, stereo=stereo, local_stereo=True)
        (geo,) = rdkit_.to_conformers(rdm, nconfs=1)
        return geo

    def method3_(gra_):
        return embed_geometry(gra_)

    # Try geometry generation methods until one works
    methods_ = [method1_, method2_, method2_]
    methods_ += [method3_] if not stereo else []
    for try_number, method_ in enumerate(methods_):
        geo = method_(gra)

        if log:
            print(f"Try {try_number}...")
            print("Raw geometry:")
            print(geom.round_(geo))

        geo = _clean_and_validate_connected_geometry(
            gra, geo, stereo=stereo, check=check, log=log
        )

        if geo is not None:
            return geo

    raise error.FailedGeometryGenerationError(f"Failed graph:\n{string(orig_gra)}")


def _clean_and_validate_connected_geometry(
    gra, geo, stereo=True, local_stereo=True, check=True, log=False
):
    """Validate and clen up a connected geometry

    :param gra: Connected molecular graph with standard keys
    :type gra: automol graph data structure
    :param geo: Molecular geometry
    :type geo: automol geom data structure
    :param stereo: Take stereochemistry into consideration? defaults to True
    :type stereo: bool, optional
    :param local_stereo: Does the graph have local stereo assignments? defaults to True
    :type local_stereo: bool, optional
    :param check: Check stereo and connectivity? defaults to True
    :type check: bool, optional
    :param log: Log information to the screen? defaults to False
    :type log: bool, optional
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


def inchi(gra, stereo=True):
    """Generate an InChI string from a molecular graph.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    smi = smiles(gra, stereo=stereo, res_stereo=False)
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


def rdkit_molecule(gra, stereo=True, label=False, label_dct=None):
    """Convert a molecular graph to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(gra))

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    :returns: the RDKit molecule
    """
    if stereo:
        assert gra == with_explicit_stereo_hydrogens(gra), (
            f"Graph is missing hydrogens needed to depict stereo:\n{gra}\n"
            f"You can add them using graph.with_explicit_stereo_hydrogens()"
        )

    rdkit_.turn_3d_visualization_off()
    rdm = rdkit_.from_graph(gra, stereo=stereo, label=label, label_dct=label_dct)
    return rdm


def svg_string(gra, stereo=True, label=False, label_dct=None, image_size=200):
    """Get an SVG string for visualizing the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
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
    rdm = rdkit_molecule(gra, stereo=stereo, label=label, label_dct=label_dct)
    svg_str = rdkit_.to_svg_string(rdm, image_size=image_size)
    return svg_str


def ipywidget(gra, stereo=True, label=False, label_dct=None, image_size=300):
    """Get an ipywidget object for visualizing the graph

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    :param image_size: The image size, defaults to 150
    :type image_size: int, optional
    :return: The widget object
    :rtype: ipywidgets.Image
    """
    svg_str = svg_string(
        gra, stereo=stereo, label=label, label_dct=label_dct, image_size=image_size
    )
    widget = ipywidgets.Image(
        value=svg_str.encode("utf-8"),
        format="svg+xml",
        width=image_size,
        height=image_size,
    )
    return widget


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


def display(gra, stereo=True, label=False, label_dct=None):
    """Display graph to IPython using the RDKit visualizer

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param stereo: Include stereochemistry information?
    :type stereo: bool
    :param label: Display the molecule with atom labels?
    :type label: bool
    :param label_dct: Atom labels, by atom key.  If `None` and `label` is
        `True`, the atom keys themselves will be used as labels.
    :param label_dct: bool
    """
    arrow_svg_str = """
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 100">
    <defs>
        <marker id="arrowhead" markerWidth="5" markerHeight="5" 
        refX="0" refY="2.5" orient="auto">
        <polygon points="0 0, 5 2.5, 0 5" />
        </marker>
    </defs>
    <line x1="20" y1="50" x2="120" y2="50" stroke="#000" 
    stroke-width="10" marker-end="url(#arrowhead)" />
    </svg>
    """

    rdkit_.turn_3d_visualization_off()
    if is_ts_graph(gra):
        rgra = ts.reactants_graph(gra, stereo=stereo)
        pgra = ts.products_graph(gra, stereo=stereo)
        rwid = ipywidget(rgra, stereo=stereo, label=label, label_dct=label_dct)
        pwid = ipywidget(pgra, stereo=stereo, label=label, label_dct=label_dct)
        arrow_widget = ipywidgets.Image(
            value=arrow_svg_str.encode("utf-8"), format="svg+xml", width=100, height=50
        )
        ipd.display(ipywidgets.HBox([rwid, arrow_widget, pwid]))
    else:
        ipd.display(
            rdkit_molecule(gra, stereo=stereo, label=label, label_dct=label_dct)
        )


def display_reaction(rgras, pgras, stereo=True):
    """Display reaction to IPython using the RDKit visualizer

    :param rgras: reactant graphs
    :param pgras: product graphs
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    """
    rdkit_.turn_3d_visualization_off()
    return ipd.display(rdkit_reaction(rgras, pgras, stereo=stereo))


# # TS geometry helpers
def ts_geometry_from_reactants(
    tsg,
    rct_geos,
    geo_idx_dct=None,
    fdist_factor=1.1,
    bdist_factor=0.9,
    max_dist_err=0.2,
    check=True,
    log=False,
):
    """Generate a TS geometry from reactants

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
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param max_dist_err: The distance convergence threshold, in angstroms
    :type max_dist_err: float, optional
    :param check: Check stereo and connectivity? defaults to True
    :type check: bool, optional
    :param log: Print optimization log?, defaults to False
    :type log: bool, optional
    :return: TS geometry
    :rtype: automol geom data structure
    """
    keys = sorted(atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )
    tsg = relabel(tsg, geo_idx_dct)

    # 0. Determine appropriate distance ranges
    dist_range_dct = ts_distance_ranges_from_reactant_geometries(
        tsg,
        rct_geos,
        fdist_factor=fdist_factor,
        bdist_factor=bdist_factor,
        angstrom=True,
    )

    # 1. Join geometries for bimolecular reactions, yielding a single starting structure
    if len(rct_geos) > 1:
        ts_geo = ts_join_reactant_geometries(
            tsg,
            rct_geos,
            fdist_factor=fdist_factor,
        )
    else:
        (ts_geo,) = rct_geos

    # 2. Convert to local stereo. (Must be done after joining, because join function
    # reverses the TS graph before localizing to handle Sn2 reactions correctly)
    tsg = to_local_stereo(tsg)

    # 3. Correct the stereochemistry against the TS graph, so it is consistent with
    # both reactants and products
    ts_geo = stereo_corrected_geometry(tsg, ts_geo, local_stereo=True)

    if log:
        print(f"Raw TS stere-corrected geometry before cleaning:\n{ts_geo}")

    # 4. Embed the TS structure, using distances from the *original* reactant
    # geometries along with forming/breaking bond distances
    ts_geo = clean_geometry(
        tsg,
        ts_geo,
        rct_geos=rct_geos,
        dist_range_dct=dist_range_dct,
        max_dist_err=max_dist_err,
        relax_angles=ts.has_reacting_ring(tsg),
        local_stereo=True,
        none_if_failed=False,
        log=log,
    )

    if log:
        print(f"TS geometry after cleaning:\n{ts_geo}")

    if check and not geometry_matches(tsg, ts_geo, local_stereo=True, log=log):
        raise error.FailedGeometryGenerationError(f"Failed TS graph:\n{string(tsg)}")

    # 6. Re-order the TS geometry to match the TS graph
    gra_key_dct = dict(map(reversed, geo_idx_dct.items()))
    ts_geo = geom.reorder(ts_geo, gra_key_dct)
    return ts_geo


def ts_join_reactant_geometries(tsg, rct_geos, geo_idx_dct=None, fdist_factor=1.1):
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
    akeys = sorted(atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(akeys)} if geo_idx_dct is None else geo_idx_dct
    )

    tsg = relabel(tsg, geo_idx_dct)
    # If there are multiple forming keys, join on the one with atom stereochemistry
    par_dct = atom_stereo_parities(tsg)
    frm_key, *_ = sorted(
        ts.forming_bond_keys(tsg),
        key=lambda bk: sum(par_dct[k] is None for k in bk),
    )

    assert (
        len(rct_geos) == 2
    ), f"This requires two reactants, but {len(rct_geos)} were given"
    len1, len2 = map(geom.count, rct_geos)
    idxs1 = tuple(range(len1))
    idxs2 = tuple(range(len1, len1 + len2))
    (fidx1,) = frm_key & set(idxs1)
    (fidx2,) = frm_key & set(idxs2)
    geo = sum(rct_geos, ())
    # 1. Find the atoms on either side of the forming bond
    fxyz1 = geom.coordinates(geo)[fidx1]
    fxyz2 = geom.coordinates(geo)[fidx2]
    # 2. Translate to place them right on top of each other
    geo = geom.translate(geo, numpy.subtract(fxyz1, fxyz2), idxs=idxs2)
    fxyz2 = geom.coordinates(geo)[fidx2]
    # 3. Find the reacting electron directions
    rvec1 = geom.ts_reacting_electron_direction(geo, tsg, fidx1)
    rvec2 = geom.ts_reacting_electron_direction(geo, tsg, fidx2)
    # 4. Rotate to align them antiparallel
    rot_ang = vector.angle(rvec1, numpy.negative(rvec2))
    rot_axis = vector.unit_perpendicular(rvec1, rvec2)
    geo = geom.rotate(geo, rot_axis, rot_ang, orig_xyz=fxyz2, idxs=idxs2)
    # 5. Translate the second reagent to the appropriate distance away
    fdist = ts.heuristic_bond_distance(
        tsg, fidx1, fidx2, fdist_factor=fdist_factor, angstrom=True
    )
    fvec = numpy.multiply(rvec1, fdist)
    geo = geom.translate(geo, fvec, idxs=idxs2, angstrom=True)

    return geo


def ts_distance_ranges_from_reactant_geometries(
    tsg,
    rct_geos,
    geo_idx_dct=None,
    fdist_factor=1.1,
    bdist_factor=0.9,
    angles=True,
    angstrom=True,
):
    """return a dictionary of distances for certain atoms in a sequence of
    reactant geometries, shifting the keys as needed

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
    :param bdist_factor: Set the breaking bond distance to this times the average
        van der Waals radius, defaults to 0.9
    :type bdist_factor: float, optional
    :param angles: include distances between the ends of bond angles?
    :type angles: bool
    :param angstrom: return the distances in angstroms?
    :type angstrom: bool
    """
    keys = sorted(atom_keys(tsg))
    geo_idx_dct = (
        {k: i for i, k in enumerate(keys)} if geo_idx_dct is None else geo_idx_dct
    )

    geo = sum(rct_geos, ())
    tsg = relabel(tsg, geo_idx_dct)
    gra = ts.reactants_graph(tsg)

    # Add measured bond distances from the reactants geometry
    bnd_keys = bond_keys(gra)
    dists = [geom.distance(geo, *k, angstrom=angstrom) for k in bnd_keys]
    dist_dct = dict(zip(bnd_keys, dists))
    # Add measured angle distances from the reactants geometry, if requested
    if angles:
        ang_keys = [frozenset({k1, k3}) for k1, _, k3 in angle_keys(gra)]
        dists = [geom.distance(geo, *k, angstrom=angstrom) for k in ang_keys]
        dist_dct = dict(zip(ang_keys, dists))

    # Add heuristic bond distances from the TS graph
    for key in ts.reacting_bond_keys(tsg):
        dist_dct[key] = ts.heuristic_bond_distance(
            tsg,
            *key,
            fdist_factor=fdist_factor,
            bdist_factor=bdist_factor,
            angstrom=angstrom,
        )

    dist_range_dct = {k: (d, d) for k, d in dist_dct.items()}
    return dist_range_dct


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
    ang = ang * phycon.DEG2RAD if degree else ang
    geo_idx_dct = (
        geo_idx_dct
        if geo_idx_dct is not None
        else {k: i for i, k in enumerate(sorted(atom_keys(gra)))}
    )

    # Keep checking for planar dihedrals and rotating the corresponding bonds
    # until none are left
    cis_keys_lst = geometry_dihedrals_near_value(
        gra, geo, 0.0, geo_idx_dct=geo_idx_dct, tol=ang
    )
    trans_keys_lst = geometry_dihedrals_near_value(
        gra, geo, numpy.pi, geo_idx_dct=geo_idx_dct, tol=ang
    )

    # Otherwise, adjust the dihedral to give it a non-planar value
    dih_dct = {}
    dih_dct.update({k: 0.0 + ang for k in cis_keys_lst})
    dih_dct.update({k: numpy.pi - ang for k in trans_keys_lst})
    for dih_keys, dih_val in dih_dct.items():
        dih_idxs = list(map(geo_idx_dct.__getitem__, dih_keys))
        geo = geom.set_dihedral_angle(
            geo,
            dih_idxs,
            dih_val,
            degree=False,
            gra=gra,
        )
    return geo
