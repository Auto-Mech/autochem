""" graph conversions
"""
import functools
import itertools

import IPython.display as ipd
import ipywidgets
import numpy

import automol.inchi.base
import automol.smiles.base
from automol import error, geom
from automol.extern import rdkit_
from automol.graph._0embed import (
    clean_geometry,
    distance_ranges_from_coordinates,
    embed_geometry,
)
from automol.graph.base import (
    amchi,
    angle_keys,
    atom_keys,
    atom_stereo_parities,
    bond_keys,
    connected_components,
    explicit,
    has_stereo,
    inchi_is_bad,
    is_ts_graph,
    isomorphism,
    linear_vinyl_corrected_geometry,
    relabel,
    set_stereo_parities,
    smiles,
    standard_keys,
    stereo_corrected_geometry,
    stereo_keys,
    stereo_parities,
    string,
    ts,
)
from automol.util import vec


# # conversions
def geometry(gra, check=True):
    """Convert a molecular graph to a molecular geometry.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol molecular geometry data structure
    """
    is_ts = is_ts_graph(gra)
    if is_ts:
        tsg = gra
        gra = ts.reactants_graph(tsg)

    gra = explicit(gra)
    gras = connected_components(gra)

    geos = [_connected_geometry(g, check=check) for g in gras]
    geos = [geom.translate(g, [50.0 * i, 0.0, 0.0]) for i, g in enumerate(geos)]
    geo = functools.reduce(geom.join, geos)

    # Work out how the geometry needs to be re-ordered to match the input graph
    all_keys = itertools.chain(*(sorted(atom_keys(g)) for g in gras))
    gra_key_dct = dict(enumerate(all_keys))

    if is_ts:
        geo_idx_dct = dict(map(reversed, gra_key_dct.items()))
        geo = ts_geometry_from_reactants(tsg, geos, geo_idx_dct=geo_idx_dct)
    else:
        geo = geom.reorder(geo, gra_key_dct)

    return geo


def _connected_geometry(gra, check=True):
    """Generate a geometry for a connected molecular graph.
    :param gra: connected molecular graph
    :type gra: automol graph data structure
    :param check: check stereo and connectivity?
    :type check: bool
    :rtype: automol molecular geometry
    """
    # Standardize the graph before doing anything else
    gra = standard_keys(gra)

    ste_keys = stereo_keys(gra)

    smi = smiles(gra, res_stereo=False)
    has_ste = has_stereo(gra)

    def _gen1():
        nonlocal smi

        rdm = rdkit_.from_smiles(smi)
        (geo,) = rdkit_.to_conformers(rdm, nconfs=1)
        return geo

    def _gen3():
        nonlocal gra

        if has_ste:
            raise ValueError

        gra_ = explicit(gra)
        geo = embed_geometry(gra_)
        return geo

    success = False
    for gen_ in (_gen1, _gen1, _gen1, _gen3):
        try:
            geo = gen_()
        except (RuntimeError, TypeError, ValueError):
            continue

        if check:
            # First, check connectivity.
            gra_ = geom.graph(geo)
            geo = linear_vinyl_corrected_geometry(gra_, geo)
            geo = clean_geometry(gra_, geo, stereo=False)
            gra_ = geom.graph(geo)

            idx_dct = isomorphism(gra_, gra, stereo=False)

            if idx_dct is None:
                continue

            # Reorder the geometry to match the input graph connectivity.
            geo = geom.reorder(geo, idx_dct)

            # If connectivity matches and there is no stereo, we are done.
            if not has_ste:
                success = True
                break

            # Otherwise, there is stereo.
            # First, try an isomorphism to see if the parities already match.
            gra_ = geom.graph(geo)
            par_dct_ = stereo_parities(gra_)
            par_dct_ = {k: (p if k in ste_keys else None) for k, p in par_dct_.items()}
            gra_ = set_stereo_parities(gra_, par_dct_)
            idx_dct = isomorphism(gra_, gra, stereo=True)

            # If connectivity and stereo match, we are done.
            if idx_dct:
                # Reorder the geometry to match the input graph
                geo = geom.reorder(geo, idx_dct)
                success = True
                break

            # If the stereo doesn't match, try a stereo correction.
            print("A", geo)
            geo = stereo_corrected_geometry(gra, geo, local_stereo=False)
            print("B", geo)
            geo = clean_geometry(gra_, geo, stereo=True)
            print("C", geo)

            # Now, re-try the isomorphism
            gra_ = geom.graph(geo)
            par_dct_ = stereo_parities(gra_)
            par_dct_ = {k: (p if k in ste_keys else None) for k, p in par_dct_.items()}
            gra_ = set_stereo_parities(gra_, par_dct_)
            idx_dct = isomorphism(gra_, gra, stereo=True)

            # If this fails, this geometry won't work. Continue
            if not idx_dct:
                continue

            # The stereo matches after correction.
            # Reorder the geometry to match the input graph, and we are done.
            geo = geom.reorder(geo, idx_dct)
            success = True
            break

    if not success:
        raise error.FailedGeometryGenerationError(f"Failed graph:\n{string(gra)}")

    return geo


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
    if stereo and (label or label_dct is not None):
        raise NotImplementedError("Stereo with display labels not implemented.")

    rdkit_.turn_3d_visualization_off()
    if stereo:
        rdm = rdkit_.from_smiles(smiles(gra, stereo=True, res_stereo=False))
    else:
        rdm = rdkit_.from_graph(gra, stereo=False, label=label, label_dct=label_dct)

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
    rxn_smi = automol.smiles.base.reaction(rsmis, psmis)
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
    rdkit_.turn_3d_visualization_off()
    ipd.display(rdkit_molecule(gra, stereo=stereo, label=label, label_dct=label_dct))


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
    max_dist_err=2e-1,
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

    # 1. Join geometries for bimolecular reactions, yielding a single starting structure
    if len(rct_geos) > 1:
        ts_geo = ts_join_reactant_geometries(
            tsg,
            rct_geos,
            fdist_factor=fdist_factor,
        )
    else:
        (ts_geo,) = rct_geos

    # 2. Correct the stereochemistry against the TS graph, so it is consistent with both
    # reactants and products
    ts_geo = stereo_corrected_geometry(tsg, ts_geo)

    # 3. Embed the TS structure, using distances from the *original* reactant geometries
    # along with forming/breaking bond distances
    dist_dct = ts_distances_from_reactant_geometries(
        tsg,
        rct_geos,
        fdist_factor=fdist_factor,
        bdist_factor=bdist_factor,
    )
    dist_range_dct = distance_ranges_from_coordinates(tsg, dist_dct, angstrom=True)
    ts_geo = clean_geometry(
        tsg,
        ts_geo,
        rct_geos=rct_geos,
        dist_range_dct=dist_range_dct,
        relax_angles=ts.has_reacting_ring(tsg),
        max_dist_err=max_dist_err,
        log=log,
    )

    # 4. Re-order the TS geometry to match the TS graph
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
    rot_ang = vec.angle(rvec1, numpy.negative(rvec2))
    rot_axis = vec.unit_perpendicular(rvec1, rvec2)
    geo = geom.rotate(geo, rot_axis, rot_ang, orig_xyz=fxyz2, idxs=idxs2)
    # 5. Translate the second reagent to the appropriate distance away
    fdist = ts.heuristic_bond_distance(
        tsg, fidx1, fidx2, fdist_factor=fdist_factor, angstrom=True
    )
    fvec = numpy.multiply(rvec1, fdist)
    geo = geom.translate(geo, fvec, idxs=idxs2, angstrom=True)

    return geo


def ts_distances_from_reactant_geometries(
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

    return dist_dct
