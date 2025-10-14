"""Level 4 geometry functions"""

import functools
import itertools
from collections.abc import Sequence
from typing import Dict, Optional

import numpy
import pyparsing as pp
from numpy.typing import ArrayLike
from pyparsing import pyparsing_common as ppc
from phydat import phycon

from .. import vmat
from ..extern import molfile, py3dmol_, rdkit_
from ..graph import base as graph_base
from ..inchi import base as inchi_base
from ..util import ZmatConv, dict_, heuristic, vector, zmat_conv
from ..zmat import base as zmat_base
from . import _0molsym
from .base import (
    central_angle,
    coordinates,
    count,
    dihedral_angle,
    distance,
    from_data,
    insert,
    is_atom,
    reorder,
    rotate,
    set_coordinates,
    subgeom,
    symbols,
    translate,
    without_dummy_atoms,
    xyz_string,
    xyz_trajectory_string,
)


# # conversions
def graph(geo, stereo=True, local_stereo=False, fix_hyper: bool = True):
    """Generate a molecular graph from the molecular geometry that has
    connectivity information and, if requested, stereochemistry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool, optional
    :param local_stereo: Return local stereo assignments? defaults to False
    :type local_stereo: bool, optional
    :param fix_hyper: Correct hypervalencies by removing the most distant neighbors?
    :rtype: automol molecular graph data structure
    """
    gra = graph_without_stereo(geo, fix_hyper=fix_hyper)
    if stereo:
        gra = graph_base.set_stereo_from_geometry(gra, geo, local_stereo=local_stereo)

    return gra


def graph_without_stereo(geo, dist_factor=None, fix_hyper: bool = True):
    """Generate a molecular graph from the molecular geometry that has
    connectivity information, but not stereochemistry

    Anything less than `dist_factor` times the max of (a.) the sum of covalent radii
    and (b.) the average van der Waals radius between two atoms will be considered
    bonded.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param dist_factor: The multiplier on the distance limit, defaults to None
    :type dist_factor: float, optional
    :param fix_hyper: Correct hypervalencies by removing the most distant neighbors?
    """
    symb_dct = dict(enumerate(symbols(geo)))
    keys = sorted(symb_dct.keys())
    bnd_keys = [
        frozenset({k1, k2})
        for k1, k2 in itertools.combinations(keys, r=2)
        if distance(geo, k1, k2, angstrom=True)
        <= heuristic.bond_distance_limit(
            symb_dct[k1], symb_dct[k2], angstrom=True, dist_factor=dist_factor
        )
    ]

    gra = graph_base.from_data(atm_symb_dct=symb_dct, bnd_keys=bnd_keys)

    # Check for dummy atoms with more than one neighbor
    dummy_keys = graph_base.atom_keys(gra, symb="X")
    nkeys_dct = graph_base.atoms_neighbor_atom_keys(gra)
    dnkeys_dct = dict_.by_key(nkeys_dct, dummy_keys)
    dnkeys_dct = dict_.filter_by_value(dnkeys_dct, lambda n: len(n) > 1)
    if dnkeys_dct:
        for dkey, dnkeys in dnkeys_dct.items():
            # Find the neighbor that is closest to 1 Angstrom away
            best_dnkey = None
            best_diff = numpy.inf
            for dnkey in dnkeys:
                diff = abs(1 - distance(geo, dkey, dnkey, angstrom=True))
                if diff < best_diff:
                    best_diff = diff
                    best_dnkey = dnkey
            # Remove bonds to all but this neighbor
            bad_dnkeys = dnkeys - {best_dnkey}
            bad_bkeys = [(dkey, k) for k in bad_dnkeys]
            gra = graph_base.remove_bonds(gra, bad_bkeys)

    # Remove hypervalencies, if requested
    if fix_hyper:
        nhyp_dct = dict_.by_value(graph_base.atom_hypervalencies(gra))
        nkeys_dct = graph_base.atoms_neighbor_atom_keys(
            graph_base.without_dummy_atoms(gra)
        )

        for key, nhyp in nhyp_dct.items():
            # Get the neighboring keys sorted by distance
            nkeys = sorted(nkeys_dct[key], key=functools.partial(distance, geo, key))
            # Remove the `nhyp` longest hypervalent bonds
            hyp_bkeys = {(key, k) for k in nkeys[-nhyp:]}
            gra = graph_base.remove_bonds(gra, hyp_bkeys)

    return gra


def connectivity_graph_deprecated(
    geo, rqq_bond_max=3.45, rqh_bond_max=2.6, rhh_bond_max=1.9
):
    """Generate a molecular graph from the molecular geometry that has
    information about bond connectivity.

    :param rqq_bond_max: maximum distance between heavy atoms
    :type rqq_bond_max: float
    :param rqh_bond_max: maximum distance between heavy atoms and hydrogens
    :type rqh_bond_max: float
    :param rhh_bond_max: maximum distance between hydrogens
    :type rhh_bond_max: float
    :rtype: automol molecular graph structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    def _distance(idx_pair):
        xyz1, xyz2 = map(xyzs.__getitem__, idx_pair)
        dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
        return dist

    def _are_bonded(idx_pair):
        sym1, sym2 = map(symbs.__getitem__, idx_pair)
        dist = _distance(idx_pair)
        return (
            False
            if "X" in (sym1, sym2)
            else (
                (dist < rqh_bond_max)
                if "H" in (sym1, sym2)
                else (
                    (dist < rhh_bond_max)
                    if (sym1 == "H" and sym2 == "H")
                    else (dist < rqq_bond_max)
                )
            )
        )

    idxs = range(len(xyzs))
    atm_symb_dct = dict(enumerate(symbs))
    bnd_keys = tuple(
        map(frozenset, filter(_are_bonded, itertools.combinations(idxs, r=2)))
    )

    gra = graph_base.from_data(atm_symb_dct=atm_symb_dct, bnd_keys=bnd_keys)

    return gra


def zmatrix(geo, gra=None):
    """Generate a corresponding Z-Matrix for a molecular geometry
    using internal autochem procedures.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A graph identifying the connectivity of this geometry
    :type gra: automol graph data structure
    :returns: automol Z-Matrix data structure
    """
    zma, _ = zmatrix_with_conversion_info(geo, gra=gra)
    return zma


def zmatrix_with_conversion_info(geo, gra=None, zc_: ZmatConv = None):
    """Generate a Z-Matrix for a molecular geometry, along with a z-matrix conversion
    data structure describing the conversion

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A graph identifying the connectivity of this geometry
    :type gra: automol graph data structure
    :param zc_: Request a specific z-matrix conversion, defaults to None
    :type zc_: ZmatConv, optional
    :returns: An automol Z-Matrix data structure and a z-matrix conversion
    :rtype: automol zmat data structure, ZmatConv
    """
    # Handle monatomics separately
    if is_atom(geo):
        symbs = symbols(geo)
        key_mat = [[None, None, None]]
        val_mat = [[None, None, None]]
        zma = zmat_base.from_data(symbs, key_mat, val_mat)
        zc_ = zmat_conv.from_zmat_data(1, {})
        return zma, zc_

    orig_gra = graph_without_stereo(geo) if gra is None else gra

    if zc_ is not None:
        # If a specific z-matrix conversion was requested, apply it to the graph and see
        # if it replicates without reordering
        gra = graph_base.apply_zmatrix_conversion(orig_gra, zc_)
        vma, zkeys = graph_base.vmat.vmatrix(gra)
        assert list(zkeys) == sorted(zkeys), (
            f"Failed to replicate z-matrix conversion:\n{zc_}\n{geo}"
        )
    else:
        # Build an initial z-matrix conversion data structure to put dummies on linear
        # atoms
        nreal = count(geo)
        lin_gkeys = graph_base.linear_atom_keys(orig_gra)
        zc0 = zmat_conv.from_geom_data(nreal, lin_gkeys)

        # Apply this z-matrix conversion to the graph
        gra = graph_base.apply_zmatrix_conversion(orig_gra, zc0)

        # Generate a v-matrix for the graph and get the z-matrix reordering
        vma, zkeys = graph_base.vmat.vmatrix(gra)

        # Re-generate the conversion info, to make sure the direction keys are present
        gkeys = zmat_conv.relabel_zmatrix_key_sequence(zc0, zkeys, dummy=False)
        orig_gkey_dct = dict(enumerate(gkeys))
        zc_ = zmat_conv.relabel(vmat.conversion_info(vma), orig_gkey_dct, "geom")

    # Apply the z-matrix conversion to the geometry and generate the z-matrix
    geo = apply_zmatrix_conversion(geo, zc_, gra=orig_gra)
    zma = zmat_base.from_geometry(vma, geo)

    return zma, zc_


def update_zmatrix(geo, zma, zc_: Optional[ZmatConv] = None):
    """Update a z-matrix from a geometry, optionally specifying the z-matrix conversion

    :param geo: The updated geometry
    :type geo: automol geom data structure
    :param zma: The original z-matrix
    :type zma: automol zmat data structure
    :param zc_: The z-matrix conversion (only needed if geometry is reordered),
        defaults to None
    :type zc_: Optional[ZmatConv], optional
    :returns: The updated z-matrix
    :rtype: automol zmat data structure
    """
    # If no conversion info was passed in, assume no reordering and get it from the
    # z-matrix
    zc_ = zmat_base.conversion_info(zma) if zc_ is None else zc_

    # 1. Apply z-matrix conversion to geometry so it matches
    geo = apply_zmatrix_conversion(geo, zc_)

    # 2. Form the updated z-matrix
    vma = zmat_base.vmatrix(zma)
    zma = zmat_base.from_geometry(vma, geo)

    return zma


def amchi(geo, stereo=True):
    """Generate an AMChI string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    ach, _ = amchi_with_numbers(geo, stereo=stereo)

    return ach


def amchi_with_numbers(geo, stereo=True, gra=None):
    """Determine the AMChI string and sort order for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param gra: molecular graph (to avoid recalculating)
    :type gra: automol molecular graph data structure
    :returns: the AMChI string and AMChI canonical sort ordering for each
        connected component (components in multi-component AMChI ordering)
    :rtype: (str, tuple[tuple[int]])
    """
    gra = graph(geo, stereo=stereo) if gra is None else gra
    ach, num_dcts = graph_base.amchi_with_numbers(gra, stereo=stereo)
    nums_lst = tuple(map(dict_.keys_sorted_by_value, num_dcts))
    return ach, nums_lst


def inchi(geo, stereo=True):
    """Generate an InChI string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    ich, _ = inchi_with_numbers(geo, stereo=stereo)

    return ich


def inchi_with_numbers(geo, stereo=True, gra=None):
    """Determine the InChI string and sort order for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param gra: molecular graph (to avoid recalculating)
    :type gra: automol molecular graph data structure
    :returns: the InChI string and InChI canonical sort ordering for each
        connected component (components in multi-component InChI ordering)
    :rtype: (str, tuple[tuple[int]])
    """
    gra = graph_without_stereo(geo) if gra is None else gra
    if not stereo:
        geo = None
    else:
        gra = graph_base.set_stereo_from_geometry(gra, geo)

    mlf, key_map_inv = molfile_with_atom_mapping(gra, geo=geo)
    rdm = rdkit_.from_molfile(mlf)
    ich, aux_info = rdkit_.to_inchi(rdm, with_aux_info=True)

    nums_lst = _parse_sort_order_from_aux_info(aux_info)
    nums_lst = tuple(tuple(map(key_map_inv.__getitem__, nums)) for nums in nums_lst)

    # Assuming the MolFile InChI works, the above code is all we need. What
    # follows is to correct cases where it fails.
    # This only appears to work sometimes, so when it doesn't, we fall back on
    # the original inchi output.
    if geo is not None:
        gra = graph_base.set_stereo_from_geometry(gra, geo)
        gra = graph_base.implicit(gra)
        sub_ichs = inchi_base.split(ich)

        failed = False

        new_sub_ichs = []
        for sub_ich, nums in zip(sub_ichs, nums_lst):
            sub_gra = graph_base.subgraph(gra, nums, stereo=True)
            sub_ich = _connected_inchi_with_graph_stereo(sub_ich, sub_gra, nums)
            if sub_ich is None:
                failed = True
                break

            new_sub_ichs.append(sub_ich)

        # If it worked, replace the InChI with our forced-stereo InChI.
        if not failed:
            ich = inchi_base.join(new_sub_ichs)
            ich = inchi_base.standard_form(ich)

    return ich, nums_lst


def _connected_inchi_with_graph_stereo(ich, gra, nums):
    """For a connected inchi/graph, check if the inchi is missing stereo; If
    so, add stereo based on the graph.

    Currently only checks for missing bond stereo, since this is all we have
    seen so far, but could be generalized.

    :param ich: the inchi string
    :param gra: the graph
    :param nums: graph indices to backbone atoms in canonical inchi order
    :type nums: tuple[int]
    """
    # First, do a check to see if the InChI is missing bond stereo
    # relative to the graph.
    ich_ste_keys = inchi_base.stereo_bonds(ich)
    our_ste_keys = graph_base.bond_stereo_keys(gra)

    miss_ich_ste_keys = inchi_base.unassigned_stereo_bonds(ich)

    if len(ich_ste_keys) > len(our_ste_keys):
        raise RuntimeError("Our code is missing stereo bonds")

    if len(ich_ste_keys) < len(our_ste_keys) or miss_ich_ste_keys:
        # Convert to implicit graph and relabel based on InChI sort
        atm_key_dct = dict(map(reversed, enumerate(nums)))
        gra = graph_base.relabel(gra, atm_key_dct)
        gra = graph_base.explicit(gra)
        exp_h_keys = graph_base.nonbackbone_hydrogen_keys(gra)
        exp_h_key_dct = {k: -k for k in exp_h_keys}
        gra = graph_base.relabel(gra, exp_h_key_dct)

        gra = graph_base.to_local_stereo(gra)

        # Translate internal stereo parities into InChI stereo parities
        # and generate the appropriate b-layer string for the InChI
        ste_dct = graph_base.bond_stereo_parities(gra)
        ste_keys = tuple(
            sorted(tuple(reversed(sorted(k))) for k in graph_base.bond_stereo_keys(gra))
        )
        blyr_strs = []
        for atm1_key, atm2_key in ste_keys:
            par = ste_dct[frozenset({atm1_key, atm2_key})]

            blyr_strs.append(f"{atm1_key + 1}-{atm2_key + 1}{'+' if par else '-'}")

        # After forming the b-layer string, generate the new InChI
        blyr_str = ",".join(blyr_strs)
        ste_dct = {"b": blyr_str}
        ich = inchi_base.standard_form(ich, ste_dct=ste_dct)

    return ich


def _parse_sort_order_from_aux_info(aux_info):
    PREFIX = pp.Suppress(
        "AuxInfo=" + ppc.integer + "/" + ppc.integer + "/" + pp.Suppress("N:")
    )
    NUMBERS = pp.Group(pp.delimitedList(ppc.integer))
    AUX_INFO = PREFIX + pp.delimitedList(NUMBERS, delim=";")

    nums_lst = tuple(map(tuple, AUX_INFO.parseString(aux_info).asList()))
    return nums_lst


def molfile_with_atom_mapping(
    gra, geo=None, dummy: bool = False, bond_order: bool = True
):
    """Generate an MOLFile from a molecular graph.
    If coordinates are passed in, they are used to determine stereo.

    :param gra: molecular graph
    :type gra: automol graph data structure
    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param dummy: Include dummy atoms?
    :param bond_order: Include bond orders?
    :returns: the MOLFile string, followed by a mapping from MOLFile atoms
        to atoms in the graph
    :rtype: (str, dict)
    """
    gra = graph_base.without_bonds_by_orders(gra, ords=[0], skip_dummies=False)
    if not dummy:
        gra = graph_base.without_dummy_atoms(gra)
        geo = None if geo is None else without_dummy_atoms(geo)

    if bond_order:
        gra = graph_base.kekule(gra)

    atm_keys = sorted(graph_base.atom_keys(gra))
    bnd_keys = list(graph_base.bond_keys(gra))
    atm_syms = dict_.values_by_key(
        graph_base.atom_symbols(gra, dummy_symbol="He"), atm_keys
    )
    atm_bnd_vlcs = dict_.values_by_key(graph_base.atom_bond_counts(gra), atm_keys)
    atm_rad_vlcs = dict_.values_by_key(
        graph_base.atom_unpaired_electrons(gra), atm_keys
    )
    bnd_ords = dict_.values_by_key(graph_base.bond_orders(gra), bnd_keys)

    atm_xyzs = None if geo is None else coordinates(geo, angstrom=True)

    mlf, key_map_inv = molfile.from_data(
        atm_keys,
        bnd_keys,
        atm_syms,
        atm_bnd_vlcs,
        atm_rad_vlcs,
        bnd_ords,
        atm_xyzs=atm_xyzs,
    )
    return mlf, key_map_inv


def chi(geo, stereo=True):
    """Generate a ChI string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    chi_, _ = chi_with_sort(geo, stereo=stereo)

    return chi_


def chi_with_sort(geo, stereo=True, gra=None):
    """Determine the ChI string and sort order for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param gra: molecular graph (to avoid recalculating)
    :type gra: automol molecular graph data structure
    :returns: the AMChI string and AMChI canonical sort ordering for each
        connected component (components in multi-component AMChI ordering)
    :rtype: (str, tuple[tuple[int]])
    """
    gra = graph(geo, stereo=stereo) if gra is None else gra

    # new implementation
    chi_, nums_lst = inchi_with_numbers(geo, stereo=stereo, gra=gra)
    if graph_base.inchi_is_bad(gra, chi_):
        chi_, nums_lst = amchi_with_numbers(geo, stereo=stereo, gra=gra)

    return chi_, nums_lst


def smiles(geo, stereo=True, res_stereo=True):
    """Generate a SMILES string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :param res_stereo: allow resonant double-bond stereo?
    :type res_stereo: bool
    :returns: the SMILES string
    :rtype: str
    """
    gra = graph(geo, stereo=stereo)
    smi = graph_base.smiles(gra, stereo=stereo, res_stereo=res_stereo)
    return smi


def rdkit_molecule(geo, gra=None, stereo=True):
    """Convert a geometry to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(geo))

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_on()
    gra = graph(geo, stereo=stereo) if gra is None else gra
    return rdkit_.from_geometry_with_graph(geo, gra)


def py3dmol_view(
    geo,
    gra=None,
    view=None,
    label: bool = False,
    image_size: int = 400,
    mode: Optional[ArrayLike] = None,
):
    """Get a py3DMol view of this molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity, defaults to None
    :type gra: automol graph data structure, optional
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    :param label: Whether to display atom labels
    :param image_size: The image size, if creating a new view, defaults to 400
    :param mode: A vibrational mode or molecular motion to visualize
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    if gra is not None and graph_base.is_ts_graph(gra):
        gra = graph_base.ts.reactants_graph(gra, stereo=False, dummy=True)

    if mode is not None:
        xyz_str = xyz_string(geo, mode=mode)
        return py3dmol_.view_molecule_from_xyz(
            xyz_str, view=view, image_size=image_size, vib=True
        )

    gra = graph(geo, stereo=False) if gra is None else gra
    mlf_str, key_dct = molfile_with_atom_mapping(gra, geo, dummy=True, bond_order=False)
    view = py3dmol_.view_molecule_from_molfile(
        mlf_str, view=view, image_size=image_size
    )

    # Add labels, if requested
    if label:
        for mlf_key, gra_key in key_dct.items():
            view.addLabel(
                gra_key,
                {
                    "backgroundOpacity": 0.0,
                    "fontColor": "black",
                    "alignment": "center",
                    "inFront": True,
                },
                {"index": mlf_key - 1},
            )

    return view


def display(
    geo,
    gra=None,
    view=None,
    label: bool = False,
    image_size=400,
    vis_bkeys: Optional[tuple[tuple[int, int]]] = None,
    mode: Optional[ArrayLike] = None,
):
    """Display molecule to IPython using the RDKit visualizer.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    :param label: Whether to display atom labels
    :param image_size: The image size, if creating a new view, defaults to 400
    :param vis_bkeys: Only visualize these bonds, by key
    :param mode: A vibrational mode or molecular motion to visualize
    """
    ts_ = gra is not None and graph_base.is_ts_graph(gra)
    if ts_:
        tsg = gra
        gra = graph_base.ts.reactants_graph(gra, stereo=False, dummy=True)

    gra = graph(geo, stereo=False) if gra is None else gra

    # If requested, visualize only a subset of the bonds by removing others
    if vis_bkeys is not None:
        all_bkeys = graph_base.bond_keys(gra)
        vis_bkeys = set(map(frozenset, vis_bkeys))
        gra = graph_base.remove_bonds(gra, all_bkeys - vis_bkeys, stereo=False)
        gra = graph_base.add_bonds(gra, vis_bkeys - all_bkeys)

    view = py3dmol_view(
        geo, gra=gra, view=view, label=label, image_size=image_size, mode=mode
    )

    if ts_:
        for frm_bkey in graph_base.ts.forming_bond_keys(tsg):
            fidx1, fidx2 = frm_bkey
            fxyz1, fxyz2 = coordinates(geo, idxs=(fidx1, fidx2), angstrom=True)
            rvec1 = ts_reacting_electron_direction(geo, tsg, fidx1)
            rvec2 = ts_reacting_electron_direction(geo, tsg, fidx2)
            view = py3dmol_.view_vector(rvec1, orig_xyz=fxyz1, view=view)
            view = py3dmol_.view_vector(rvec2, orig_xyz=fxyz2, view=view)

    return view.show()


def display_trajectory(geos: Sequence[object], image_size=400):
    """Display molecule trajectory to IPython using the RDKit visualizer.

    :param geo: molecular geometry
    :param image_size: The image size, if creating a new view, defaults to 400
    """
    xyz_str = xyz_trajectory_string(geos)
    view = py3dmol_.view_molecule_from_xyz_trajectory(xyz_str, image_size=image_size)
    return view.show()


# # derived properties
def is_connected(geo):
    """Determine if all atoms in geometry are completely connected.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: bool
    """
    comps = graph_base.connected_components(graph(geo, stereo=False))
    return len(comps) == 1


def linear_atoms(geo, gra=None, tol=5.0):
    """find linear atoms in a geometry (atoms with 180 degree bond angle)

    :param geo: the geometry
    :type geo: automol geometry data structure
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :type gra: automol graph data structure
    :param tol: the tolerance threshold for linearity, in degrees
    :type tol: float
    :rtype: tuple(int)
    """

    gra = graph_without_stereo(geo) if gra is None else gra
    ngb_idxs_dct = graph_base.atoms_neighbor_atom_keys(gra)

    lin_idxs = []
    for idx in range(count(geo)):
        nidxs = ngb_idxs_dct[idx]
        if len(nidxs) >= 2:
            for nidx1, nidx2 in itertools.combinations(nidxs, 2):
                ang = central_angle(geo, nidx1, idx, nidx2, degree=True)
                if numpy.abs(ang - 180.0) < tol:
                    lin_idxs.append(idx)

    lin_idxs = tuple(lin_idxs)

    return lin_idxs


def closest_unbonded_atom_distances(
    geo,
    idx,
    gra=None,
    angstrom=False,
    incl_idxs: tuple = (),
    excl_second_degree: bool = True,
    dist_frac: float = 0.2,
) -> Dict[int, float]:
    """For a specific atom in a geometry, find the closest unbonded atoms, along with
    their distances

    Finds the closest unbonded atom and then finds any other atoms within some fraction
    of the minimum distance

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx: The atom to do this for
    :type idx: int
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :type gra: automol graph data structure
    :param angstrom: parameter to control conversion to Angstrom
    :type angstrom: bool
    :param incl_idxs: Specific atoms to include, whether they are bonded or not
    :type incl_idxs: tuple, optional
    :param excl_second_degree: Exclude second-degree neighboring atoms? defaults to True
    :type excl_second_degree: bool, optional
    :param dist_frac: Atoms within this fraction of the minimum distance will be
        included long with the closest atom; defaults to 0.2
    :type dist_frac: float, optional
    :returns: The closest unbonded atom indices along with their distances
    :rtype: Dict[int, float]
    """
    gra = graph_without_stereo(geo) if gra is None else gra
    idxs = graph_base.atom_keys(gra, excl_symbs=["X"])
    idxs -= graph_base.atom_neighbor_atom_keys(
        gra, idx, include_self=True, second_degree=excl_second_degree
    )
    idxs |= set(incl_idxs)

    dist_dct = {i: distance(geo, idx, i, angstrom=angstrom) for i in idxs}
    min_dist = min(dist_dct.values())

    dist_thresh = min_dist * (1.0 + dist_frac)
    dist_dct = dict_.by_value(dist_dct, lambda d: d < dist_thresh)
    return dist_dct


def could_be_forming_bond(geo, idx1: int, idx2: int, gra=None) -> bool:
    """Check whether two atoms could be forming a bond in the geometry

    Checks that each one is among the closest unbonded atoms to the other

    :param geo: A molecular geometry
    :type geo: automol geom data structure
    :param idx1: The first atom index
    :type idx1: int
    :param idx2: The second atom index
    :type idx2: int
    :param gra: A graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :type gra: automol graph data structure
    :return: `True` if it is, `False` if it isn't
    :rtype: bool
    """
    dist_dct1 = closest_unbonded_atom_distances(geo, idx1, gra=gra, incl_idxs=[idx2])
    dist_dct2 = closest_unbonded_atom_distances(geo, idx2, gra=gra, incl_idxs=[idx1])
    return idx1 in dist_dct2 and idx2 in dist_dct1


def ts_reacting_electron_direction(geo, tsg, key) -> vector.Vector:
    """Identify the direction of a reacting electron on a bond-forming atom

    Forming bond direction accounts for atom stereochemistry in this atom.

    :param geo: A geometry aligned to the TS graph
    :type geo: automol geom data structure
    :param tsg: TS graph
    :type tsg: automol graph data structure
    :param key: Key for the atom, which must be part of a forming bond
    :type key: int
    :returns: A vector indicating the direction
    :rtype: vec.Vector
    """
    frm_key = next((k for k in graph_base.ts.forming_bond_keys(tsg) if key in k), None)
    assert frm_key is not None, f"Atom {key} is not forming a bond in this graph:{tsg}"

    # Get the normal vector
    pkeys = graph_base.ts.reacting_atom_plane_keys(tsg, key)
    pxyzs = coordinates(geo, idxs=pkeys)
    zvec = vector.best_unit_perpendicular(pxyzs)

    # Get the direction information:
    #   1. xkey: a bond key giving an x direction
    #   2. ykey: a bond key giving an y direction
    #   3. phi: a rotational angle
    # The electron direction is obtained by rotating the x direction by phi around a z
    # axis of a right-handed coordinate system (happens in the `else` below)
    xkey, ykey, phi = graph_base.ts.reacting_electron_direction(tsg, key)
    # `None` indicates that the electron is perpendicular to the plane, along the normal
    # vector
    if xkey is None:
        rvec = zvec
    # Otherwise, the electron is in-plane, given by rotating an existing bond direction
    # by `phi`
    #   1. do pi rotations by simply reversing the direction of xkey and normalizing
    elif numpy.allclose(phi, numpy.pi):
        xxyz1, xxyz2 = coordinates(geo, idxs=xkey)
        rvec = -vector.unit_norm(numpy.subtract(xxyz2, xxyz1))
    #   2. do other rotations by forming a right-handed coordinate system and rotating
    #   in the x-y plane in the direction of y
    else:
        xxyz1, xxyz2 = coordinates(geo, idxs=xkey)
        xvec = numpy.subtract(xxyz2, xxyz1)
        if ykey is not None:
            yxyz1, yxyz2 = coordinates(geo, idxs=ykey)
            yvec = numpy.subtract(yxyz2, yxyz1)
            zvec = vector.flip_if_left_handed(xvec, yvec, zvec)
        xvec = vector.orthogonalize(zvec, xvec, normalize=True)
        rot_ = vector.rotator(zvec, phi)
        rvec = rot_(xvec)

    # Make sure the direction matches atom stereochemistry
    # Reverse the TS graph before checking stereo, so that Sn2 reactions will be
    # corrected as well (otherwise, it will be checked against the breaking bond, which
    # should already be in place)
    tsg = graph_base.ts.reverse(tsg)
    tsg = graph_base.to_local_stereo(tsg)
    apar_dct = dict_.filter_by_value(
        graph_base.atom_stereo_parities(tsg), lambda x: x is not None
    )
    if key in apar_dct:
        # Create a dummy geometry with the attacking neighbor at this position
        (xyz,) = coordinates(geo, idxs=(key,))
        (key_,) = frm_key - {key}
        xyz_ = numpy.add(xyz, rvec)
        geo_ = set_coordinates(geo, {key_: xyz_})

        # Evaluate the parity of this configuration
        par = graph_base.geometry_atom_parity(tsg, geo_, key)

        # If it doesn't match, reverse the direction, so the atom will be attacked from
        # the other side
        if par != apar_dct[key]:
            rvec = numpy.negative(rvec)

    return rvec


def external_symmetry_factor(geo, chiral_center=True):
    """Obtain the external symmetry factor for a geometry using MolSym

    If requested, divides by the enantiomeric factor

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: float
    """

    if is_atom(geo):
        ext_sym_fac = 1.0
    else:
        pg_obj = _0molsym.point_group_from_geometry(geo)
        ext_sym_fac = _0molsym.point_group_symmetry_number(pg_obj)
        if _0molsym.point_group_is_chiral(pg_obj) and chiral_center:
            ext_sym_fac *= 0.5

    return ext_sym_fac


# # derived operations
def apply_zmatrix_conversion(geo, zc_: ZmatConv, gra=None, dist: float = 1.0):
    """Apply a z-matrix conversion to this geometry, inserting dummy atoms and
    reordering as described in a z-matrix conversion data structure

    :param geo: A molecular geometry
    :type geo: automol geom data structure
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :param gra: A molecular graph, defaults to None
    :type gra: automol graph data structure, optional
    :param dist: The distance of the dummy atom to its parent atom, in Angstroms,
        defaults to 1.0
    :type dist: float, optional
    :returns: The transformed molecular geometry
    :rtype: automol geom data structure
    """
    lin_keys = zmat_conv.linear_atom_keys(zc_, "geom")
    dir_key_dct = dict(zmat_conv.dummy_source_keys(zc_, "geom"))
    gra = graph_without_stereo(geo) if gra is None else gra

    # Partition parent atoms into adjacent segments
    lin_seg_dct = graph_base.linear_segment_cap_keys(gra, lin_keys=lin_keys)
    # lin_seg_keys_lst = graph_base.linear_segments_atom_keys(gra, lin_keys=lin_keys)

    # Build an initial z-matrix conversion data structure to describe the insertion
    nreal = count(geo)
    lin_gkeys = list(itertools.chain(*lin_seg_dct.keys()))
    zc0 = zmat_conv.from_geom_data(nreal, lin_gkeys)

    # Identify a perpendicular direction for each segment and insert the dummy atoms
    # (Parent atom indices don't change, since the dummy atoms are added to the end)
    xyzs = coordinates(geo, angstrom=True)
    for seg_keys, cap_keys in lin_seg_dct.items():
        key1, key2 = cap_keys
        if key1 == key2:
            key1, *_ = seg_keys
            key2 = dir_key_dct[key1]
        xyz1, xyz2 = coordinates(geo, idxs=(key1, key2))
        dir_vec = vector.unit_norm(numpy.subtract(xyz2, xyz1))
        dum_vec = vector.arbitrary_unit_perpendicular(dir_vec)
        for idx in seg_keys:
            xyz = numpy.add(xyzs[idx], numpy.multiply(dist, dum_vec))
            geo = insert(geo, "X", xyz, angstrom=True)

    # Reorder the geometry to match the desired conversion
    geo = reorder(geo, zmat_conv.isomorphism(zc0, zc_))
    return geo


def undo_zmatrix_conversion(geo, zc_: ZmatConv):
    """Undo a z-matrix conversion of this geometry, removing dummy atoms and
    reversing any reordering

    :param geo: A molecular geometry
    :type geo: automol geom data structure
    :param zc_: A z-matrix conversion
    :type zc_: ZmatConv
    :returns: The transformed molecular geometry
    :rtype: automol geom data structure
    """
    rel_dct = zmat_conv.relabel_dict(zc_, "zmat")
    idxs = zmat_conv.real_keys(zc_, "zmat")
    idxs = sorted(idxs, key=rel_dct.__getitem__)
    return subgeom(geo, idxs)


def set_distance(
    geo,
    dist_idxs,
    dist_val,
    angstrom=True,
    gra=None,
):
    """Set a particular bond distance in the geometry to a new value

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param dist_idxs: A pair of indices, identifying the bond distance to be adjusted;
        the second in the pair is the atom that will be moved, along with atoms
        connected to it
    :type dist_idxs: Tuple[int, int]
    :param dist_val: The desired distance coordinate value
    :type dist_val: float
    :param angstrom: are distances in angstrom? If not, assume bohr.
    :type angstrom: bool
    :param gra: molecular graph for tracking connectivity (will be generated if None)
    :type gra: automol molecular graph data structure
    :rtype: automol geometry data structure
    """
    assert len(dist_idxs) == 2
    idx1, idx2 = dist_idxs
    gra = gra if gra is not None else graph_without_stereo(geo)

    # For TS graphs, cut the branch of a reacting ring at the appropriate place
    # (It would be better to have a more systematic approach to this to handle all rings
    # in the graph...)
    ts_rng_keys = graph_base.vmat.ts_zmatrix_starting_ring_keys(gra)
    if ts_rng_keys is not None:
        drop_bkey = (ts_rng_keys[0], ts_rng_keys[-1])
        gra = graph_base.remove_bonds(gra, [drop_bkey])

    dist_val = dist_val if not angstrom else dist_val * phycon.ANG2BOHR
    idxs = graph_base.branch_atom_keys(gra, idx1, idx2)

    xyzs = coordinates(geo)
    bvec0 = numpy.subtract(xyzs[idx2], xyzs[idx1])
    bvec = dist_val * bvec0 / numpy.linalg.norm(bvec0)
    disp = bvec - bvec0
    geo = translate(geo, disp, idxs=idxs)

    return geo


def set_central_angle(
    geo,
    ang_idxs,
    ang_val,
    degree=True,
    plane_idx=None,
    gra=None,
):
    """Set a particular bond angle in the geometry to a new value

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ang_idxs: A triplet of indices, identifying the central angle to be adjusted;
        the last in the triplet is the atom that will be moved, along with atoms
        connected to it
    :type ang_idxs: Tuple[int, int, int]
    :param ang_val: The desired angle coordinate value
    :type ang_val: float
    :param degree: are angles in degrees? If not, assume radians.
    :type degree: bool
    :param plane_idx: An additional atom that will be used to define the plane of
        rotation, in the event that the angle is initially linear
    :param gra: molecular graph for tracking connectivity (will be generated if None)
    :type gra: automol molecular graph data structure
    :rtype: automol geometry data structure
    """
    idx1, idx2, idx3 = ang_idxs
    gra = gra if gra is not None else graph_without_stereo(geo)
    ang_val = ang_val if not degree else ang_val * phycon.DEG2RAD
    idxs = graph_base.branch_atom_keys(gra, idx2, idx3)
    xyzs = coordinates(geo)
    ang0 = central_angle(geo, idx3, idx2, idx1)

    tol = 2.0 * phycon.DEG2RAD
    lin = numpy.pi

    ang_diff = ang_val - ang0
    # If idx-idx1-idx2 is not linear, use the normal to this plane,
    # otherwise, use the normal to the idx1-idx2-idx3 plane
    if not numpy.abs(ang0) < tol or numpy.abs(ang0 - lin) < tol:
        axis = vector.unit_perpendicular(xyzs[idx3], xyzs[idx1], orig_xyz=xyzs[idx2])
    else:
        axis = vector.unit_perpendicular(
            xyzs[idx2], xyzs[plane_idx], orig_xyz=xyzs[idx2]
        )
    # I don't know how to figure out which way to rotate, so just try both
    # and see which one works
    for dang in [-ang_diff, ang_diff]:
        geo_ = rotate(geo, axis, dang, orig_xyz=xyzs[idx2], idxs=idxs)
        ang_out = central_angle(geo_, idx3, idx2, idx1)
        abs_diff = numpy.abs(
            numpy.mod(ang_val, 2 * numpy.pi) - numpy.mod(ang_out, 2 * numpy.pi)
        )
        ang_comp = numpy.abs(numpy.pi - numpy.abs(abs_diff - numpy.pi))
        if ang_comp < tol:
            geo = geo_
            break

    return geo


def set_dihedral_angle(
    geo,
    dih_idxs,
    dih_val,
    degree=True,
    gra=None,
):
    """Change the z-matrix coordinates of a given atom, shifting those
    connected to it accordingly.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param dih_idxs: A quartet of indices, identifying the dihedral angle to be
        adjusted; the last in the quarted is the atom that will be moved, along with
        atoms connected to it
    :type dih_idxs: Tuple[int, int, int]
    :param dih_val: The desired angle coordinate value
    :type dih_val: float
    :param degree: are angles in degrees? If not, assume radians.
    :type degree: bool
    :param gra: molecular graph for tracking connectivity (will be generated if None)
    :type gra: automol molecular graph data structure
    :rtype: automol geometry data structure
    """
    idx1, idx2, idx3, idx4 = dih_idxs
    gra = gra if gra is not None else graph_without_stereo(geo)
    dih_val = dih_val if not degree else dih_val * phycon.DEG2RAD
    idxs = graph_base.branch_atom_keys(gra, idx3, idx4)
    xyzs = coordinates(geo)
    dih0 = dihedral_angle(geo, idx4, idx3, idx2, idx1)

    ddih = dih_val - dih0
    axis = numpy.subtract(xyzs[idx3], xyzs[idx2])
    geo = rotate(geo, axis, ddih, orig_xyz=xyzs[idx3], idxs=idxs)

    return geo


# interfaces
def ase_atoms(geo):
    """Get an ASE Atoms object from a molecular geometry

    :param geo: A molecular geometry
    :return: The ASE Atoms object
    """
    from ase import Atoms

    return Atoms(symbols=symbols(geo), positions=coordinates(geo, angstrom=True))


def from_ase_atoms(atms_obj):
    """Read a molecular geometry from an ASE Atoms object

    :param atms_obj: An ASE Atoms object
    :return: The molecular geometry
    """
    symbs = atms_obj.get_chemical_symbols()
    xyzs = atms_obj.get_positions()
    return from_data(symbs, xyzs, angstrom=True)
