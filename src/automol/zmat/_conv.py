""" Level 4 Z-Matrix functions
"""

import itertools

import numpy

from .. import geom, util
from ..graph import base as graph_base
from ..util import ZmatConv
from .base import (
    conversion_info,
    dihedral_angle_coordinates,
    distance_coordinates,
    dummy_keys,
    dummy_source_dict,
    key_matrix,
    name_matrix,
    neighbor_keys,
    string,
    symbols,
    value_matrix,
)
from .base import (
    keys as zmatrix_keys,
)


# # conversions
def graph(zma, stereo: bool = True, dummy: bool = False, zmat_bonds: bool = False):
    """Convert a Z-Matrix to a molecular graph.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param stereo: Calculate stereochemistry?
    :param dummy: parameter to include dummy atoms
    :param zmat_bonds: Include only bonds that are in the z-matrix?
    :rtype: (automol molecular geometry data structure, dict[int, int])
    """
    geo = geometry(zma, dummy=True)
    if not dummy:
        geo = geom.without_dummy_atoms(geo)

    gra = geom.graph(geo, stereo=stereo)

    if zmat_bonds:
        all_bkeys = graph_base.bond_keys(gra)
        zma_bkeys = set(map(frozenset, distance_coordinates(zma).values()))
        gra = graph_base.remove_bonds(gra, all_bkeys - zma_bkeys, stereo=stereo)
        gra = graph_base.add_bonds(gra, zma_bkeys - all_bkeys)

    return gra


def graph_without_stereo(zma, dummy=False, dist_factor=None):
    """Convert the Z-Matrix to a molecular graph that has connectivity information, but
    not stereochemistry

    Anything less than `dist_factor` times the max of (a.) the sum of covalent radii
    and (b.) the average van der Waals radius between two atoms will be considered
    bonded.

    :param zma: Z-Matrix
    :type zma: automol zmat data structure
    :param dist_factor: The multiplier on the distance limit, defaults to None
    :type dist_factor: float, optional
    """

    geo = geometry(zma, dummy=dummy)
    gra = geom.graph_without_stereo(geo, dist_factor=dist_factor)

    return gra


def geometry(zma, dummy=False, zc_: ZmatConv = None):
    """Convert a Z-Matrix to a molecular geometry that includes dummy atoms.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param dummy: include dummy atoms in the geometry?
    :type dummy: bool
    :param zc_: Restore the original geometry before conversion by reversing the
        corresponding z-matrix conversion; defaults to None
    :type zc_: ZmatConv, optional
    :returns: automol molecular geometry data structure
    """

    syms = symbols(zma)

    natms = len(syms)
    key_mat = key_matrix(zma)
    val_mat = value_matrix(zma)

    xyzs = numpy.zeros((natms, 3))

    for key in range(1, natms):
        vals = val_mat[key][: min(key, 3)]
        keys = key_mat[key][: min(key, 3)]
        ref_xyzs = xyzs[list(keys)]
        xyz = util.vector.from_internals(*itertools.chain(*zip(vals, ref_xyzs)))
        xyzs[key] = xyz

    geo = geom.from_data(syms, xyzs)

    if zc_ is not None:
        geo = geom.undo_zmatrix_conversion(geo, zc_)
    elif not dummy:
        geo = geom.without_dummy_atoms(geo)

    return geo


def geometry_with_conversion_info(zma, zc_: ZmatConv = None):
    """Convert a Z-Matrix to a molecular geometry, along with a z-matrix conversion
    data structure describing the conversion

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :returns: automol molecular geometry data structure
    """
    zc_ = conversion_info(zma) if zc_ is None else zc_
    geo = geometry(zma, dummy=False, zc_=zc_)
    return geo, zc_


def rdkit_molecule(zma, gra=None, stereo=True):
    """Convert a z-matrix to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(zma))

    :param zma: molecular z-matrix
    :type zma: automol z-matrix data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :returns: the RDKit molecule
    """
    geo = geometry(zma, dummy=True)
    return geom.rdkit_molecule(geo, gra=gra, stereo=stereo)


def py3dmol_view(zma, gra=None, view=None, image_size=400):
    """Get a py3DMol view of this molecular z-matrix

    :param zma: molecular z-matrix
    :type zma: automol z-matrix data structure
    :param gra: A molecular graph, describing the connectivity, defaults to None
    :type gra: automol graph data structure, optional
    :param image_size: The image size, if creating a new view, defaults to 400
    :type image_size: int, optional
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    :return: A 3D view containing the molecule
    :rtype: py3Dmol.view
    """
    geo = geometry(zma, dummy=True)
    return geom.py3dmol_view(geo, gra=gra, view=view, image_size=image_size)


def display(zma, gra=None, vis_dists: bool = False, image_size: int = 400, view=None):
    """Display molecule to IPython using the RDKit visualizer

    :param zma: molecular z-matrix
    :type zma: automol z-matrix data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    :param vis_dists: Visualize the distance coordinates, setting these as the "bonds"
        in the visualization
    :param image_size: The image size, if creating a new view, defaults to 400
    :param view: An existing 3D view to append to, defaults to None
    :type view: py3Dmol.view, optional
    """
    geo = geometry(zma, dummy=True)
    vis_bkeys = list(distance_coordinates(zma).values()) if vis_dists else None
    return geom.display(
        geo, gra=gra, view=view, image_size=image_size, vis_bkeys=vis_bkeys
    )


# # derived properties
def distance(zma, key1, key2, angstrom=False):
    """Measure the distance between two atoms defined in a Z-Matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param key1: key of atom 1 in the pair to be measured
    :type key1: int
    :param key2: key of atom 2 in the pair to be measured
    :type key2: int
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    """
    geo = geometry(zma, dummy=True)
    return geom.distance(geo, key1, key2, angstrom=angstrom)


def central_angle(zma, key1, key2, key3, degree=False):
    """Measure the angle inscribed by three atoms in a Z-Matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param key1: key of atom 1 in the triplet to be measured
    :type key1: int
    :param key2: key of atom 2 in the triplet to be measured
    :type key2: int
    :param key3: key of atom 3 in the triplet to be measured
    :type key3: int
    :param degree: parameter to control radian->degree conversion
    :type degree: bool
    """
    geo = geometry(zma, dummy=True)
    return geom.central_angle(geo, key1, key2, key3, degree=degree)


def dihedral_angle(zma, key1, key2, key3, key4, degree=False):
    """Measure the angle inscribed by four atoms in a molecular geometry.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param key1: key of atom 1 in the quartet to be measured
    :type key1: int
    :param key2: key of atom 2 in the quartet to be measured
    :type key2: int
    :param key3: key of atom 3 in the quartet to be measured
    :type key3: int
    :param key4: key of atom 4 in the quartet to be measured
    :type key4: int
    :param degree: parameter to control radian->degree conversion
    :type degree: bool
    """
    geo = geometry(zma, dummy=True)
    return geom.dihedral_angle(geo, key1, key2, key3, key4, degree=degree)


# # torsions
def torsion_coordinate_name(zma, key1, key2, zgra=None):
    """Obtain the name for dihedral coordinate about a torsion axis
    (rotational bond).

    :param zma: the z-matrix
    :type zma: automol Z-Matrix data structure
    :param key1: the first key in the torsion axis (rotational bond)
    :type key1: int
    :param key2: the second key in the torsion axis (rotational bond)
    :type key2: int
    :param gra: an automol graph data structure, aligned to the z-matrix;
        used to check connectivity when necessary
    :rtype: str
    """

    key = torsion_leading_atom(zma, key1, key2, zgra=zgra)
    name_mat = name_matrix(zma)
    name = name_mat[key][-1]

    return name


def torsion_leading_atom(zma, key1, key2, zgra=None):
    """Obtain the leading atom for a torsion coordinate about a torsion axis.

    The leading atom is the atom whose dihedral defines the torsional
    coordinate, which must always be the first dihedral coordinate
    for this bond.

    A bond is properly decoupled if all other dihedrals along this
    bond depend on the leading atom.

    :param zma: the z-matrix
    :type zma: automol Z-Matrix data structure
    :param key1: the first key in the torsion axis (rotational bond)
    :type key1: int
    :param key2: the second key in the torsion axis (rotational bond)
    :type key2: int
    :param gra: an automol graph data structure, aligned to the z-matrix;
        used to check connectivity when necessary
    :rtype: int
    """

    key_mat = key_matrix(zma)
    krs1 = [(key, row) for key, row in enumerate(key_mat) if row[:2] == (key1, key2)]
    krs2 = [(key, row) for key, row in enumerate(key_mat) if row[:2] == (key2, key1)]

    lead_key_candidates = []

    for krs in (krs1, krs2):
        if krs:
            keys, rows = zip(*krs)
            start_key = keys[0]
            assert all(row[-1] == start_key for row in rows[1:]), (
                "Torsion coordinate along bond "
                f"{key1:d}-{key2:d} not decoupled:\n"
                f"{string(zma, one_indexed=False)}"
            )
            if rows[0][-1] is not None:
                lead_key_candidates.append(start_key)

    if not lead_key_candidates:
        lead_key = None
    elif len(lead_key_candidates) == 1:
        lead_key = lead_key_candidates[0]
    else:
        # If we get to this point, then the z-matrix includes dihedrals across
        # the key1-key2 bond in both directions and we have to choose which
        # dihedral to use. This mans there will be two lead_key_candidates.
        zgra = graph(zma) if zgra is None else zgra

        # Let key0 be the lead key and let (key1, key2, key3) be its key row in
        # the z-matrix. For the torsion coordinate, key0-key1-key2-key3 should
        # all be connected in a line. For subsidiary dihedral coordinates, key3
        # will be connected to key1 instead of key2.
        # A simple solution is therefore to choose the lead key based on
        # whether or not key2 and key3 are connected, which is what this code
        # does.
        bnd_keys = graph_base.bond_keys(zgra)
        lead_key = next(
            (k for k in lead_key_candidates if frozenset(key_mat[k][-2:]) in bnd_keys),
            None,
        )

        # If that fails, choose the key that appears earlier. It's possible
        # that it would be better to choose the later one, in which case we
        # would replace the min() here with a max().
        if lead_key is None:
            lead_key = min(lead_key_candidates)

    return lead_key


# repulsion energy
def has_low_relative_repulsion_energy(
    zma, ref_zma, model: str = "exp6", thresh=40.0
) -> bool:
    """Identify whether a z-matrix has low repulsion energy relative to a reference
    z-matrix

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param ref_zma: A reference z-matrix to compare to
    :type ref_zma: automol zmat data structure
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :param thresh: Threshold for excess repulsion energy (kcal/mol), defaults to 40.0
    :type thresh: float, optional
    :return: `True` if it does, `False` if it doesn't
    :rtype: bool
    """
    geo = geometry(zma)
    ref_geo = geometry(ref_zma)
    return geom.has_low_relative_repulsion_energy(
        geo, ref_geo, model=model, thresh=thresh
    )


def total_repulsion_energy(zma, model: str = "exp6") -> float:
    """Calculate the z-matrix's total repulsion energy using a model potential

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param model: The model potential to use, "exp6" (default) or "lj_12_6"
    :type model: str, optional
    :return: The repulsion energy
    :rtype: float
    """
    geo = geometry(zma)
    return geom.total_repulsion_energy(geo, model=model)


def dummy_coordinate_names(
    zma: object, with_lin_frag_dih: bool = True
) -> tuple[str, ...]:
    """Obtain names of all coordinates associated with dummy atoms
    defined in the V-Matrix.

    When one of two fragments is linear, the dihedral off of the dummy atom controlling
    their relative orientation is a bad coordinate and should be frozen along with other
    dummy coordinates. To include such coordinates in the list of dummy coordinates, set
    `with_lin_frag_dih = True`.

    :param vma: V-Matrix
    :param with_lin_frag_dih: Return linear fragment dihedrals along with the other
        dummy coordinates?
    :return: Name of coordinates
    """
    name_mat = numpy.array(name_matrix(zma))
    dum_keys = dummy_keys(zma)
    dum_names = []
    for dum_key in dum_keys:
        for col_idx in range(3):
            dum_name = next(
                filter(lambda x: x is not None, name_mat[dum_key:, col_idx])
            )
            dum_names.append(dum_name)

    # If requested, add in linear fragment dihedrals, if present
    if with_lin_frag_dih:
        geo = geometry(zma, dummy=True)
        keys = zmatrix_keys(zma)
        dum_keys = dummy_keys(zma)
        src_dct = dummy_source_dict(zma, dir_=False)
        nkeys_dct = neighbor_keys(zma)
        dih_dct = dihedral_angle_coordinates(zma)

        for dum_key in dum_keys:
            src_key = src_dct.get(dum_key)
            src_nkeys = nkeys_dct.get(src_key) | {src_key}
            frag1_keys = (frozenset(keys[:src_key]) | src_nkeys) - {dum_key}
            frag2_keys = (frozenset(keys[src_key:]) | src_nkeys) - {dum_key}
            has_lin_frag = any(
                geom.is_linear(geo, idxs=ks) for ks in (frag1_keys, frag2_keys)
            )
            if has_lin_frag:
                bad_dih_name = next(
                    (n for n, c in dih_dct.items() if c[-1] == dum_key), None
                )
                if bad_dih_name not in dum_names and bad_dih_name is not None:
                    dum_names.append(bad_dih_name)

    dum_names = tuple(dum_names)

    return dum_names
