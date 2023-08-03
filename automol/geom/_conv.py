""" Level 4 geometry functions
"""

import itertools

import numpy
from phydat import phycon

import automol.amchi.base
import automol.graph
import automol.inchi.base
import automol.zmat.base
from automol import util
from automol.extern import rdkit_
from automol.geom import _pyx2z
from automol.geom.base import (
    central_angle,
    coordinates,
    count,
    dihedral_angle,
    distance,
    from_subset,
    insert,
    is_atom,
    move_atom,
    rotate,
    symbols,
    translate,
)


# # conversions
def graph(geo, stereo=True):
    """Generate a molecular graph from the molecular geometry that has
    information about bond connectivity and if requested, stereochemistry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: automol molecular graph data structure
    """
    gra = connectivity_graph(geo)
    if stereo:
        gra = automol.graph.set_stereo_from_geometry(gra, geo)

    return gra


def connectivity_graph(geo, rqq_bond_max=3.45, rqh_bond_max=2.6, rhh_bond_max=1.9):
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
            else (dist < rqh_bond_max)
            if "H" in (sym1, sym2)
            else (dist < rhh_bond_max)
            if (sym1 == "H" and sym2 == "H")
            else (dist < rqq_bond_max)
        )

    idxs = range(len(xyzs))
    atm_symb_dct = dict(enumerate(symbs))
    bnd_keys = tuple(
        map(frozenset, filter(_are_bonded, itertools.combinations(idxs, r=2)))
    )

    gra = automol.graph.from_data(atm_symb_dct=atm_symb_dct, bnd_keys=bnd_keys)

    return gra


def zmatrix(geo, ts_bnds=()):
    """Generate a corresponding Z-Matrix for a molecular geometry
    using internal autochem procedures.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    :returns: automol Z-Matrix data structure
    """
    zma, _, _ = zmatrix_with_conversion_info(geo, ts_bnds=ts_bnds)
    return zma


def zmatrix_with_conversion_info(geo, ts_bnds=()):
    """Generate a corresponding Z-Matrix for a molecular geometry
    using internal autochem procedures.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    :returns: automol Z-Matrix data structure, Z-Matrix atom ordering, and
        a dictionary mapping linear atoms onto their associated dummy atoms
    """

    if ts_bnds:
        raise NotImplementedError

    if is_atom(geo):
        symbs = symbols(geo)
        key_mat = [[None, None, None]]
        val_mat = [[None, None, None]]
        zma = automol.zmat.base.from_data(symbs, key_mat, val_mat)
        zma_keys = [0]
        dummy_key_dct = {}
    else:
        geo, dummy_key_dct = insert_dummies_on_linear_atoms(geo)
        gra = connectivity_graph(geo)
        bnd_keys = tuple(dummy_key_dct.items())
        ord_dct = {k: 0 for k in bnd_keys}
        gra = automol.graph.add_bonds(gra, bnd_keys, ord_dct=ord_dct)
        vma, zma_keys = automol.graph.vmat.vmatrix(gra)
        geo = from_subset(geo, zma_keys)
        zma = automol.zmat.base.from_geometry(vma, geo)

    return zma, zma_keys, dummy_key_dct


def x2z_zmatrix(geo, ts_bnds=()):
    """Generate a corresponding Z-Matrix for a molecular geometry
    using x2z interface.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    """

    if is_atom(geo):
        symbs = automol.geom.base.symbols(geo)
        key_mat = [[None, None, None]]
        val_mat = [[None, None, None]]
        zma = automol.zmat.base.from_data(symbs, key_mat, val_mat)
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        zma = _pyx2z.to_zmatrix(x2m)
    zma = automol.zmat.base.standard_form(zma)

    return zma


def amchi(geo, stereo=True):
    """Generate an AMChI string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    ach, _ = amchi_with_sort(geo, stereo=stereo)

    return ach


def amchi_with_sort(geo, stereo=True, gra=None):
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
    ach, ach_idx_dcts = automol.graph.amchi_with_indices(gra, stereo=stereo)
    nums_lst = tuple(map(util.dict_.keys_sorted_by_value, ach_idx_dcts))
    return ach, nums_lst


def inchi(geo, stereo=True):
    """Generate an InChI string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param stereo: parameter to include stereochemistry information
    :type stereo: bool
    :rtype: str
    """
    ich, _ = inchi_with_sort(geo, stereo=stereo)

    return ich


def inchi_with_sort(geo, stereo=True, gra=None):
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
    gra = connectivity_graph(geo) if gra is None else gra
    if not stereo:
        geo = None
        geo_idx_dct = None
    else:
        geo_idx_dct = dict(enumerate(range(count(geo))))
    ich, nums_lst = automol.graph.inchi_with_sort_from_geometry(
        gra=gra, geo=geo, geo_idx_dct=geo_idx_dct
    )

    return ich, nums_lst


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

    # old implementation
    # if automol.graph.has_resonance_bond_stereo(gra):
    #     chi_, nums_lst = amchi_with_sort(geo, stereo=stereo, gra=gra)
    # else:
    #     chi_, nums_lst = inchi_with_sort(geo, stereo=stereo, gra=gra)
    #     # If the InChI has mobile hydrogens, revert back to AMChI
    #     if automol.amchi.base.has_mobile_hydrogens(chi_):
    #         chi_, nums_lst = amchi_with_sort(geo, stereo=stereo, gra=gra)

    # new implementation
    chi_, nums_lst = inchi_with_sort(geo, stereo=stereo, gra=gra)
    if automol.graph.inchi_is_bad(gra, chi_):
        chi_, nums_lst = amchi_with_sort(geo, stereo=stereo, gra=gra)

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
    smi = automol.graph.base.smiles(gra, stereo=stereo, res_stereo=res_stereo)
    return smi


def rdkit_molecule(geo, gra=None):
    """Convert a geometry to an RDKit molecule.

    This is mainly useful for quick visualization with IPython, which can
    be done as follows:
    >>> from IPython.display import display
    >>> display(rdkit_molecule(geo))

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    :returns: the RDKit molecule
    """
    rdkit_.turn_3d_visualization_on()
    gra = graph(geo) if gra is None else gra
    return rdkit_.from_geometry_with_graph(geo, gra)


def view(geo, gra=None, image_size=400, vec_xyz=None, vec_xyz_origin=(0, 0, 0)):
    """Get a py3DMol view of this molecular geometry

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity, defaults to None
    :type gra: automol graph data structure, optional
    :param image_size: The image size in pixels, defaults to 400
    :type image_size: int, optional
    :param vec_xyz: Optionally, show a vector on the view, defaults to None
    :type vec_xyz: Tuple[float, float, float], optional
    :param vec_xyz_origin: Origin for `vec_xyz`, defaults to (0, 0, 0)
    :type vec_xyz_origin: Tuple[float, float, float], optional
    :return: A viewer object
    :rtype: py3DMol.GLViewer
    """
    rdm = rdkit_molecule(geo, gra=gra)
    view = rdkit_.to_3d_view(rdm, image_size=image_size)
    if vec_xyz is not None:
        view.addArrow(
            {
                "start": dict(zip("xyz", vec_xyz_origin)),
                "end": dict(zip("xyz", numpy.add(vec_xyz, vec_xyz_origin))),
            }
        )
    return view


def display(geo, gra=None, image_size=400):
    """Display molecule to IPython using the RDKit visualizer

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: A molecular graph, describing the connectivity
    :type gra: automol graph data structure
    """
    view = view(geo, gra=gra, image_size=image_size)
    return view.show()


# # derived properties
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

    gra = connectivity_graph(geo) if gra is None else gra
    ngb_idxs_dct = automol.graph.atoms_neighbor_atom_keys(gra)

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


def closest_unbonded_atoms(geo, gra=None):
    """Determine which pair of unbonded atoms in a molecular geometry
    are closest together.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :type gra: automol graph data structure
    :rtype: (frozenset(int), float)
    """

    gra = connectivity_graph(geo) if gra is None else gra
    atm_keys = automol.graph.atom_keys(gra)
    bnd_keys = automol.graph.bond_keys(gra)
    poss_bnd_keys = set(map(frozenset, itertools.combinations(atm_keys, r=2)))

    # The set of candidates includes all unbonded pairs of atoms
    cand_bnd_keys = poss_bnd_keys - bnd_keys

    min_bnd_key = None
    min_dist_val = 1000.0
    for bnd_key in cand_bnd_keys:
        dist_val = distance(geo, *bnd_key)
        if dist_val < min_dist_val:
            min_dist_val = dist_val
            min_bnd_key = bnd_key

    return min_bnd_key, min_dist_val


def external_symmetry_factor(geo, chiral_center=True):
    """Obtain the external symmetry factor for a geometry using x2z interface
    which determines the initial symmetry factor and then divides by the
    enantiomeric factor.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: float
    """

    if is_atom(geo):
        ext_sym_fac = 1.0
    else:
        oriented_geom = _pyx2z.to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
        if oriented_geom.is_enantiomer() and chiral_center:
            ext_sym_fac *= 0.5

    return ext_sym_fac


def x2z_torsion_coordinate_names(geo, ts_bnds=()):
    """Generate a list of torsional coordinates using x2z interface. These
    names corresond to the Z-Matrix generated using the same algorithm.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    :rtype: tuple(str)
    """

    symbs = symbols(geo)
    if len(symbs) == 1:
        names = ()
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        names = _pyx2z.zmatrix_torsion_coordinate_names(x2m)

        zma = _pyx2z.to_zmatrix(x2m)
        name_dct = automol.zmat.base.standard_names(zma)
        names = tuple(map(name_dct.__getitem__, names))

    return names


def x2z_atom_ordering(geo, ts_bnds=()):
    """Generate a dictionary which maps the order of atoms from the input
    molecular geometry to the order of atoms of the resulting Z-Matrix
    that is generated by the x2z interface.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param ts_bnds: keys for the breaking/forming bonds in a TS
    :type ts_bnds: tuple(frozenset(int))
    :rtype: dict[int: int]
    """

    symbs = symbols(geo)
    if len(symbs) == 1:
        idxs = {0: 0}
    else:
        x2m = _pyx2z.from_geometry(geo, ts_bnds=ts_bnds)
        idxs = _pyx2z.zmatrix_atom_ordering(x2m)

    return idxs


# # derived operations
def insert_dummies_on_linear_atoms(geo, lin_idxs=None, gra=None, dist=1.0, tol=5.0):
    """Insert dummy atoms over linear atoms in the geometry.

    :param geo: the geometry
    :type geo: automol molecular geometry data structure
    :param lin_idxs: the indices of the linear atoms; if None, indices are
        automatically determined from the geometry based on the graph
    :type lin_idxs: tuple(int)
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :type gra: automol molecular graph data structure
    :param dist: distance of dummy atom from the linear atom, in angstroms
    :type dist: float
    :param tol: the tolerance threshold for linearity, in degrees
    :type tol: float
    :returns: geometry with dummy atoms inserted, along with a dictionary
        mapping the linear atoms onto their associated dummy atoms
    :rtype: automol molecular geometry data structure
    """

    lin_idxs = linear_atoms(geo) if lin_idxs is None else lin_idxs
    gra = connectivity_graph(geo) if gra is None else gra

    dummy_ngb_idxs = set(automol.graph.dummy_atoms_neighbor_atom_key(gra).values())
    atoms_already_have_dummys = dummy_ngb_idxs & set(lin_idxs)
    assert not atoms_already_have_dummys, (
        "Attempting to add dummy atoms on atoms that already have them"
        f"{atoms_already_have_dummys}"
    )

    ngb_idxs_dct = automol.graph.atoms_sorted_neighbor_atom_keys(gra)

    xyzs = coordinates(geo, angstrom=True)

    def _perpendicular_direction(idxs):
        """find a nice perpendicular direction for a series of linear atoms"""
        triplets = []
        for idx in idxs:
            for n1idx in ngb_idxs_dct[idx]:
                for n2idx in ngb_idxs_dct[n1idx]:
                    if n2idx != idx:
                        ang = central_angle(geo, idx, n1idx, n2idx, degree=True)
                        if numpy.abs(ang - 180.0) > tol:
                            triplets.append((idx, n1idx, n2idx))

        if triplets:
            idx1, idx2, idx3 = min(triplets, key=lambda x: x[1:])
            xyz1, xyz2, xyz3 = map(xyzs.__getitem__, (idx1, idx2, idx3))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            r23 = util.vec.unit_direction(xyz2, xyz3)
            direc = util.vec.orthogonalize(r12, r23, normalize=True)
        else:
            if len(idxs) > 1:
                idx1, idx2 = idxs[:2]
            else:
                (idx1,) = idxs
                idx2 = ngb_idxs_dct[idx1][0]
                # idx2, = ngb_idxs_dct[idx1]

            xyz1, xyz2 = map(xyzs.__getitem__, (idx1, idx2))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            for i in range(3):
                disp = numpy.zeros((3,))
                disp[i] = -1.0
                alt = numpy.add(r12, disp)
                direc = util.vec.unit_perpendicular(r12, alt)
                if numpy.linalg.norm(direc) > 1e-2:
                    break

        return direc

    # partition the linear atoms into adjacent groups, to be handled together
    lin_idxs_lst = sorted(
        map(
            sorted,
            util.equivalence_partition(lin_idxs, lambda x, y: x in ngb_idxs_dct[y]),
        )
    )

    dummy_key_dct = {}

    for idxs in lin_idxs_lst:
        direc = _perpendicular_direction(idxs)
        for idx in idxs:
            xyz = numpy.add(xyzs[idx], numpy.multiply(dist, direc))
            dummy_key_dct[idx] = count(geo)

            geo = insert(geo, "X", xyz, angstrom=True)

    return geo, dummy_key_dct


def insert_dummies(geo, dummy_key_dct, gra=None, dist=1.0, tol=5.0):
    """Insert dummy atoms over atoms in a geometry in a particular order.

    :param geo: the geometry
    :type geo: automol molecular geometry data structure
    :param dummy_key_dct: the linear atoms and the desired positions of the
        dummy atoms for each; linear atom indexing should follow what they
        *will* be after the dummy atoms are moved to the appropriate
        positions
    :param dummy_key_dct: dict
    :param dist: distance of dummy atom from the linear atom, in angstroms
    :type dist: float
    :param tol: the tolerance threshold for linearity, in degrees
    :type tol: float
    :returns: geometry with dummy atoms inserted, along with a dictionary
        mapping the linear atoms onto their associated dummy atoms
    :rtype: automol molecular geometry data structure
    """
    if dummy_key_dct:
        lin_keys, dum_keys = zip(*sorted(dummy_key_dct.items(), key=lambda x: x[1]))
        dum_keys = numpy.array(list(dum_keys))
        lin_idxs = [k - sum(k > dum_keys) for k in lin_keys]
        geo, orig_dummy_key_dct = insert_dummies_on_linear_atoms(
            geo, lin_idxs=lin_idxs, gra=gra, dist=dist, tol=tol
        )

        for lin_idx, lin_key in zip(lin_idxs, lin_keys):
            orig_idx = orig_dummy_key_dct[lin_idx]
            new_idx = dummy_key_dct[lin_key]
            geo = move_atom(geo, orig_idx, new_idx)

    return geo


def change_zmatrix_row_values(
    geo,
    idx,
    dist=None,
    idx1=None,
    ang=None,
    idx2=None,
    dih=None,
    idx3=None,
    angstrom=True,
    degree=True,
    gra=None,
):
    """Change the z-matrix coordinates of a given atom, shifting those
    connected to it accordingly.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx: the atom to be shifted
    :type idx: int
    :param dist: the distance coordinate; if None, use the current distance
    :type dist: float
    :param idx1: the atom used to specify the distance coordinate
    :type idx1: int
    :param ang: the central angle coordinate
    :type ang: float
    :param idx2: the atom used to specify the central angle coordinate
    :type idx2: int
    :param dih: the dihedral angle coordinate
    :type ang: float
    :param idx3: the atom used to specify the dihedral angle coordinate
    :type idx3: int
    :param angstrom: are distances in angstrom? If not, assume bohr.
    :type angstrom: bool
    :param degree: are angles in degrees? If not, assume radians.
    :type degree: bool
    :param gra: molecular graph for tracking connectivity (will be
        generated if None)
    :type gra: automol molecular graph data structure
    :rtype: automol geometry data structure
    """
    gra = gra if gra is not None else connectivity_graph(geo)
    dist = (
        dist
        if dist is not None or idx1 is None
        else (distance(geo, idx, idx1, angstrom=angstrom))
    )
    ang = (
        ang
        if ang is not None or idx2 is None
        else (central_angle(geo, idx, idx1, idx2, degree=degree))
    )
    dih = (
        dih
        if dih is not None or idx3 is None
        else (dihedral_angle(geo, idx, idx1, idx2, idx3, degree=degree))
    )

    dist = dist if not angstrom else dist * phycon.ANG2BOHR
    ang = ang if not degree else ang * phycon.DEG2RAD
    dih = dih if not degree else dih * phycon.DEG2RAD
    tol = 2.0 * phycon.DEG2RAD
    lin = numpy.pi

    idxs = automol.graph.branch_atom_keys(gra, idx1, idx)

    # Note that the coordinates will change throughout, but the change
    # shouldn't affect the operations so we can work in terms of the inital set
    xyzs = coordinates(geo)

    if idx1 is not None:
        # First, adjust the distance
        vec0 = numpy.subtract(xyzs[idx], xyzs[idx1])
        vec = dist * vec0 / numpy.linalg.norm(vec0)
        disp = vec - vec0
        geo = translate(geo, disp, idxs=idxs)

    if idx2 is not None:
        assert idx1 is not None, "idx1 is required if changing an angle"
        # Second, adjust the angle
        ang0 = central_angle(geo, idx, idx1, idx2)
        ang_diff = ang - ang0
        # If idx-idx1-idx2 is not linear, use the normal to this plane,
        # otherwise, use the normal to the idx1-idx2-idx3 plane
        if not numpy.abs(ang0) < tol or numpy.abs(ang0 - lin) < tol:
            axis = util.vec.unit_perpendicular(
                xyzs[idx], xyzs[idx2], orig_xyz=xyzs[idx1]
            )
        else:
            axis = util.vec.unit_perpendicular(
                xyzs[idx1], xyzs[idx3], orig_xyz=xyzs[idx1]
            )
        # I don't know how to figure out which way to rotate, so just try both
        # and see which one works
        for dang in [-ang_diff, ang_diff]:
            geo_ = rotate(geo, axis, dang, orig_xyz=xyzs[idx1], idxs=idxs)
            ang_out = central_angle(geo_, idx, idx1, idx2)
            abs_diff = numpy.abs(
                numpy.mod(ang, 2 * numpy.pi) - numpy.mod(ang_out, 2 * numpy.pi)
            )
            ang_comp = numpy.abs(numpy.pi - numpy.abs(abs_diff - numpy.pi))
            if ang_comp < tol:
                geo = geo_
                break

    if idx3 is not None:
        assert idx1 is not None, "idx1 is required if changing a dihedral"
        assert idx2 is not None, "idx2 is required if changing a dihedral"
        # Third, adjust the dihedral
        dih0 = dihedral_angle(geo, idx, idx1, idx2, idx3)
        ddih = dih - dih0
        axis = numpy.subtract(xyzs[idx1], xyzs[idx2])
        geo = rotate(geo, axis, ddih, orig_xyz=xyzs[idx1], idxs=idxs)

    return geo
