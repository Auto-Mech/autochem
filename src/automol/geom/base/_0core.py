"""
Core functions defining the geometry data type
"""

import functools
import itertools
from collections.abc import Sequence
from typing import List, Optional

import more_itertools as mit
import numpy
import pyparsing as pp
from numpy.typing import ArrayLike
from pyparsing import pyparsing_common as ppc

from phydat import phycon, ptab

from ... import form, util

AXIS_DCT = {"x": 0, "y": 1, "z": 2}

CHAR = pp.Char(pp.alphas)
SYMBOL = pp.Combine(CHAR + pp.Opt(CHAR))
XYZ = ppc.fnumber * 3
XYZ_LINE = SYMBOL("symb") + XYZ("xyz") + pp.Optional(XYZ)("mode_xyz")


# # constructors
def from_data(symbs, xyzs, angstrom=False, check=True):
    """Build a geometry data structure from atomic symbols and coordinates.

    format:
        geo = (((sym1, (xcoord1, ycoord1, zcoord1)),...
               ((symn, (xcoordn, ycoordn, zcoordn)))

    :param symbs: atomic symbols of the atoms
    :type symbs: tuple(str)
    :param xyzs: xyz coordinates of the atoms
    :type xyzs: tuple(float)
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    :param check: check that argument values make sense?
    :type check: bool
    """

    symbs = list(map(ptab.to_symbol, symbs))
    natms = len(symbs)

    xyzs = numpy.array(xyzs, dtype=float)
    if check:
        assert numpy.ndim(xyzs) == 2 and numpy.shape(xyzs) == (natms, 3)

    xyzs = xyzs if not angstrom else numpy.multiply(xyzs, phycon.ANG2BOHR)
    xyzs = [tuple(map(float, xyz[:3])) for xyz in xyzs]
    geo = tuple(zip(symbs, xyzs))

    return geo


def subgeom(geo, idxs):
    """Generate a new molecular geometry from a subset of the atoms in an
    input geometry.

    (Rename this and put it under operations?)

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idxs: indices representing the subset of atoms
    :type idxs: tuple(int)
    :rtype: automol moleculer geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    symbs = list(map(symbs.__getitem__, idxs))
    xyzs = list(map(xyzs.__getitem__, idxs))

    return from_data(symbs, xyzs)


# # getters
def symbols(geo, idxs: List[int] = None) -> List[str]:
    """Obtain the atomic symbols atoms in the molecular geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idxs: indices of atoms to obtain information for
    :type idxs: List[int]
    :returns: The list of atomic symbols
    :rtype: List[str]
    """
    if geo:
        symbs, _ = zip(*geo)
    else:
        symbs = ()

    return symbs if idxs is None else tuple(map(symbs.__getitem__, idxs))


def coordinates(geo, idxs=None, angstrom=False):
    """Obtain the Cartesian coordinates of atoms in the molecular geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idxs: indexs of atoms to obtain information for
    :type idxs: tuple(int)
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    :rtype: tuple(tuple(float))
    """

    idxs = list(range(count(geo))) if idxs is None else idxs
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    xyzs = xyzs if not angstrom else numpy.multiply(xyzs, phycon.BOHR2ANG)
    xyzs = tuple(map(tuple, map(xyzs.__getitem__, idxs)))

    return xyzs


# # setters
def set_coordinates(geo, xyz_dct, angstrom=False):
    """Set coordinate values for the molecular geometry,
    using a dictionary by index.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param xyz_dct: new xyz values for a set of indices
    :type xyz_dct: dict[int: tuple(float)]
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    natms = len(symbs)
    assert all(idx in range(natms) for idx in xyz_dct)

    conv = phycon.ANG2BOHR if angstrom else 1.0
    xyzs = [
        numpy.multiply(xyz_dct[idx], conv) if idx in xyz_dct else xyz
        for idx, xyz in enumerate(xyzs)
    ]

    return from_data(symbs, xyzs)


# # I/O
def string(geo, angstrom: bool = True, mode: Optional[ArrayLike] = None):
    """Write a molecular geometry to a string:
       symb1  xyz1 xyz2 xyz3
       symbn  xyzn xyzn xyzn

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param angstrom: parameter to control coordinate conversion to Angstrom
    :param mode: A vibrational mode or molecular motion to visualize
    :rtype: str
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)

    natms = len(symbs)

    symb_strs = [f"{s:2s}" for s in symbs]
    xyz_strs = [f"{x[0]:10.6f} {x[1]:10.6f} {x[2]:10.6f}" for x in xyzs]
    lines = [" ".join([s, x]) for s, x in zip(symb_strs, xyz_strs)]

    # If requested, include the mode in the string
    if mode is not None:
        mode = numpy.reshape(mode, (natms, 3))
        mode_strs = [f"{m[0]:10.6f} {m[1]:10.6f} {m[2]:10.6f}" for m in mode]
        lines = ["  ".join([s, x]) for s, x in zip(lines, mode_strs)]

    geo_str = "\n".join(lines)

    return geo_str


def xyz_string(geo, comment="", mode: Optional[ArrayLike] = None) -> str:
    """Write a molecular geometry to a string:
       natom
       comment
       symb1  xyz1 xyz2 xyz3
       symbn  xyzn xyzn xyzn

    :param geo: molecular geometry
    :param comment: string to place in the comment line of string
    :param mode: A vibrational mode or molecular motion to visualize
    :rtype: str
    """
    geo_str = string(geo, angstrom=True, mode=mode)
    xyz_str = f" {count(geo):d}\n{comment:s}\n{geo_str:s}"
    return xyz_str


def xyz_trajectory_string(geos, comments=None):
    """Write a series of molecular geometries to trajectory file which
    is a string that collated by several xyz-file format geometry strings.

    :param geo_lst: list of molecular geometries
    :type geo_lst: tuple(automol geometry data structure)
    :param comments: list of comments for each of the molecular geometries
    :type comments: tuple(str)
    :rtype: str
    """
    ngeos = len(geos)
    comments = [""] * ngeos if comments is None else comments

    assert len(comments) == ngeos
    xyz_strs = [xyz_string(g, comment=c) for g, c in zip(geos, comments)]

    return "\n".join(xyz_strs)


def from_string(geo_str, angstrom=True):
    """Read a Cartesian molecular geometry from a string comprised
    of just the atomic symbols and coordinates.

    :param geo_str: string containing the geometry
    :type geo_str: str
    :param angstrom: parameter to control coordinate conversion to Angstrom
    :type angstrom: bool
    :rtype: automol geometry data structure
    """
    geo, *_ = parse_xyz_block(geo_str, angstrom=angstrom)
    return geo


def from_xyz_string(xyz_str):
    """Read a Cartesian molecular geometry from a string that matches the
    format of a string of a standard .xyz file.

    :param xyz_str: string obtained from reading the .xyz file
    :type xyz_str: str
    :rtype: automol geometry data structure
    """
    geo, *_ = parse_xyz_block(xyz_str, angstrom=True)
    return geo


def from_xyz_string_with_mode(xyz_str):
    """Read a Cartesian molecular geometry from a string that matches the
    format of a string of a standard .xyz file.

    :param xyz_str: string obtained from reading the .xyz file
    :type xyz_str: str
    :rtype: automol geometry data structure
    """
    geo, _, mode = parse_xyz_block(xyz_str, angstrom=True)
    return geo, mode


def xyz_string_comment(xyz_str):
    """Read the comment line of a string of a standard .xyz file.

    :param xyz_str: string obtained from reading the .xyz file
    :type xyz_str: str
    :rtype: str
    """
    _, comment, _ = parse_xyz_block(xyz_str, angstrom=True)
    return comment


def from_xyz_trajectory_string(geo_str: str):
    """Read a series of molecular geometries from a trajectory string.

    :param geo_str: string containing the geometry
    :type geo_str: str
    :rtype: (tuple(automol geometry data structure), tuple(str))
    """
    geo_str = geo_str.strip()

    natm_str, *_ = geo_str.splitlines()
    natm = int(natm_str)

    blocks = ["\n".join(b) for b in mit.chunked(geo_str.splitlines(), natm+2)]
    geos = []
    comments = []
    for block in blocks:
        geo, comment, *_ = parse_xyz_block(block, angstrom=True)
        geos.append(geo)
        comments.append(comment)

    return tuple(zip(geos, comments, strict=True))


def parse_xyz_block(
    xyz_block: str, *, angstrom: bool = True
) -> tuple[object, str | None, numpy.ndarray | None]:
    """Parse an XYZ block string.

    :param xyz_block: XYZ block string
    :param angstrom: Whether the units are angstroms
    :return: Geometry, comment (if present), and mode (if present)
    """
    lines = xyz_block.strip().splitlines()

    # If line 1 is an integer, this is a standard XYZ file
    natm = None
    comment = None
    if lines[0].strip().isdigit():
        natm_line, comment, *lines = lines
        natm = int(natm_line)
        lines = lines[:natm]

    symbs = []
    xyzs = []
    mode = []
    for line in lines:
        res = XYZ_LINE.parse_string(line).as_dict()
        if res:
            symbs.append(res.get("symb"))
            xyzs.append(res.get("xyz"))
            mode.append(res.get("mode_xyz"))

    if any(mode) and not all(mode):
        msg = f"Invalid XYZ string with mode:\n{xyz_block}"
        raise ValueError(msg)

    geo = from_data(symbs, xyzs, angstrom=angstrom)
    mode = numpy.array(mode) if all(mode) else None
    return geo, comment, mode


def yaml_data(geo) -> list:
    """A yaml-friendly data format for the geometry

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    """
    geo = round_(geo)
    symbs = symbols(geo)
    xyzs = coordinates(geo)
    geo_yml = [[symb, *map(float, xyz)] for symb, xyz in zip(symbs, xyzs)]
    return geo_yml


def from_yaml_data(geo_yml):
    """Put a yaml-formatted geometry back into standard format

    :param geo_yml: A yaml-formatted molecular geometry
    :type geo_yml: list
    """
    symbs = [row[0] for row in geo_yml]
    xyzs = [row[1:] for row in geo_yml]
    return from_data(symbs, xyzs)


# # validation
def is_valid(geo):
    """Assess if a molecular geometry has the proper data structure.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: bool
    """

    ret = hasattr(geo, "__iter__")
    if ret:
        ret = all(hasattr(obj, "__len__") and len(obj) == 2 for obj in geo)
        if ret:
            symbs, xyzs = zip(*geo)
            try:
                from_data(symbs, xyzs)
            except AssertionError:
                ret = False

    return ret


# # conversions
def formula(geo):
    """Generate a stoichiometric formula dictionary from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :type: dict[str: int]
    """

    symbs = symbols(geo)
    fml = util.formula_from_symbols(symbs)

    return fml


def formula_string(geo):
    """Generate a stoichiometric formula string from a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :type: dict[str: int]
    """

    fml = formula(geo)
    fml_str = form.string(fml)

    return fml_str


# # properties
def count(geo):
    """Obtain the number of rows of the molecular geometry, which corresponds
    to the number of atoms in the geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :rtype: int
    """
    return len(geo)


def is_atom(geo):
    """Determine if the molecular geometry corresponds to an atom.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: bool
    """

    return count(geo) == 1


def is_diatomic(geo):
    """Determine if the molecular geometry corresponds to diatomic species..

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: bool
    """

    return count(geo) == 2


def is_linear(geo, tol=2.0 * phycon.DEG2RAD, idxs: Sequence[int] | None = None):
    """Determine if the molecular geometry is linear.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param tol: tolerance of bond angle(s) for determing linearity
    :type tol: float
    :rtype: bool
    """
    geo = geo if idxs is None else subgeom(geo, idxs)

    ret = True

    if len(geo) == 1:
        ret = False
    elif len(geo) == 2:
        ret = True
    else:
        keys = range(len(symbols(geo)))
        for key1, key2, key3 in mit.windowed(keys, 3):
            cangle = numpy.abs(central_angle(geo, key1, key2, key3))
            if not (numpy.abs(cangle) < tol or numpy.abs(cangle - numpy.pi) < tol):
                ret = False
    return ret


def atom_count(geo, symb, match=True):
    """Count the number of atoms of a particular type in the geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param symb: atomic symbol
    :type symb: str
    :param match: count symbols that match the type?
    :type match: bool
    :rtype: tuple(int)
    """
    return len(atom_indices(geo, symb, match=match))


def electron_count(geo) -> int:
    """Count the number of electrons in the geometry

    :param geo: A molecular geometry
    :return: The number of electrons
    """
    return form.electron_count(formula(geo))


def atom_indices(geo, symb, match=True):
    """Obtain the indices of a atoms of a particular type in the geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param match: grab idxs that match given atom type
    :param symb: atomic symbol
    :type symb: str
    :param match: obtain indices of symbols that match the type?
    :type match: bool
    :rtype: tuple(int)
    """

    symbs = symbols(geo)
    idxs = tuple()
    for idx, symb_ in enumerate(symbs):
        if symb_ == symb and match:
            idxs += (idx,)
        elif symb_ != symb and not match:
            idxs += (idx,)

    return idxs


def dummy_atom_indices(geo):
    """Obtain the indices of a atoms of a particular type in the geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :rtype: tuple(int)
    """

    symbs = symbols(geo)
    dummy_idxs = [idx for idx, symb in enumerate(symbs) if not ptab.to_number(symb)]

    return tuple(dummy_idxs)


def masses(geo, amu=True):
    """Build a list of the atomic masses that corresponds to the list
    of atomic sybmols of a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(float)
    """

    symbs = symbols(geo)
    amas = list(map(ptab.to_mass, symbs))

    if not amu:
        amas = numpy.multiply(amas, phycon.AMU2EMASS)

    amas = tuple(amas)

    return amas


def mass_weight_vector(geo) -> numpy.ndarray:
    """Get the mass-weighting vector for a geometry.

    :param geo: The geometry
    :return: The mass-weighting vector
    """
    return numpy.sqrt(numpy.repeat(masses(geo), 3))


def total_mass(geo, amu=True):
    """Calculate the total mass of a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(float)
    """
    return sum(masses(geo, amu=amu))


def center_of_mass(geo):
    """Determine the center-of-mass for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: tuple(float)
    """

    xyzs = coordinates(geo)
    amas = masses(geo)
    cm_xyz = tuple(
        sum(numpy.multiply(xyz, ama) for xyz, ama in zip(xyzs, amas)) / sum(amas)
    )

    return cm_xyz


def mass_centered(geo):
    """Generate a new geometry where the coordinates of the input geometry
    have been translated to the center-of-mass.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :rtype: tuple(float)
    """
    return translate(geo, numpy.negative(center_of_mass(geo)))


def reduced_mass(geo1, geo2):
    """Calculate the reduced mass for two species.

    :param geo1: geometry of species 1 (Bohr)
    :type geo1: list(float)
    :param geo2: geometry of species 2 (Bohr)
    :type geo2: list(float)
    :return: reduced mass (amu)
    :rtype: float
    """

    mass1 = total_mass(geo1)
    mass2 = total_mass(geo2)

    return (mass1 * mass2) / (mass1 + mass2)


def inertia_tensor(geo, amu=True):
    """Build the moment-of-inertia tensor for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(tuple(float))
    """

    geo = mass_centered(geo)
    amas = masses(geo, amu=amu)
    xyzs = coordinates(geo)
    ine = tuple(
        map(
            tuple,
            sum(
                ama * (numpy.vdot(xyz, xyz) * numpy.eye(3) - numpy.outer(xyz, xyz))
                for ama, xyz in zip(amas, xyzs)
            ),
        )
    )

    return ine


def moments_of_inertia(geo, amu: bool = True, drop_null: bool = False):
    """Calculate the moments of inertia along the xyz axes
    (these not sorted in to A,B,C).

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(tuple(float))
    """
    moms, _ = rotational_analysis(geo, amu=amu, drop_null=drop_null)
    return moms


def principal_axes(geo, drop_null: bool = False):
    """Determine the principal axes of rotation for a molecular geometry.

    :param geo: molecular geometry
    :param amu: parameter to control electron mass -> amu conversion
    :rtype: tuple(tuple(float))
    """
    _, rot_axes = rotational_analysis(geo, drop_null=drop_null)
    return rot_axes


def aligned_to_principal_axes(geo):
    """Project the geometry onto principal axes of rotation, to standardize its
    orientation.

    :param geo: The geometry
    :return: The geometry with standard orientation
    """
    return transform_by_matrix(mass_centered(geo), principal_axes(geo), trans=False)


def rotational_analysis(
    geo, drop_null: bool = False, amu: bool = True
) -> tuple[tuple[float, ...], numpy.ndarray]:
    """Do a rotational analysis, generating the moments of inertia and principal axes.

    :param geo: A molecular geometry
    :param drop_null: Drop null rotations from the analysis? (1 if linear, 3 if atom)
    :param amu: Use atomic mass units instead of electron mass unit?, defaults to True
    :return: The eigenvalues and eigenvetors of the inertia tensor
    """
    ine = inertia_tensor(geo, amu=amu)
    eig_vals, eig_vecs = numpy.linalg.eigh(ine)
    # Sort the eigenvectors and values
    eig_vals = tuple(eig_vals)
    # Drop null rotations, if requested
    ndrop = 0 if not drop_null else 3 if is_atom(geo) else 1 if is_linear(geo) else 0
    assert numpy.allclose(eig_vals[:ndrop], 0.0), f"ndrop={ndrop} eig_vals={eig_vals}"
    return eig_vals[ndrop:], eig_vecs[:, ndrop:]


def translational_normal_modes(geo, mass_weight: bool = True) -> numpy.ndarray:
    """Calculate the translational normal modes.

    :param geo: The geometry
    :param mass_weight: Return mass-weighted normal modes?
    :return: The translational normal modes, as a 3N x 3 matrix
    """
    natms = count(geo)
    trans_coos = numpy.tile(numpy.eye(3), (natms, 1))
    if mass_weight:
        trans_coos *= mass_weight_vector(geo)[:, numpy.newaxis]
    return trans_coos


def rotational_normal_modes(geo, mass_weight: bool = True) -> numpy.ndarray:
    """Calculate the rotational normal modes.

    For each atom, the direction of rotation around a given axis is given as the cross
    product of its position with the axis.

    :param geo: The geometry
    :param mass_weight: Return mass-weighted normal modes?
    :return: The rotational normal modes, as a 3N x (3 or 2) matrix
    """
    _, rot_axes = rotational_analysis(geo, drop_null=True)
    xyzs = coordinates(mass_centered(geo))
    rot_coos = []
    for rot_axis in rot_axes.T:
        rot_coo = numpy.concatenate([numpy.cross(xyz, rot_axis) for xyz in xyzs])
        rot_coos.append(rot_coo)
    rot_coos = numpy.transpose(rot_coos)

    if mass_weight:
        rot_coos *= mass_weight_vector(geo)[:, numpy.newaxis]
    return rot_coos


def rotational_constants(geo, amu=True, drop_null: bool = False):
    """Calculate the rotational constants.
    (these not sorted in to A,B,C).

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(float)
    """
    moms = moments_of_inertia(geo, amu=amu, drop_null=drop_null)
    cons = numpy.divide(1.0, moms) / 4.0 / numpy.pi / phycon.SOL
    cons = tuple(cons)
    return cons


# # geometric measurements
def distance(geo, idx1, idx2, angstrom=False):
    """Measure the distance between two atoms in a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx1: index of atom 1 in the pair to be measured
    :type idx1: int
    :param idx2: index of atom 2 in the pair to be measured
    :type idx2: int
    :param angstrom: parameter to control conversion to Angstrom
    :type angstrom: bool
    :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    dist = util.vector.distance(xyz1, xyz2)
    dist *= phycon.BOHR2ANG if angstrom else 1
    return dist


def central_angle(geo, idx1, idx2, idx3, degree=False):
    """Measure the angle inscribed by three atoms in a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx1: index of atom 1 in the triplet to be measured
    :type idx1: int
    :param idx2: index of atom 2 in the triplet to be measured
    :type idx2: int
    :param idx3: index of atom 3 in the triplet to be measured
    :type idx3: int
    :param degree: parameter to control conversion to degree
    :type degree: bool
    :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    ang = util.vector.central_angle(xyz1, xyz2, xyz3)
    ang *= phycon.RAD2DEG if degree else 1
    return ang


def dihedral_angle(geo, idx1, idx2, idx3, idx4, degree=False):
    """Measure the angle inscribed by three atoms in a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx1: index of atom 1 in the quartet to be measured
    :type idx1: int
    :param idx2: index of atom 2 in the quartet to be measured
    :type idx2: int
    :param idx3: index of atom 3 in the quartet to be measured
    :type idx3: int
    :param idx4: index of atom 4 in the quartet to be measured
    :type idx4: int
    :param degree: parameter to control conversion to degree
    :type degree: bool
    :rtype: float
    """

    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    xyz4 = xyzs[idx4]
    dih = util.vector.dihedral_angle(xyz1, xyz2, xyz3, xyz4)
    dih *= phycon.RAD2DEG if degree else 1

    return dih


def measure(geo, coo, angstrom=False, degree=False):
    """Measure a coordinate value for distance, central angle, or dihedral angle

    :param geo: A molecular geometry
    :param coo: The indices defining the coordinate
        2 for distance, 3 for central angle, 4 for dihedral angle
    :param angstrom: Give distances in angstrom?, defaults to False
    :param degree: Give angles in degrees?, defaults to False
    :return: The measured value
    """
    assert 2 <= len(coo) <= 4, f"Invalid coordinate: {coo}"

    if len(coo) == 2:
        return distance(geo, *coo, angstrom=angstrom)
    if len(coo) == 3:
        return central_angle(geo, *coo, degree=degree)
    if len(coo) == 4:
        return dihedral_angle(geo, *coo, degree=degree)

    return None


def zmatrix_row_values(
    geo, idx, idx1=None, idx2=None, idx3=None, angstrom=True, degree=True
):
    """Get coordinate values for a row in a Z-Matrix.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx:
    :type idx: int
    :param idx1:
    :type idx1: int
    :type idx2: int
    :type idx3: int
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    :param degree: parameter to control conversion to degree
    :type degree: bool
    """

    val_row = [None, None, None]

    if idx1 is not None:
        val_row[0] = distance(geo, idx, idx1, angstrom=angstrom)

    if idx2 is not None:
        assert idx1 is not None
        val_row[1] = central_angle(geo, idx, idx1, idx2, degree=degree)

    if idx3 is not None:
        assert idx2 is not None
        val_row[2] = dihedral_angle(geo, idx, idx1, idx2, idx3, degree=degree)

    val_row = tuple(val_row)

    return val_row


# # binary functions
def join(geo1, geo2, dist_cutoff=3.0 * phycon.ANG2BOHR, theta=0.0, phi=0.0):
    """Join two molecular geometries together where the intermolecular
    separation and orientation can be specified.

    :param geo1: Molecular geometry 1
    :type geo1: automol molecular geometry data structure
    :param geo2: Molecular geometry 2
    :type geo2: automol molecular geometry data structure
    :param dist_cutoff: The closest allowable intermolecular separation
    :type: dist_cutoff: float
    :param theta: theta angle for intermolecular orientation
    :type theta: float
    :param phi: phi angle for intermolecular orientation
    :type phi: float
    :return: The combined geometry
    """

    if not geo1:
        symbs = symbols(geo2)
        xyzs = coordinates(geo2)
    elif not geo2:
        symbs = symbols(geo1)
        xyzs = coordinates(geo1)
    else:
        orient_vec = numpy.array(
            [
                numpy.sin(theta) * numpy.cos(phi),
                numpy.sin(theta) * numpy.sin(phi),
                numpy.cos(theta),
            ]
        )
        neg_orient_vec = -1.0 * orient_vec

        # get the correct distance apart
        geo1 = mass_centered(geo1)
        geo2 = mass_centered(geo2)
        ext1 = max(numpy.vdot(orient_vec, xyz) for xyz in coordinates(geo1))
        ext2 = max(numpy.vdot(neg_orient_vec, xyz) for xyz in coordinates(geo2))

        cm_dist = ext1 + dist_cutoff + ext2
        dist_grid = numpy.arange(cm_dist, 0.0, -0.1)
        for dist in dist_grid:
            trans_geo2 = translate(geo2, orient_vec * dist)
            min_dist = minimum_distance(geo1, trans_geo2)
            if numpy.abs(min_dist - dist_cutoff) < 0.1:
                break

        geo2 = trans_geo2

        # now, join them together
        symbs = symbols(geo1) + symbols(geo2)
        xyzs = coordinates(geo1) + coordinates(geo2)

    return from_data(symbs, xyzs)


def join_sequence(geos):
    """Join a sequence of molecular geometries

    :param geos: A sequence of molecular geometries
    :return: The combined geometry
    """
    return functools.reduce(join, geos)


def minimum_distance(geo1, geo2):
    """get the minimum distance between atoms in geo1 and those in geo2

    :param geo1: molecular geometry 1
    :type geo1: automol molecular geometry data structure
    :param geo2: molecular geometry 2
    :type geo2: automol molecular geometry data structure
    :rtype: float
    """

    xyzs1 = coordinates(geo1)
    xyzs2 = coordinates(geo2)
    return min(
        util.vector.distance(xyz1, xyz2)
        for xyz1, xyz2 in itertools.product(xyzs1, xyzs2)
    )


def permutation(geo, ref_geo, thresh=1e-4):
    """Determine the permutation of one geometry that reproduces another
    (if there isn't one -- the geometries are not aligned, return None).

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param ref_geo: molecular geometry
    :type ref_geo: automol molecular geometry data structure
    :param thresh: theshold for assessing if permutation exists
    :type thresh: float
    :rtype: tuple(int)
    """

    natms = count(geo)
    symbs = symbols(geo)
    xyzs = coordinates(geo)

    perm_idxs = [None] * natms
    for idx, (symb, xyz) in enumerate(zip(symbs, xyzs)):
        # Loop over atoms in the reference geometry with the same symbol
        ref_idxs = atom_indices(ref_geo, symb=symb)
        ref_xyzs = coordinates(ref_geo, idxs=ref_idxs)
        perm_idx = next(
            (
                ref_idx
                for ref_idx, ref_xyz in zip(ref_idxs, ref_xyzs)
                if util.vector.distance(xyz, ref_xyz) < thresh
            ),
            None,
        )
        perm_idxs[idx] = perm_idx

    perm_idxs = tuple(perm_idxs)

    if any(perm_idx is None for perm_idx in perm_idxs):
        perm_idxs = None

    return perm_idxs


# # add/remove/rearrange atoms
def insert(geo, symb, xyz, idx=None, angstrom=False):
    """Insert an atom into a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param symb: symbol of atom to add
    :type symb: str
    :param xyz: xyz coordinates of atom to add
    :type xyz: tuple(float)
    :param idx: index of geometry to place atom
    :type idx: int
    :rtype: automol geometry date structure
    """

    symbs = list(symbols(geo))
    xyzs = list(coordinates(geo, angstrom=angstrom))

    idx = idx if idx is not None else len(symbs)

    symbs.insert(idx, symb)
    xyzs.insert(idx, xyz)

    return from_data(symbs, xyzs, angstrom=angstrom)


def remove(geo, idxs=()):
    """Remove atoms by from the molecular geometry by their indices.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idxs: indices of atoms to remove
    :type idxs: tuple(int)
    :rtype: automol molecular geometry data structure
    """
    return tuple(row for i, row in enumerate(geo) if i not in idxs)


def without_dummy_atoms(geo):
    """Return a copy of the molecular geometry without dummy atoms."""

    symbs = symbols(geo)
    non_dummy_idxs = [idx for idx, symb in enumerate(symbs) if ptab.to_number(symb)]

    return subgeom(geo, non_dummy_idxs)


def reorder(geo, idx_dct):
    """Reorder the atoms of a molecular geometry using
    the mapping of an input dictionary.

    :param geo: The geometry
    :param idx_dct: The new order of the atoms, by index
    :type idx_dct: dict
    :rtype: automol geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)

    idxs = [idx for idx, _ in sorted(idx_dct.items(), key=lambda x: x[1])]
    assert len(symbs) == len(xyzs) == len(idxs)

    symbs = [symbs[idx] for idx in idxs]
    xyzs = [xyzs[idx] for idx in idxs]

    return from_data(symbs, xyzs)


def move_atom(geo, idx1, idx2):
    """Move an atom to a different position in the geometry

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idx1: index of the atom to be moved
    :type idx1: int
    :param idx2: new position that the atom should be moved to
    :type idx2: int
    :returns: the transformed geometry
    :rtype: molecular geometry
    """
    symbs = list(symbols(geo))
    xyzs = list(coordinates(geo))
    symbs.insert(idx2, symbs.pop(idx1))
    xyzs.insert(idx2, xyzs.pop(idx1))
    return from_data(symbs, xyzs)


def swap_coordinates(geo, idx1, idx2):
    """Swap the order of the coordinates of two atoms in a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idx1: index for one atom to swap coordinates
    :type idx1: int
    :param idx2: index for one atom to swap coordinates
    :type idx2: int
    :rtype: molecular geometry
    """

    geo = [list(x) for x in geo]
    geo[idx1], geo[idx2] = geo[idx2], geo[idx1]
    geo_swp = tuple(tuple(x) for x in geo)

    return geo_swp


# # transformations
def round_(geo, decimals=6):
    """Round the coordinates of this geometry to a certain number of decimal
    places

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param decimals: the number of decimals to round to
    :type decimals: int
    :rtype: automol molecular geometry data structure
    """
    symbs = symbols(geo)
    xyzs = numpy.round(coordinates(geo), decimals=decimals)
    return from_data(symbs, xyzs)


def displace(geo, disp_xyzs, angstrom=False):
    """Displace a molecular geometry along a particular mode.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param disp_xyzs: Displacement coordinates
    :param idxs: indices of atoms whose coordinates are to be translated
    :type idxs: tuple(int)
    :param angstrom: whether or not the translation is in angstrom
    :type angstrom: bool
    :rtype: automol molecular geometry data structure
    """
    disp_xyzs = numpy.reshape(disp_xyzs, (-1, 3))
    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)
    xyzs = numpy.add(xyzs, disp_xyzs)
    return from_data(symbs, xyzs, angstrom=angstrom)


def translate(geo, xyz, idxs=None, angstrom=False):
    """Translate the coordinates of a molecular geometry along
    a three-dimensiona vector.
    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param xyz: vector to translate along
    :type xyz: tuple(float)
    :param idxs: indices of atoms whose coordinates are to be translated
    :type idxs: tuple(int)
    :param angstrom: whether or not the translation is in angstrom
    :type angstrom: bool
    :rtype: automol molecular geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)

    natm = len(symbs)
    xyz = list(xyz)
    if idxs is None:
        disp = [xyz] * natm
    else:
        disp = [xyz if i in idxs else [0.0, 0.0, 0.0] for i in range(natm)]

    xyzs = numpy.add(xyzs, disp)
    return from_data(symbs, xyzs, angstrom=angstrom)


def translate_along_matrix(geo, disp_mat, angstrom=False):
    """Translate the coordinates of a molecular geometry along
    a matrix.
    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param xyz: vector to translate along
    :type xyz: tuple(float)
    :param idxs: indices of atoms whose coordinates are to be translated
    :type idxs: tuple(int)
    :param angstrom: whether or not the translation is in angstrom
    :type angstrom: bool
    :rtype: automol molecular geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)
    xyzs = numpy.add(xyzs, disp_mat)

    return from_data(symbs, xyzs, angstrom=angstrom)


def perturb(geo, atm_idx, pert_xyz):
    """Perturb the position of one atom by
    changing the value of an xyz coord by some amount.
    """

    # Get the xyz coordinates of the atom to perturb
    atm_coords = list(coordinates(geo)[atm_idx])

    # Get the perturbed set of atomic coordinates
    for idx, val in enumerate(pert_xyz):
        atm_coords[idx] += val
    pert_dct = {atm_idx: atm_coords}

    # Perturb the coordinates of the atom
    pert_geo = set_coordinates(geo, pert_dct)

    return pert_geo


def rotate(geo, axis, angle, orig_xyz=(0.0, 0.0, 0.0), idxs=None, degree=False):
    """Rotate the coordinates of a molecular geometry about
    an axis by a specified angle. A set of `idxs` can be supplied
    to transform a subset of coordinates.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param axis: axis to rotate about
    :type axis: tuple(float)
    :param angle: angle of rotation
    :type angle: float
    :param orig_xyz: xyz coordinates of the origin
    :type orig_xyz: tuple(float)
    :param idxs: indices of atoms whose coordinates are to be rotated
    :type idxs: tuple(int)
    :param degree: is the rotation angle in degrees? If not, radians.
    :type degree: bool
    """
    func = util.vector.rotator(axis, angle, orig_xyz=orig_xyz, degree=degree)
    return transform(geo, func, idxs=idxs)


def transform(geo, func, idxs=None):
    """Transform the coordinates of a geometry by a function.
    A set of `idxs` can be supplied to transform a subset of coordinates.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param func: transformation function
    :type func: function object
    :param idxs: indices representing the subset of atoms
    :type idxs: tuple(int)
    """

    idxs = list(range(count(geo))) if idxs is None else idxs
    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = [func(xyz) if idx in idxs else xyz for idx, xyz in enumerate(xyzs)]

    return from_data(symbs, xyzs)


def transform_by_matrix(geo, mat, trans: bool = True):
    """Transform the coordinates of a molecular geometry by multiplying
    it by some input transformation matrix.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param mat: transformation matrix
    :type mat: tuple(tuple(float))
    :param trans: Whether to transpose the matrix for multiplication as r . M^T
    :rtype: automol moleculer geometry data structure
    """
    if trans:
        mat = numpy.transpose(mat)

    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.dot(xyzs, mat)

    return from_data(symbs, xyzs)


def reflect_coordinates(geo, idxs=None, axes=("x",)):
    """Reflect a specified set of coordinates of a molecular geometry
    about some each of the requested axes.

    A set of `idxs` can be supplied to transform a subset of coordinates.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idxs: indices of atoms whose coordinates are to be reflected
    :type idxs: tuple(int)
    :param axes: axes to reflect about
    :type axes: tuple(str)
    :rtype: automol geometry data structure
    """

    idxs = list(range(count(geo))) if idxs is None else idxs

    # check input
    assert all(idx < len(geo) for idx in idxs)
    assert all(axis in ("x", "y", "z") for axis in axes)

    # get coords
    coords = coordinates(geo)

    # convert x,y,z to nums
    axes = [AXIS_DCT[axis] for axis in axes]

    # build set atom dct with relected coords
    reflect_dct = {}
    for idx in idxs:
        coord_lst = list(coords[idx])
        for axis in axes:
            coord_lst[axis] *= -1.0
        reflect_dct[idx] = coord_lst

    # Reflect coords with dct
    geo_reflected = set_coordinates(geo, reflect_dct)

    return geo_reflected


def shift_atom_position(geo, idx1, idx2):
    """Move the atom at position idx1 to idx2, shifting all other atoms.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param idx1: index of atom 1 in the pair to be measured
    :type idx1: int
    :param idx2: index of atom 2 in the pair to be measured
    :type idx2: int
    :rtype: automol geometry data structure
    """

    # Get the coordinates at idx1 that are to be moved
    geo = [list(x) for x in geo]
    moving_coords = geo[idx1]

    # move the coordinates to idx2
    geo.remove(moving_coords)
    geo.insert(idx2, moving_coords)
    geo_move = tuple(tuple(x) for x in geo)

    return geo_move


def align_to_atoms(geo, idxs: Sequence[int], oop_idxs: tuple[int, int] | None = None):
    """Align geometry to a set of 3 atoms.

    The geometry is aligned to 3 atoms specified by `idxs` as follows:
        Origin: Position of atom 1
        X axis: Unit vector along atom 1 -> atom 2
        Y axis: Unit vector toward atom 1 -> atom 3, perpendicular to X
        Z axis: Unit vector completing a right-handed coordinate system (X cross Y)

    Optionally, one can specify an additional pair of atoms, `oop_idxs`, to
    define the out-of-plane direction, reflecting the geometry if necessary. If
    atom 1 -> atom 2 does not point out-of-plane (along the +Z Axis) for this
    pair, then the geometry will be reflected through the XY plane so that it
    does.

    :param geo: Geometry
    :param idxs: Indices of up to 3 atoms for alignment
    :param oop_idxs: Indices of 2 atoms for reflection
    :return: Geometry
    """
    if not 1 <= len(idxs) <= 4:
        msg = f"Expected 1-4 indices, but got {len(idxs)}: {idxs}"
        raise ValueError(msg)

    # 0. Translate to origin
    origin, *axes = map(numpy.array, coordinates(geo, idxs=idxs))
    geo = translate(geo, -origin)
    axes = [v - origin for v in axes]

    # 1. Define x-axis from atom 1
    x_axis = axes.pop(0) if axes else [1, 0, 0]
    x_axis = util.vector.unit_norm(x_axis)

    # 2. Define y-axis from atom 2
    y_axis = axes.pop(0) if axes else util.vector.arbitrary_unit_perpendicular(x_axis)
    y_axis = util.vector.orthogonalize(x_axis, y_axis, normalize=True)

    # 3. Define z-axis from atom 3
    z_axis = util.vector.unit_perpendicular(x_axis, y_axis)

    # Transform into alignment
    change_basis = numpy.column_stack((x_axis, y_axis, z_axis)).T
    geo = transform_by_matrix(geo, mat=change_basis)

    # If requested, reflect to match out-of-plane vector
    if oop_idxs is not None:
        pos1, pos2 = map(numpy.array, coordinates(geo, idxs=oop_idxs))
        z_comp = numpy.vdot(pos2 - pos1, [0, 0, 1])
        if z_comp < 0:
            geo = reflect_coordinates(geo, axes=["z"])

    return geo
