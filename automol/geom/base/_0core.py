"""
    Core functions defining the geometry data type
"""

import itertools

import more_itertools as mit
import numpy
import pyparsing as pp
from pyparsing import pyparsing_common as ppc
from phydat import phycon, ptab

import automol.formula
from automol import util

AXIS_DCT = {"x": 0, "y": 1, "z": 2}

CHAR = pp.Char(pp.alphas)
SYMBOL = pp.Combine(CHAR + pp.Opt(CHAR))
XYZ_LINE = pp.Group(SYMBOL + pp.Group(ppc.fnumber * 3))
XYZ_LINES = pp.delimitedList(XYZ_LINE, delim=pp.lineEnd())


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


def from_subset(geo, idxs):
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
def symbols(geo, idxs=None):
    """Obtain the atomic symbols atoms in the molecular geometry.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param idxs: indexs of atoms to obtain information for
    :type idxs: tuple(int)
    :rtype: tuple(str)
    """

    idxs = list(range(count(geo))) if idxs is None else idxs

    if geo:
        symbs, _ = zip(*geo)
    else:
        symbs = ()

    symbs = tuple(symb for idx, symb in enumerate(symbs) if idx in idxs)
    return symbs


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
    :rtype: automol molecular graph data structure
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
def string(geo, angstrom=True):
    """Write a molecular geometry to a string:
       symb1  xyz1 xyz2 xyz3
       symbn  xyzn xyzn xyzn

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param angstrom: parameter to control coordinate conversion to Angstrom
    :type angstrom: bool
    :rtype: str
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)

    natms = len(symbs)
    assert len(xyzs) == natms

    geo_str = "\n".join(
        f"{symb:2s} {xyz[0]:10.6f} {xyz[1]:10.6f} {xyz[2]:10.6f}"
        for symb, xyz in zip(symbs, xyzs)
    )

    return geo_str


def xyz_string(geo, comment=""):
    """Write a molecular geometry to a string:
       natom
       comment
       symb1  xyz1 xyz2 xyz3
       symbn  xyzn xyzn xyzn

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param comment: string to place in the comment line of string
    :type comment: str
    :rtype: str
    """
    geo_str = string(geo, angstrom=True)
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
    rows = XYZ_LINES.parseString(geo_str).asList()

    symbs, xyzs = zip(*rows) if rows else ([], [])

    geo = from_data(symbs, xyzs, angstrom=angstrom)

    return geo


def from_xyz_string(xyz_str):
    """Read a Cartesian molecular geometry from a string that matches the
    format of a string of a standard .xyz file.

    :param xyz_str: string obtained from reading the .xyz file
    :type xyz_str: str
    :rtype: automol geometry data structure
    """
    lines = xyz_str.splitlines()
    natms = int(lines.pop(0).strip())
    lines.pop(0)

    geo_str = "\n".join(lines)
    geo = from_string(geo_str, angstrom=True)

    assert natms == count(geo), f"XYZ string with inconsistent count: {xyz_str}"

    return geo


def xyz_string_comment(xyz_str):
    """Read the comment line of a string of a standard .xyz file.

    :param xyz_str: string obtained from reading the .xyz file
    :type xyz_str: str
    :rtype: str
    """
    return xyz_str.splitlines()[1].strip()


def from_xyz_trajectory_string(geo_str):
    """Read a series of molecular geometries from a trajectory string.

    :param geo_str: string containing the geometry
    :type geo_str: str
    :rtype: (tuple(automol geometry data structure), tuple(str))
    """

    def _blocks(lst, size):
        """Split list into parts of size n"""
        split_lst = []
        for i in range(0, len(lst), size):
            split_lst.append(lst[i: i + size])
        return split_lst

    # Split the lines for iteration
    geo_lines = [line for line in geo_str.splitlines() if line != ""]

    # Get the number of atoms used to partition the trajectory file
    line1 = geo_lines[0].strip()
    natoms = int(line1)

    geoms, comments = tuple(), tuple()
    for block in _blocks(geo_lines, natoms + 2):
        comments += (block[1],)
        geoms += (from_string("\n".join(block[2:])),)

    return tuple(zip(geoms, comments))


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
    fml_str = automol.formula.string(fml)

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


def is_linear(geo, tol=2.0 * phycon.DEG2RAD):
    """Determine if the molecular geometry is linear.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param tol: tolerance of bond angle(s) for determing linearity
    :type tol: float
    :rtype: bool
    """

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


def principal_axes(geo, amu=True):
    """Determine the principal axes of rotation for a molecular geometry.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(tuple(float))
    """

    ine = inertia_tensor(geo, amu=amu)
    _, paxs = numpy.linalg.eigh(ine)
    paxs = tuple(map(tuple, paxs))

    return paxs


def moments_of_inertia(geo, amu=True):
    """Calculate the moments of inertia along the xyz axes
    (these not sorted in to A,B,C).

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(tuple(float))
    """

    ine = inertia_tensor(geo, amu=amu)
    moms, _ = numpy.linalg.eigh(ine)
    moms = tuple(moms)

    return moms


def rotational_constants(geo, amu=True):
    """Calculate the rotational constants.
    (these not sorted in to A,B,C).

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param amu: parameter to control electron mass -> amu conversion
    :type amu: bool
    :rtype: tuple(float)
    """

    moms = moments_of_inertia(geo, amu=amu)
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
    dist = util.vec.distance(xyz1, xyz2)
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
    ang = util.vec.central_angle(xyz1, xyz2, xyz3)
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
    dih = util.vec.dihedral_angle(xyz1, xyz2, xyz3, xyz4)
    dih *= phycon.RAD2DEG if degree else 1

    return dih


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

    :param geo1: molecular geometry 1
    :type geo1: automol molecular geometry data structure
    :param geo2: molecular geometry 2
    :type geo2: automol molecular geometry data structure
    :param dist_cutoff: threshhold for center-of-mass distance
    :type: dist_cutoff: float
    :param theta: theta angle for intermolecular orientation
    :type theta: float
    :param phi: phi angle for intermolecular orientation
    :type phi: float
    :rtype: automol molecular geometry data structure
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
        util.vec.distance(xyz1, xyz2) for xyz1, xyz2 in itertools.product(xyzs1, xyzs2)
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
                if util.vec.distance(xyz, ref_xyz) < thresh
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

    return from_subset(geo, non_dummy_idxs)


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


def rotate(geo, axis, angle, orig_xyz=None, idxs=None, degree=False):
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
    angle = angle if not degree else angle * phycon.DEG2RAD

    func = util.vec.rotator(axis, angle, orig_xyz=orig_xyz)

    return transform(geo, func, idxs=idxs)


def euler_rotate(geo, theta, phi, psi):
    """Rotate the coordinates of a molecular geometry about
    the three Euler angles.

    :param geo: molecular geometry
    :type geo: automol geometry data structure
    :param theta: angle to rotate about z-axis
    :type theta: float
    :param phi: angle to rotate about x'-axis
    :type phi: float
    :param psi: angle to rotate about z'-axis
    :type psi: float
    :rtype: automol geometry data structure
    """

    mat = util.mat.euler_rotation_matrix(theta, phi, psi)

    return transform_by_matrix(geo, mat)


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


def transform_by_matrix(geo, mat):
    """Transform the coordinates of a molecular geometry by multiplying
    it by some input transfomration matrix.

    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :param mat: transformation matrix
    :type mat: tuple(tuple(float))
    :rtype: automol moleculer geometry data structure
    """

    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.dot(xyzs, numpy.transpose(mat))

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
