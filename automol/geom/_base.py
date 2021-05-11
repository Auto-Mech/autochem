""" cartesian geometries
"""

import itertools
from phydat import phycon, ptab
import autoread as ar
import autowrite as aw
from automol.graph._graph_dep import atom_keys
from automol.graph._graph_dep import bond_keys
from automol.graph.geom import symbols as _symbols
from automol.graph.geom import coordinates as _coordinates
from automol.graph.geom import count as _count
from automol.convert.geom import from_subset as _from_subset
from automol.convert.geom import connected as _connected
from automol.convert.geom import components_graph as _components_graph
from automol.convert.geom import connectivity_graph as _connectivity_graph
from automol.convert.geom import without_dummy_atoms as _without_dummy_atoms
from automol.convert.geom import central_angle as _central_angle
from automol.convert.geom import distance as _distance
from automol.convert.geom import dihedral_angle as _dihedral_angle
from automol.convert.geom import linear_atoms as _linear_atoms
import automol.create as _create
from automol import util


# constructors
def from_subset(geo, idxs):
    """ Generate a new molecular geometry from a subset of the atoms in an
        input geometry.

        (Rename this and put it under operations?)

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indices representing the subset of atoms
        :type idxs: tuple(int)
        :rtype: automol moleculer geometry data structure
    """
    return _from_subset(geo, idxs)


# getters
def symbols(geo, idxs=None):
    """ Obtain the atomic symbols atoms in the molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indexs of atoms to obtain information for
        :type idxs: tuple(int)
        :rtype: tuple(str)
    """
    return _symbols(geo, idxs=idxs)


def coordinates(geo, idxs=None, angstrom=False):
    """ Obtain the Cartesian coordinates of atoms in the molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indexs of atoms to obtain information for
        :type idxs: tuple(int)
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :rtype: tuple(tuple(float))
    """
    return _coordinates(geo, idxs=idxs, angstrom=angstrom)


def count(geo):
    """ Obtain the number of rows of the molecular geometry, which corresponds to
        the number of atoms in the geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: int
    """
    return _count(geo)


def atom_count(geo, symb, match=True):
    """ Count the number of atoms of a particular type in the geometry.

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
    """ Obtain the indices of a atoms of a particular type in the geometry.

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
    """ Obtain the indices of a atoms of a particular type in the geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: tuple(int)
    """

    symbs = symbols(geo)
    dummy_idxs = [idx for idx, symb in enumerate(symbs)
                  if not ptab.to_number(symb)]

    return tuple(dummy_idxs)


def components_graph(geo, stereo=True):
    """ Generate a list of molecular graphs where each element is a graph that
        consists of fully connected (bonded) atoms. Stereochemistry is included
        if requested.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: automol molecular graph data structure
    """
    return _components_graph(geo, stereo=stereo)


def connected(geo, stereo=True):
    """ Determine if all atoms in geometry are completely connected.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param stereo: parameter to include stereochemistry information
        :type stereo: bool
        :rtype: bool
    """
    return _connected(geo, stereo=stereo)


# validation
def is_valid(geo):
    """ Assess if a molecular geometry has the proper data structure.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: bool
    """

    ret = hasattr(geo, '__iter__')
    if ret:
        ret = all(hasattr(obj, '__len__') and len(obj) == 2 for obj in geo)
        if ret:
            symbs, xyzs = zip(*geo)
            try:
                _create.geom.from_data(symbs, xyzs)
            except AssertionError:
                ret = False

    return ret


# setters
def set_coordinates(geo, xyz_dct):
    """ Set coordinate values for the molecular geometry,
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

    xyzs = [xyz_dct[idx] if idx in xyz_dct else xyz
            for idx, xyz in enumerate(xyzs)]

    return _create.geom.from_data(symbs, xyzs)


def without_dummy_atoms(geo):
    """ Return a copy of the molecular geometry without dummy atoms.
    """

    return _without_dummy_atoms(geo)


# I/O
def from_string(geo_str, angstrom=True):
    """ Read a Cartesian molecular geometry from a string comprised
        of just the atomic symbols and coordinates.

        :param geo_str: string containing the geometry
        :type geo_str: str
        :param angstrom: parameter to control coordinate conversion to Angstrom
        :type angstrom: bool
        :rtype: automol geometry data structure
    """

    symbs, xyzs = ar.geom.read(geo_str)
    geo = _create.geom.from_data(symbs, xyzs, angstrom=angstrom)

    return geo


def from_xyz_string(xyz_str):
    """ Read a Cartesian molecular geometry from a string that matches the
        format of a string of a standard .xyz file.

        :param xyz_str: string obtained from reading the .xyz file
        :type xyz_str: str
        :rtype: automol geometry data structure
    """

    symbs, xyzs = ar.geom.read_xyz(xyz_str)
    geo = _create.geom.from_data(symbs, xyzs, angstrom=True)

    return geo


def xyz_string_comment(xyz_str):
    """ Read the comment line of a string of a standard .xyz file.

        :param xyz_str: string obtained from reading the .xyz file
        :type xyz_str: str
        :rtype: str
    """
    return xyz_str.splitlines()[1].strip()


def string(geo, angstrom=True):
    """ Write a molecular geometry to a string:
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
    geo_str = aw.geom.write(symbs=symbs, xyzs=xyzs)

    return geo_str


def xyz_string(geo, comment=None):
    """ Write a molecular geometry to a string:
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

    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=True)
    geo_str = aw.geom.write_xyz(symbs=symbs, xyzs=xyzs, comment=comment)

    return geo_str


def xyz_trajectory_string(geo_lst, comments=None):
    """ Write a series of molecular geometries to trajectory file which
        is a string that collated by several xyz-file format geometry strings.

        :param geo_lst: list of molecular geometries
        :type geo_lst: tuple(automol geometry data structure)
        :param comments: list of comments for each of the molecular geometries
        :type comments: tuple(str)
        :rtype: str
    """

    symbs_lst = [symbols(geo) for geo in geo_lst]
    xyzs_lst = [coordinates(geo, angstrom=True) for geo in geo_lst]
    assert len(set(symbs_lst)) == 1, (
        'set {} \n symbol lists {}'.format(
            set(symbs_lst), symbs_lst))
    symbs = symbs_lst[0]
    xyz_traj_str = aw.geom.write_xyz_trajectory(symbs, xyzs_lst,
                                                comments=comments)
    return xyz_traj_str


def from_xyz_trajectory_string(geo_str):
    """ Read a series of molecular geometries from a trajectory string.

        :param geo_str: string containing the geometry
        :type geo_str: str
        :rtype: (tuple(automol geometry data structure), tuple(str))
    """

    def _blocks(lst, size):
        """ Split list into parts of size n"""
        split_lst = []
        for i in range(0, len(lst), size):
            split_lst.append(lst[i:i + size])
        return split_lst

    # Split the lines for iteration
    geo_lines = [line for line in geo_str.splitlines()
                 if line != '']

    # Get the number of atoms used to partition the trajectory file
    line1 = geo_lines[0].strip()
    natoms = int(line1)

    geoms, comments = tuple(), tuple()
    for block in _blocks(geo_lines, natoms+2):
        comments += (block[1],)
        geoms += (from_string('\n'.join(block[2:])),)

    return tuple(zip(geoms, comments))


# geometric properties
def distance(geo, idx1, idx2, angstrom=False):
    """ Measure the distance between two atoms in a molecular geometry.

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

    return _distance(
        geo, idx1, idx2, angstrom=angstrom)


def central_angle(geo, idx1, idx2, idx3, degree=False):
    """ Measure the angle inscribed by three atoms in a molecular geometry.

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

    return _central_angle(
        geo, idx1, idx2, idx3, degree=degree)


def dihedral_angle(geo, idx1, idx2, idx3, idx4, degree=False):
    """ Measure the angle inscribed by three atoms in a molecular geometry.

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

    return _dihedral_angle(
        geo, idx1, idx2, idx3, idx4, degree=degree)


def zmatrix_row_values(geo, idx, idx1=None, idx2=None, idx3=None,
                       angstrom=True, degree=True):
    """ Get coordinate values for a row in a Z-Matrix.

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


def linear_atoms(geo, gra=None, tol=5.):
    """ find linear atoms in a geometry (atoms with 180 degree bond angle)

        :param geo: the geometry
        :type geo: automol geometry data structure
        :param gra: the graph describing connectivity; if None, a connectivity
            graph will be generated using default distance thresholds
        :type gra: automol graph data structure
        :param tol: the tolerance threshold for linearity, in degrees
        :type tol: float
        :rtype: tuple(int)
    """
    return _linear_atoms(geo, gra=gra, tol=tol)


def closest_unbonded_atoms(geo, gra=None):
    """ Determine which pair of unbonded atoms in a molecular geometry
        are closest together.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param gra: the graph describing connectivity; if None, a connectivity
            graph will be generated using default distance thresholds
        :type gra: automol graph data structure
        :rtype: (frozenset(int), float)
    """

    gra = _connectivity_graph(geo) if gra is None else gra
    atm_keys = atom_keys(gra)
    bnd_keys = bond_keys(gra)
    poss_bnd_keys = set(map(frozenset, itertools.combinations(atm_keys, r=2)))

    # The set of candidates includes all unbonded pairs of atoms
    cand_bnd_keys = poss_bnd_keys - bnd_keys

    min_bnd_key = None
    min_dist_val = 1000.
    for bnd_key in cand_bnd_keys:
        dist_val = distance(geo, *bnd_key)
        if dist_val < min_dist_val:
            min_dist_val = dist_val
            min_bnd_key = bnd_key

    return min_bnd_key, min_dist_val
