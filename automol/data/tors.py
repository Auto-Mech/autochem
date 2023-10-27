"""Implements a data structure for encoding information about a single torsion

Torsions are used to construct Rotor data structures
"""
import dataclasses
from typing import List, Optional, Tuple

import numpy
import yaml
from phydat import phycon

from automol import zmat
from automol.util import ZmatConv, zmat_conv

Axis = Tuple[int, int]
Groups = Tuple[List[int], List[int]]
Grid = List[float]


@dataclasses.dataclass
class Torsion:
    """Encodes information for a single torsion, which is one component of a rotor

    :param name: The z-matrix coordinate name
    :type name: str
    :param axis: The pair of atom keys defining the rotational axis
    :type axis: Tuple[int, int]
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symmetry: The rotational symmetry number of the torsion
    :type symmetry: int
    """

    name: str
    axis: Axis
    groups: Groups
    symmetry: int


# Torsion functions
# # Constructors
def torsion_from_data(name, axis, groups, symm) -> Torsion:
    """Construct a torsion from data

    :param name: The z-matrix coordinate name
    :type name: str
    :param axis: The pair of atom keys defining the rotational axis
    :type axis: Tuple[int, int]
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symm: The rotational symmetry number of the torsion
    :type symm: int
    :return: The torsion data structure
    :rtype: Torsion
    """
    assert len(axis) == 2
    assert len(groups) == 2
    return Torsion(
        name=str(name),
        axis=tuple(axis),
        groups=tuple(map(tuple, groups)),
        symmetry=int(symm),
    )


# # Getters
def torsion_name(tor: Torsion) -> str:
    """Get the coordinate name of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :return: The coordinate name
    :rtype: str
    """
    return tor.name


def torsion_axis(
    tor: Torsion, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> Axis:
    """Get the rotational axis of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsion rotational axis
    :rtype: str
    """
    axis = tor.axis
    if key_typ == "geom":
        axis = zmat_conv.geometry_keys(zc_, axis)
    return axis


def torsion_groups(
    tor: Torsion, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> Groups:
    """Get the rotational groups of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsion rotational groups
    :rtype: str
    """
    groups = tor.groups
    if key_typ == "geom":
        groups = zmat_conv.geometry_keys(zc_, groups)
    return groups


def torsion_symmetry(tor: Torsion) -> int:
    """Get the rotational symmetry of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :return: The rotational symmetry
    :rtype: int
    """
    return tor.symmetry


def torsion_grid(
    tor: Torsion, zma, span=2 * numpy.pi, increment=30 * phycon.DEG2RAD
) -> Grid:
    """Get the coordinate grid for a torsion

    :param tor: A torsion
    :type tor: Torsion
    :param zma: The z-matrix associated with this torsion
    :type zma: automol zmat data structure
    :param span: The angular span of the grid, in radians
    :type span: float
    :param increment: The grid increment, in radians
    :type increment: float
    :return: The coordinate grid
    :rtype: Grid
    """
    symm = tor.symmetry
    # [0, 30, 60, 90, ...] << in degrees
    grid = numpy.arange(0, span / symm, increment)
    # Start from the equilibrium value
    grid += zmat.value(zma, tor.name)
    return tuple(map(float, grid))


# # Transformations
def torsion_with_geometry_indices(tor: Torsion, zc_: ZmatConv) -> Torsion:
    """Given a z-matrix torsion and a z-matrix conversion, return a torsion with
    geometry indices

    That is, the axis and group indices will skip dummy atoms

    :param tor: A torsion data structure
    :type tor: Torsion
    :param zc_: A z-matrix conversion data structure
    :type zc_: ZmatConv
    :return: A torsion data structure using geometry indices
    :rtype: Torsion
    """
    return torsion_from_data(
        name=torsion_name(tor),
        axis=torsion_axis(tor, key_typ="geom", zc_=zc_),
        groups=torsion_groups(tor, key_typ="geom", zc_=zc_),
        symm=torsion_symmetry(tor),
    )


# Torsion List functions
def torsions_string(tor_lst: List[Torsion], one_indexed: bool = True) -> str:
    """Write a list of torsions to a string

    :param tor_lst: A list of torsions
    :type tor_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    tor_yml_dct = torsions_yaml_data(tor_lst, one_indexed=one_indexed)
    tor_str = yaml.dump(tor_yml_dct, sort_keys=False)
    return tor_str


def torsions_from_string(tor_str: str, one_indexed: bool = True) -> List[Torsion]:
    """Write a list of torsions to a string

    :param tor_str: A string representations of the torsions, as a flattened list
    :type tor_str: str
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A list of torsions
    :rtype: List[Torsion]
    """
    tor_yml_dct = yaml.load(tor_str, Loader=yaml.FullLoader)
    tor_lst = torsions_from_yaml_data(tor_yml_dct, one_indexed=one_indexed)
    return tor_lst


def torsions_yaml_data(tor_lst: List[Torsion], one_indexed: bool = True) -> dict:
    """Write a list of torsions to a yaml-formatted dictionary

    :param tor_lst: A list of torsions
    :type tor_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    shift = 1 if one_indexed else 0

    tor_yml_dct = {}
    for tor in tor_lst:
        axis1, axis2 = (k + shift for k in tor.axis)
        groups = ([k + shift for k in g] for g in tor.groups)
        group1, group2 = ("-".join(map(str, g)) if len(g) > 1 else g[0] for g in groups)
        tor_yml_dct[tor.name] = {
            "axis1": axis1,
            "group1": group1,
            "axis2": axis2,
            "group2": group2,
            "symmetry": tor.symmetry,
        }
    return tor_yml_dct


def torsions_from_yaml_data(tor_yml_dct: dict, one_indexed: bool = True) -> dict:
    """Read a list of torsions out of a yaml-formatted torsion dictionary

    :param tor_lst: A list of torsions
    :type tor_lst: List[Torsion]
    :param one_indexed: Is this a one-indexed string? defaults to True
    :type one_indexed: bool, optional
    :returns: A string representations of the torsions, as a flattened list
    :rtype: str
    """
    shift = -1 if one_indexed else 0

    tor_lst = []
    for name, vals_dct in tor_yml_dct.items():
        raw_axis = list(map(vals_dct.__getitem__, ["axis1", "axis2"]))
        raw_groups = list(map(vals_dct.__getitem__, ["group1", "group2"]))
        raw_groups = [
            [g] if isinstance(g, int) else map(int, g.split("-")) for g in raw_groups
        ]

        tor = torsion_from_data(
            name=name,
            axis=[k + shift for k in raw_axis],
            groups=[[k + shift for k in g] for g in raw_groups],
            symm=vals_dct["symmetry"],
        )

        tor_lst.append(tor)

    return tuple(tor_lst)
