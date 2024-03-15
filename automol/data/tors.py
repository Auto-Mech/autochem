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
DihCoord = Tuple[int, int, int, int]
Groups = Tuple[List[int], List[int]]
Grid = List[float]


@dataclasses.dataclass
class Torsion:
    """Encodes information for a single torsion, which is one component of a rotor

    :param name: The z-matrix coordinate name
    :type name: str
    :param coordinate: The z-matrix keys defining the torsion coordinate
    :type coordinate: DihKey
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symmetry: The rotational symmetry number of the torsion
    :type symmetry: int
    """

    name: str
    coordinate: DihCoord
    groups: Groups
    symmetry: int


# Torsion functions
# # Constructors
def from_data(name_: str, coo: DihCoord, grps: Groups, symm: int) -> Torsion:
    """Construct a torsion from data

    :param name_: The z-matrix coordinate name
    :type name_: str
    :param coo: The z-matrix keys defining the torsion coordinate
    :type coo: DihKey
    :param groups: The sets of atoms keys defining the rotational groups
    :type groups: Tuple[List[int], List[int]]
    :type symm: The rotational symmetry number of the torsion
    :type symm: int
    :return: The torsion data structure
    :rtype: Torsion
    """
    assert len(coo) == 4, f"Invalid torsion coordinate: {coo}"
    assert len(grps) == 2, f"Invalid torsion groups: {grps}"
    return Torsion(
        name=str(name_),
        coordinate=tuple(coo),
        groups=tuple(map(tuple, grps)),
        symmetry=int(symm),
    )


# # Getters
def name(tor: Torsion) -> str:
    """Get the coordinate name of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :return: The coordinate name
    :rtype: str
    """
    return tor.name


def coordinate(
    tor: Torsion, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None
) -> DihCoord:
    """Get the z-matrix keys defining the torsion coordinate

    :param tor: A torsion
    :type tor: Torsion
    :param key_typ: The type of keys to return, "zmat" (default) or "geom"
    :type key_typ: str, optional
    :param zc_: Z-matrix conversion info, to avoid re-calculation, defaults to None
    :type zc_: Optional[ZmatConv], optional
    :return: The torsion rotational keys
    :rtype: str
    """
    coo = tor.coordinate
    if key_typ == "geom":
        coo = zmat_conv.relabel_zmatrix_key_sequence(zc_, coo, dummy=True)
    return coo


def axis(tor: Torsion, key_typ: str = "zmat", zc_: Optional[ZmatConv] = None) -> Axis:
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
    ks_ = coordinate(tor, key_typ=key_typ, zc_=zc_)
    return ks_[1:3]


def groups(
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
    grps = tor.groups
    if key_typ == "geom":
        grps = zmat_conv.relabel_zmatrix_key_sequence(zc_, grps)
    return grps


def symmetry(tor: Torsion) -> int:
    """Get the rotational symmetry of a torsion

    :param tor: A torsion
    :type tor: Torsion
    :return: The rotational symmetry
    :rtype: int
    """
    return tor.symmetry


def span(tor: Torsion) -> float:
    """Get the angular span of a torsion, based on the symmetry number

    :param tor: A torsion
    :type tor: Torsion
    :return: The angular span of the torsion, 2 pi / symmetry
    :rtype: float
    """
    return 2 * numpy.pi / symmetry(tor)


def grid(tor: Torsion, zma, increment: float=30 * phycon.DEG2RAD) -> Grid:
    """Get the coordinate grid for a torsion

    :param tor: A torsion
    :type tor: Torsion
    :param zma: The z-matrix associated with this torsion
    :type zma: automol zmat data structure
    :param increment: The grid increment, in radians
    :type increment: float
    :return: The coordinate grid
    :rtype: Grid
    """
    # [0, 30, 60, 90, ...] << in degrees
    vals = numpy.arange(0, span(tor), increment)
    # Start from the equilibrium value
    vals += zmat.value(zma, name(tor))
    return tuple(map(float, vals))


# # Transformations
def with_geometry_indices(tor: Torsion, zc_: ZmatConv) -> Torsion:
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
    return from_data(
        name_=name(tor),
        coo=coordinate(tor, key_typ="geom", zc_=zc_),
        grps=groups(tor, key_typ="geom", zc_=zc_),
        symm=symmetry(tor),
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
        coo = [k + shift for k in coordinate(tor)]
        axis1, axis2 = (k + shift for k in axis(tor))
        grps = ([k + shift for k in g] for g in groups(tor))
        grp1, grp2 = ("-".join(map(str, g)) if len(g) > 1 else g[0] for g in grps)
        coo_str = "-".join("*" if k is None else str(k) for k in coo)
        tor_yml_dct[name(tor)] = {
            "axis1": axis1,
            "group1": grp1,
            "axis2": axis2,
            "group2": grp2,
            "symmetry": symmetry(tor),
            "coordinate": coo_str,
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
    for name_, vals_dct in tor_yml_dct.items():
        ax_ = list(map(vals_dct.__getitem__, ["axis1", "axis2"]))
        ax_ = [k + shift for k in ax_]
        grps = list(map(vals_dct.__getitem__, ["group1", "group2"]))
        grps = [[g] if isinstance(g, int) else map(int, g.split("-")) for g in grps]
        grps = [[k + shift for k in g] for g in grps]
        coo = vals_dct.get("coordinate", None)
        if coo is None:
            coo = [None, *ax_, None]
        else:
            coo = [None if s == "*" else int(s) + shift for s in coo.split("-")]

        tor = from_data(
            name_=name_,
            coo=coo,
            grps=grps,
            symm=vals_dct["symmetry"],
        )

        tor_lst.append(tor)

    return tuple(tor_lst)
